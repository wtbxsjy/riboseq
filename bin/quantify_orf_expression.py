#!/usr/bin/env python3
"""
从 RiboseQC P-site bedgraph 为每个 unified ORF 提取 per-sample 表达量 (reads + pN)。

不依赖 gencode 或 AMP 特定格式 —— 接受标准的 unified ORF metadata 作为输入。

用法:
  quantify_orf_expression.py \
      --orf-meta unified_orfs.metadata.tsv \
      --psites-dir riboseqc/ \
      --sample-pattern "*_P_sites_plus.bedgraph" \
      --output expression_summary.tsv \
      --min-ocs 0.0 \
      --max-orfs 100 \
      --workers 4
"""

import argparse
import os
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd


class ORFExpressionQuantifier:
    """从 RiboseQC P-site bedgraph 查询每个 unified ORF 的定量表达值。"""

    def __init__(
        self,
        orf_meta_path,
        psites_dir,
        sample_pattern="*_P_sites_plus.bedgraph",
        orf_confidence_path=None,
        min_ocs=0.0,
        max_orfs=0,
        workers=4,
        chrom_prefix="",
    ):
        self.orf_meta_path = Path(orf_meta_path)
        self.psites_dir = Path(psites_dir)
        self.sample_pattern = sample_pattern
        self.orf_confidence_path = (
            Path(orf_confidence_path) if orf_confidence_path else None
        )
        self.min_ocs = min_ocs
        self.max_orfs = max_orfs
        self.workers = workers
        self.chrom_prefix = chrom_prefix

        self.orf_df = None
        self.samples = []
        self.confidence_df = None

    # ── Loading ──────────────────────────────────────────────────────────

    def load_orf_metadata(self):
        """加载 unified ORF metadata 并提取坐标列。"""
        print(f"加载 ORF metadata: {self.orf_meta_path}")
        self.orf_df = pd.read_csv(self.orf_meta_path, sep="\t")
        print(f"  加载了 {self.orf_df.shape[0]} 个 ORF, {self.orf_df.shape[1]} 列")

        # 验证必需列
        required = ["orf_id", "chrom", "start", "end", "strand"]
        missing = [c for c in required if c not in self.orf_df.columns]
        if missing:
            raise ValueError(f"Metadata 缺少必需列: {missing}")

        # 标准化类型
        self.orf_df["chrom"] = self.orf_df["chrom"].astype(str)
        self.orf_df["start"] = self.orf_df["start"].astype(int)
        self.orf_df["end"] = self.orf_df["end"].astype(int)
        self.orf_df["strand"] = self.orf_df["strand"].astype(str)

        # 应用染色体前缀
        if self.chrom_prefix:
            self.orf_df["chrom"] = self.chrom_prefix + self.orf_df["chrom"]

        # 限制 ORF 数量 (即使没有 confidence 文件)
        if self.max_orfs > 0 and self.orf_df.shape[0] > self.max_orfs:
            self.orf_df = self.orf_df.head(self.max_orfs)
            print(f"  限制为前 {self.max_orfs} 个 ORF (无 OCS 排序)")

    def load_confidence(self):
        """加载 ORF confidence 数据并过滤。"""
        if self.orf_confidence_path is None or not self.orf_confidence_path.exists():
            print("未提供 OCS 文件, 跳过置信度过滤")
            return
            return

        print(f"加载 ORF confidence: {self.orf_confidence_path}")
        self.confidence_df = pd.read_csv(self.orf_confidence_path, sep="\t")

        if "ocs" not in self.confidence_df.columns:
            print("  WARNING: OCS 列不存在, 跳过过滤")
            self.confidence_df = None
            return

        # 合并 OCS 到 orf_df
        self.orf_df = self.orf_df.merge(
            self.confidence_df[["orf_id", "ocs", "tier"]],
            on="orf_id",
            how="left",
        )

        # 过滤低 OCS
        before = self.orf_df.shape[0]
        self.orf_df = self.orf_df[self.orf_df["ocs"] >= self.min_ocs]
        print(f"  OCS >= {self.min_ocs}: {before} → {self.orf_df.shape[0]} ORFs")

        # 限制数量
        if self.max_orfs > 0 and self.orf_df.shape[0] > self.max_orfs:
            # 按 OCS 降序取 top N
            self.orf_df = self.orf_df.sort_values("ocs", ascending=False).head(
                self.max_orfs
            )
            print(f"  限制为 top {self.max_orfs} ORFs")

    def discover_samples(self):
        """扫描 bedgraph 目录发现样本名称。"""
        pattern = self.sample_pattern
        files = sorted(self.psites_dir.glob(pattern))
        if not files:
            # 尝试不带 _plus 后缀
            files = sorted(self.psites_dir.glob("*_P_sites_plus.bedgraph"))

        for f in files:
            name = f.name
            # 从 "Ribo_11_P_sites_plus.bedgraph" → "Ribo_11"
            for suffix in ["_P_sites_plus.bedgraph", "_P_sites_minus.bedgraph"]:
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
                    break
            if name and name not in self.samples:
                self.samples.append(name)

        print(f"发现 {len(self.samples)} 个样本: {', '.join(self.samples[:10])}...")

    # ── Bedgraph query ───────────────────────────────────────────────────

    @staticmethod
    def query_bedgraph_region(bedgraph_file, chrom, start, end):
        """用 awk 快速查询 bedgraph 中指定区间的行。"""
        if not Path(bedgraph_file).exists():
            return None

        try:
            cmd = [
                "awk",
                "-F\\t",
                f'$1=="{chrom}" && $2<{end} && $3>{start} {{print $0}}',
                str(bedgraph_file),
            ]
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0 or not result.stdout.strip():
                return None

            lines = result.stdout.strip().split("\n")
            data = []
            for line in lines:
                if not line:
                    continue
                parts = line.split("\t")
                try:
                    data.append(
                        [parts[0], int(parts[1]), int(parts[2]), float(parts[3])]
                    )
                except (ValueError, IndexError):
                    pass

            if data:
                return pd.DataFrame(
                    data, columns=["chrom", "start", "end", "value"]
                )
            return None
        except subprocess.TimeoutExpired:
            return None
        except Exception:
            return None

    @staticmethod
    def calculate_pn(values):
        """计算 pN 值: max(coverage) / mean(coverage)。"""
        if len(values) == 0 or np.sum(values) == 0:
            return 0.0
        pn = np.max(values) / (np.mean(values) + 1e-10)
        return round(float(pn), 4)

    def query_orf_single_sample(self, orf_row, sample):
        """查询单个 ORF 在单个样本中的 P-site 数据。"""
        chrom = orf_row["chrom"]
        start = orf_row["start"]
        end = orf_row["end"]
        strand = orf_row["strand"]

        # 选择 bedgraph 文件
        if strand == "+":
            bg_file = self.psites_dir / f"{sample}_P_sites_plus.bedgraph"
        else:
            bg_file = self.psites_dir / f"{sample}_P_sites_minus.bedgraph"

        if not bg_file.exists():
            return {"reads": 0, "mean_cov": 0.0, "max_cov": 0, "pN": 0.0}

        df = self.query_bedgraph_region(bg_file, chrom, start, end)
        if df is None or len(df) == 0:
            return {"reads": 0, "mean_cov": 0.0, "max_cov": 0, "pN": 0.0}

        values = df["value"].values
        total_reads = int(np.sum(values))

        # 对于 P-site bedgraph, value 是 P-site count
        # 需要计算: total P-sites in ORF region
        # bedgraph 记录的是 "每个区间的 P-site 数量"
        # 区间长度 × value 的总和
        total_psites = int(
            np.sum((df["end"].values - df["start"].values) * values)
        )

        return {
            "reads": total_psites,
            "mean_cov": round(float(np.mean(values)), 2),
            "max_cov": int(np.max(values)),
            "pN": self.calculate_pn(values),
        }

    # ── Main processing ──────────────────────────────────────────────────

    def process_single_orf(self, orf_row):
        """处理单个 ORF: 查询所有样本并返回一行结果。"""
        result = {
            "orf_id": orf_row["orf_id"],
            "chrom": orf_row["chrom"],
            "start": orf_row["start"],
            "end": orf_row["end"],
            "strand": orf_row["strand"],
        }
        if "gene_id" in orf_row:
            result["gene_id"] = orf_row["gene_id"]
        if "ocs" in orf_row:
            result["ocs"] = orf_row.get("ocs", np.nan)
        if "tier" in orf_row:
            result["tier"] = orf_row.get("tier", "")

        total_reads = 0
        n_expressed = 0
        for sample in self.samples:
            q = self.query_orf_single_sample(orf_row, sample)
            result[f"{sample}_reads"] = q["reads"]
            result[f"{sample}_pN"] = q["pN"]
            total_reads += q["reads"]
            if q["reads"] > 0:
                n_expressed += 1

        result["total_reads"] = total_reads
        result["n_expressed_samples"] = n_expressed
        return result

    def run(self):
        """执行完整的表达量提取流程。"""
        t0 = time.time()

        # 1. 加载
        self.load_orf_metadata()
        self.load_confidence()
        self.discover_samples()

        if not self.samples:
            print("ERROR: 未发现任何样本 bedgraph 文件")
            return None

        if self.orf_df.shape[0] == 0:
            print("ERROR: 没有 ORF 需要处理 (过滤太严格?)")
            return None

        # 2. 处理
        n_orfs = self.orf_df.shape[0]
        print(f"\n开始提取 {n_orfs} 个 ORF 的 P-site 表达量...")
        print(f"  使用 {self.workers} 个 worker, {len(self.samples)} 个样本")

        results = []
        rows = [row for _, row in self.orf_df.iterrows()]

        with ProcessPoolExecutor(max_workers=self.workers) as exe:
            futures = {exe.submit(self.process_single_orf, row): i for i, row in enumerate(rows)}
            for fut in as_completed(futures):
                try:
                    res = fut.result()
                    results.append(res)
                except Exception as e:
                    idx = futures[fut]
                    print(f"  ERROR ORF idx={idx}: {e}")

        # 3. 保存
        result_df = pd.DataFrame(results)

        # 确保列顺序一致
        static_cols = [
            "orf_id", "chrom", "start", "end", "strand",
        ]
        extra_cols = []
        for c in ["gene_id", "ocs", "tier"]:
            if c in result_df.columns:
                extra_cols.append(c)
        sample_cols = []
        for s in self.samples:
            sample_cols.append(f"{s}_reads")
            sample_cols.append(f"{s}_pN")
        tail_cols = ["total_reads", "n_expressed_samples"]

        ordered_cols = (
            [c for c in static_cols + extra_cols if c in result_df.columns]
            + sample_cols
            + tail_cols
        )
        result_df = result_df[ordered_cols]

        elapsed = time.time() - t0
        print(f"\n完成! 提取了 {result_df.shape[0]} 个 ORF 的表达量")
        print(f"  耗时 {elapsed:.1f}s ({elapsed / max(n_orfs, 1):.3f}s per ORF)")
        print(f"  输出: {result_df.shape[0]} rows × {result_df.shape[1]} cols")

        return result_df


# ── CLI ───────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="从 RiboseQC P-site bedgraph 提取 per-ORF 表达量"
    )
    parser.add_argument(
        "--orf-meta", required=True, help="Unified ORF metadata TSV"
    )
    parser.add_argument(
        "--psites-dir", required=True, help="RiboseQC P-site bedgraph 目录"
    )
    parser.add_argument(
        "--sample-pattern",
        default="*_P_sites_plus.bedgraph",
        help="样本匹配 glob 模式 (默认: *_P_sites_plus.bedgraph)",
    )
    parser.add_argument(
        "--orf-confidence", default=None, help="ORF confidence TSV (可选)"
    )
    parser.add_argument(
        "--min-ocs",
        type=float,
        default=0.0,
        help="最低 OCS 阈值 (默认 0.0)",
    )
    parser.add_argument(
        "--max-orfs",
        type=int,
        default=0,
        help="最多处理 ORF 数 (0=全部, 按 OCS 降序取)",
    )
    parser.add_argument(
        "--workers", type=int, default=4, help="并行 worker 数 (默认 4)"
    )
    parser.add_argument(
        "--chrom-prefix",
        default="",
        help="染色体名前缀 (如 'chr', 空字符串表示不添加)",
    )
    parser.add_argument(
        "--output", required=True, help="输出文件路径"
    )
    args = parser.parse_args()

    quantifier = ORFExpressionQuantifier(
        orf_meta_path=args.orf_meta,
        psites_dir=args.psites_dir,
        sample_pattern=args.sample_pattern,
        orf_confidence_path=args.orf_confidence,
        min_ocs=args.min_ocs,
        max_orfs=args.max_orfs,
        workers=args.workers,
        chrom_prefix=args.chrom_prefix,
    )

    result_df = quantifier.run()

    if result_df is not None and len(result_df) > 0:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_path, sep="\t", index=False)
        print(f"\n保存结果至: {output_path}")

        # 打印统计
        print("\n=== 表达量摘要 ===")
        print(f"ORF 总数: {len(result_df)}")
        for s in quantifier.samples[:5]:
            col = f"{s}_reads"
            if col in result_df.columns:
                nonzero = (result_df[col] > 0).sum()
                total = result_df[col].sum()
                print(f"  {s}: {nonzero} non-zero, total={total:,} reads")

    else:
        print("WARNING: 没有输出结果")
        sys.exit(1)


if __name__ == "__main__":
    main()
