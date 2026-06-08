#!/usr/bin/env python3
"""
从 RiboseQC coverage bedgraph (已 RPM 归一化) 计算每个 ORF 的 RPKM 和 TPM。

coverage bedgraph 的 value 字段已经是 RPM (reads per million mapped reads),
因此不需要额外获取 library size。

用法:
  calc_orf_rpkm_tpm.py \
      --expression expression_summary.tsv \
      --coverage-dir riboseqc/ \
      --sample-pattern "*_coverage_plus.bedgraph" \
      --output expression_rpkm_tpm.tsv \
      --workers 4
"""

import argparse
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd


class ORFRPKMCalculator:
    """从 RiboseQC coverage bedgraph 计算 ORF 的 RPKM/TPM。"""

    def __init__(
        self,
        expression_path,
        coverage_dir,
        sample_pattern="*_coverage_plus.bedgraph",
        workers=4,
        chrom_prefix="",
    ):
        self.expression_path = Path(expression_path)
        self.coverage_dir = Path(coverage_dir)
        self.sample_pattern = sample_pattern
        self.workers = workers
        self.chrom_prefix = chrom_prefix

        self.expr_df = None
        self.samples = []

    def load_expression(self):
        """加载表达量文件 (含 ORF 坐标)。"""
        print(f"加载表达量: {self.expression_path}")
        self.expr_df = pd.read_csv(self.expression_path, sep="\t")
        print(f"  {self.expr_df.shape[0]} ORFs × {self.expr_df.shape[1]} cols")

        for col in ["orf_id", "chrom", "start", "end", "strand"]:
            if col not in self.expr_df.columns:
                raise ValueError(f"表达量文件缺少必需列: {col}")

    def discover_coverage_samples(self):
        """扫描 coverage bedgraph 目录发现样本。"""
        # 优先用 expression 文件中的 reads 列推断
        read_cols = [c for c in self.expr_df.columns if c.endswith("_reads")]
        if read_cols:
            self.samples = [c.replace("_reads", "") for c in read_cols]
            print(f"从表达量文件中发现 {len(self.samples)} 个样本")
            return

        # 回退: 扫描 coverage bedgraph 文件
        files = sorted(self.coverage_dir.glob(self.sample_pattern))
        if not files:
            files = sorted(self.coverage_dir.glob("*_coverage_plus.bedgraph"))

        for f in files:
            name = f.name
            for suffix in ["_coverage_plus.bedgraph", "_coverage_minus.bedgraph"]:
                if name.endswith(suffix):
                    name = name[: -len(suffix)]
                    break
            if name and name not in self.samples:
                self.samples.append(name)

        print(f"从 bedgraph 文件发现 {len(self.samples)} 个样本")

    @staticmethod
    def query_coverage_region(bedgraph_file, chrom, start, end):
        """用 awk 查询 coverage bedgraph 中指定区间的总 RPM。"""
        if not Path(bedgraph_file).exists():
            return 0.0

        try:
            cmd = [
                "awk",
                "-F\\t",
                (
                    f'$1=="{chrom}" && $2<{end} && $3>{start} '
                    f"{{sum += $4 * ($3-$2); count += $3-$2}} "
                    f'END {{printf "%.6f\\t%d", sum, count}}'
                ),
                str(bedgraph_file),
            ]
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split("\t")
                return float(parts[0])
            return 0.0
        except Exception:
            return 0.0

    def process_single_orf_coverage(self, orf_row):
        """查询单个 ORF 在所有样本中的 coverage RPM 总和。"""
        result = {"orf_id": orf_row["orf_id"]}
        chrom = str(orf_row["chrom"])
        start = int(orf_row["start"])
        end = int(orf_row["end"])
        strand = str(orf_row["strand"])

        for sample in self.samples:
            if strand == "+":
                bg_file = (
                    self.coverage_dir / f"{sample}_coverage_plus.bedgraph"
                )
            else:
                bg_file = (
                    self.coverage_dir / f"{sample}_coverage_minus.bedgraph"
                )

            total_rpm = self.query_coverage_region(bg_file, chrom, start, end)
            result[f"{sample}_coverage_rpm"] = round(total_rpm, 4)

        return result

    def calculate_rpkm_tpm(self, coverage_df):
        """从 coverage RPM 值计算 RPKM 和 TPM。"""
        result_df = coverage_df.copy()

        # ORF 长度 (kb)
        # 需要从 expression 文件中获取坐标
        coord_lookup = {}
        for _, row in self.expr_df.iterrows():
            coord_lookup[row["orf_id"]] = (
                int(row["start"]),
                int(row["end"]),
            )

        orf_lengths = []
        for _, row in result_df.iterrows():
            s, e = coord_lookup.get(row["orf_id"], (0, 0))
            orf_lengths.append(abs(e - s))
        result_df["orf_length"] = orf_lengths
        result_df["orf_length_kb"] = [l / 1000.0 for l in orf_lengths]

        # RPKM: RPM_total / length_kb
        for sample in self.samples:
            rpm_col = f"{sample}_coverage_rpm"
            rpkm_col = f"{sample}_rpkm"
            if rpm_col in result_df.columns:
                result_df[rpkm_col] = (
                    result_df[rpm_col] / result_df["orf_length_kb"]
                )
                # handle zero-length edge case
                result_df[rpkm_col] = result_df[rpkm_col].fillna(0.0)
                result_df[rpkm_col] = result_df[rpkm_col].replace(
                    [np.inf, -np.inf], 0.0
                )

        # TPM: normalize so sum(TPM) = 1e6 per sample
        for sample in self.samples:
            rpkm_col = f"{sample}_rpkm"
            tpm_col = f"{sample}_tpm"
            if rpkm_col in result_df.columns:
                rpkm_sum = result_df[rpkm_col].sum()
                if rpkm_sum > 0:
                    result_df[tpm_col] = (
                        result_df[rpkm_col] / rpkm_sum * 1e6
                    )
                else:
                    result_df[tpm_col] = 0.0

        return result_df

    def run(self):
        """执行完整流程。"""
        t0 = time.time()

        self.load_expression()
        self.discover_coverage_samples()

        if not self.samples:
            print("ERROR: 未发现任何样本")
            return None

        # 并行提取 coverage RPM
        n_orfs = self.expr_df.shape[0]
        print(f"\n开始并行提取 {n_orfs} 个 ORF 的 coverage RPM...")

        results = []
        rows = [row for _, row in self.expr_df.iterrows()]

        with ProcessPoolExecutor(max_workers=self.workers) as exe:
            futures = {
                exe.submit(self.process_single_orf_coverage, row): i
                for i, row in enumerate(rows)
            }
            for fut in as_completed(futures):
                try:
                    res = fut.result()
                    results.append(res)
                except Exception as e:
                    idx = futures[fut]
                    print(f"  ERROR ORF idx={idx}: {e}")

        coverage_df = pd.DataFrame(results)

        # 计算 RPKM/TPM
        print("计算 RPKM 和 TPM...")
        result_df = self.calculate_rpkm_tpm(coverage_df)

        # 合并原始表达量数据中的坐标信息
        merge_cols = ["orf_id"]
        for c in ["chrom", "start", "end", "strand", "gene_id", "ocs", "tier"]:
            if c in self.expr_df.columns:
                merge_cols.append(c)
        result_df = self.expr_df[merge_cols].merge(
            result_df, on="orf_id", how="left"
        )

        elapsed = time.time() - t0
        print(f"\n完成! 耗时 {elapsed:.1f}s")
        print(f"  输出: {result_df.shape[0]} rows × {result_df.shape[1]} cols")

        return result_df


def main():
    parser = argparse.ArgumentParser(
        description="从 RiboseQC coverage bedgraph 计算 per-ORF RPKM/TPM"
    )
    parser.add_argument(
        "--expression", required=True, help="Phase 1 表达量 TSV"
    )
    parser.add_argument(
        "--coverage-dir", required=True, help="RiboseQC coverage bedgraph 目录"
    )
    parser.add_argument(
        "--sample-pattern",
        default="*_coverage_plus.bedgraph",
        help="样本匹配模式",
    )
    parser.add_argument(
        "--workers", type=int, default=4, help="并行 worker 数"
    )
    parser.add_argument(
        "--output", required=True, help="输出 RPKM/TPM TSV"
    )
    args = parser.parse_args()

    calc = ORFRPKMCalculator(
        expression_path=args.expression,
        coverage_dir=args.coverage_dir,
        sample_pattern=args.sample_pattern,
        workers=args.workers,
    )

    result_df = calc.run()

    if result_df is not None and len(result_df) > 0:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_path, sep="\t", index=False)
        print(f"\n保存结果至: {output_path}")

        # 统计
        rpkm_cols = [c for c in result_df.columns if c.endswith("_rpkm")]
        if rpkm_cols:
            print("\n=== RPKM 统计 ===")
            for col in rpkm_cols[:5]:
                vals = result_df[col]
                print(
                    f"  {col}: mean={vals.mean():.2f}, "
                    f"median={vals.median():.2f}, "
                    f"max={vals.max():.2f}, "
                    f"nonzero={(vals > 0).sum()}"
                )
    else:
        print("WARNING: 没有输出")
        sys.exit(1)


if __name__ == "__main__":
    main()
