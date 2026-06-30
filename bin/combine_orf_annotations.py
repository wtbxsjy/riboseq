#!/usr/bin/env python3
"""
合并 ORF 注释: OCS + tier + 分类 + 表达量 + TE counts → 完整 ORF 注释表。

用法:
  combine_orf_annotations.py \
      --orf-meta unified_orfs.metadata.tsv \
      --orf-confidence orf_confidence.tsv \
      --expression expression_summary.tsv \
      --rpkm-tpm expression_rpkm_tpm.tsv \
      --te-counts merged_counts.tsv \
      --classification-dir orf_classification/ \
      --output orf_expression_combined.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_optional_tsv(path, label):
    """加载可选的 TSV 文件。"""
    if path is None:
        print(f"  {label}: 未提供, 跳过")
        return None
    p = Path(path)
    if not p.exists():
        print(f"  {label}: 文件不存在 ({p}), 跳过")
        return None
    df = pd.read_csv(p, sep="\t")
    print(f"  {label}: {df.shape[0]} rows × {df.shape[1]} cols")
    return df


def load_classification_files(classification_dir):
    """加载分类结果文件。"""
    d = Path(classification_dir)
    results = {}

    # GENCODE classification
    gencode_file = d / "gencode" / "gencode_results.orfs.out"
    if gencode_file.exists():
        try:
            df = pd.read_csv(gencode_file, sep="\t")
            # 提取 orf_id + orf_biotype
            if "orf_id" in df.columns and "orf_biotype" in df.columns:
                results["gencode"] = df[["orf_id", "orf_biotype"]]
                print(f"  GENCODE: {len(results['gencode'])} rows")
            else:
                print(f"  GENCODE: 缺少 orf_id/orf_biotype 列")
        except Exception as e:
            print(f"  GENCODE: 读取失败 ({e})")

    # ORFquant classification
    orfquant_file = d / "orfquant_classification.tsv"
    if orfquant_file.exists():
        try:
            df = pd.read_csv(orfquant_file, sep="\t")
            key_cols = ["orf_id"]
            cat_cols = [
                c for c in df.columns if "category" in c.lower()
            ]
            keep = key_cols + cat_cols
            if len(keep) > 1:
                results["orfquant"] = df[keep]
                print(f"  ORFquant: {len(results['orfquant'])} rows")
        except Exception as e:
            print(f"  ORFquant: 读取失败 ({e})")

    # ORF-type classification
    orftype_file = d / "orf_type" / "orftype_classification.tsv"
    if orftype_file.exists():
        try:
            df = pd.read_csv(orftype_file, sep="\t")
            if "orf_id" in df.columns:
                # 取 orf_id + classification 列
                cls_col = next(
                    (c for c in df.columns if "class" in c.lower()), None
                )
                if cls_col:
                    results["orftype"] = df[["orf_id", cls_col]]
                else:
                    results["orftype"] = df
                print(f"  ORF-type: {len(results['orftype'])} rows")
        except Exception as e:
            print(f"  ORF-type: 读取失败 ({e})")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="合并 ORF 注释为完整表"
    )
    parser.add_argument("--orf-meta", help="Unified ORF metadata")
    parser.add_argument("--orf-confidence", default=None, help="ORF confidence TSV")
    parser.add_argument("--expression", default=None, help="Expression summary TSV")
    parser.add_argument("--rpkm-tpm", default=None, help="RPKM/TPM TSV")
    parser.add_argument("--te-counts", default=None, help="TE merged counts TSV")
    parser.add_argument(
        "--classification-dir", default=None, help="orf_classification/ 目录"
    )
    parser.add_argument("--output", required=True, help="输出合并 TSV")
    args = parser.parse_args()

    print("=== 加载输入文件 ===")

    # 决定基础表: 优先用表达量 (如果有), 否则用 metadata
    expression_df = load_optional_tsv(args.expression, "Expression")
    rpkm_df = load_optional_tsv(args.rpkm_tpm, "RPKM/TPM")
    used_as_base = set()  # 记录哪些文件已作为基础表, 避免重复 merge

    if expression_df is not None:
        merged = expression_df.copy()
        used_as_base.add("expression")
        meta = load_optional_tsv(args.orf_meta, "ORF metadata")
        if meta is not None:
            extra_cols = ["orf_id"] + [c for c in meta.columns if c not in merged.columns and c != "orf_id"]
            merged = merged.merge(meta[extra_cols], on="orf_id", how="left")
        print(f"  基础表 (expression + meta): {merged.shape[0]} rows × {merged.shape[1]} cols")
    elif rpkm_df is not None:
        merged = rpkm_df.copy()
        used_as_base.add("rpkm_tpm")
        meta = load_optional_tsv(args.orf_meta, "ORF metadata")
        if meta is not None:
            extra_cols = ["orf_id"] + [c for c in meta.columns if c not in merged.columns and c != "orf_id"]
            merged = merged.merge(meta[extra_cols], on="orf_id", how="left")
        print(f"  基础表 (rpkm + meta): {merged.shape[0]} rows × {merged.shape[1]} cols")
    else:
        merged = load_optional_tsv(args.orf_meta, "ORF metadata")
        if merged is None:
            print("ERROR: 至少需要 --orf-meta, --expression, 或 --rpkm-tpm")
            sys.exit(1)

    # 合并 OCS (只保留需要的列, 避免重复)
    conf = load_optional_tsv(args.orf_confidence, "ORF confidence")
    if conf is not None:
        ocs_cols = ["orf_id"]
        for c in [
            "ocs", "tier", "s_translation", "s_agreement", "s_coverage",
            "s_periodicity", "s_readlevel", "detecting_tools", "n_detecting",
        ]:
            if c in conf.columns:
                ocs_cols.append(c)
        # 避免重复列
        existing = set(merged.columns)
        ocs_cols = [c for c in ocs_cols if c not in existing or c == "orf_id"]
        merged = merged.merge(conf[ocs_cols], on="orf_id", how="left")
        print(f"  合并后: {merged.shape[0]} rows × {merged.shape[1]} cols")

    # 合并表达量 (仅在未作为基础表时)
    if "expression" not in used_as_base and expression_df is not None:
        keep = ["orf_id"] + [
            c for c in expression_df.columns if c.endswith("_reads") or c.endswith("_pN")
        ]
        keep = [c for c in keep if c not in merged.columns or c == "orf_id"]
        merged = merged.merge(expression_df[keep], on="orf_id", how="left")

    # 合并 RPKM/TPM (仅在未作为基础表时)
    if "rpkm_tpm" not in used_as_base and rpkm_df is not None:
        keep = ["orf_id"] + [
            c for c in rpkm_df.columns
            if c.endswith("_rpkm") or c.endswith("_tpm") or c.endswith("_coverage_rpm")
        ]
        keep = [c for c in keep if c not in merged.columns or c == "orf_id"]
        merged = merged.merge(rpkm_df[keep], on="orf_id", how="left")

    # 合并 TE counts
    te = load_optional_tsv(args.te_counts, "TE counts")
    if te is not None:
        # TE 的 gene_id 对应 ORF ID; 加 te_ 前缀避免与表达量列名冲突
        te_renamed = te.rename(columns={"gene_id": "orf_id"})
        # 只保留 orf_id + 数值列 (移除 TE 文件可能重复的坐标列)
        te_cols = ["orf_id"] + [
            c for c in te_renamed.columns
            if c != "orf_id"
        ]
        te_renamed = te_renamed[te_cols].add_prefix("te_")
        te_renamed = te_renamed.rename(columns={"te_orf_id": "orf_id"})
        merged = merged.merge(te_renamed, on="orf_id", how="left")

    # 合并分类
    if args.classification_dir:
        print("\n=== 加载分类结果 ===")
        clf = load_classification_files(args.classification_dir)
        for source, df in clf.items():
            # 避免列名冲突
            rename_map = {}
            for c in df.columns:
                if c == "orf_id":
                    continue
                if c in merged.columns:
                    rename_map[c] = f"{source}_{c}"
            df = df.rename(columns=rename_map)
            keep = ["orf_id"] + [c for c in df.columns if c != "orf_id"]
            keep = [c for c in keep if c not in merged.columns or c == "orf_id"]
            merged = merged.merge(df[keep], on="orf_id", how="left")

    # 保存
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, sep="\t", index=False)
    print(f"\n=== 合并完成 ===")
    print(f"  输出: {output_path}")
    print(f"  {merged.shape[0]} rows × {merged.shape[1]} cols")
    print(f"  文件大小: {output_path.stat().st_size / 1024 / 1024:.1f} MB")

    # 快速统计
    if "ocs" in merged.columns:
        n_high = (merged["ocs"] >= 0.7).sum()
        n_medium = ((merged["ocs"] >= 0.4) & (merged["ocs"] < 0.7)).sum()
        print(f"  High (OCS≥0.7): {n_high}, Medium (0.4≤OCS<0.7): {n_medium}")

    read_cols = [c for c in merged.columns if c.endswith("_reads")]
    if read_cols:
        n_any = (merged[read_cols].sum(axis=1) > 0).sum()
        print(f"  ORFs with P-site reads: {n_any}/{len(merged)}")


if __name__ == "__main__":
    main()
