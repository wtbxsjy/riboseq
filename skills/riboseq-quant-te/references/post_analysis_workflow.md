# post_analysis 两阶段过滤 + 批量 ggRibo（run_all.sh）

位置：`scripts/post_analysis/`（run_all.sh、config_template.yaml、README.md、
01_prelim_analysis.qmd、02_psite_filtering.qmd、03_generate_ggribo.qmd、
_functions.R、compute_psite_fast.py、compute_psite_purity.py）。
真实运行目录：`run/rice/post_analysis/`（project_config_with_orfquant.yaml +
output_with_orfquant_strandfix/ 是当前最新结果）。

## 运行

```bash
bash scripts/post_analysis/run_all.sh <project_config.yaml>
```

- 配置 YAML 作为唯一位置参数；`{result_dir}` 占位符自动展开
- 样本**不用显式指定**——compute_psite_fast.py 按 `*_P_sites_plus.bedgraph`
  文件名自动发现
- psite_purity.tsv 已存在则 Step2 跳过（"delete to recompute"）

## config_template.yaml 全键

```yaml
project:   {name, species, label}        # species 是 ggRibo 用学名
input:
  result_dir, metadata, confidence, expression, rpkm_tpm, unified_bed,
  gencode, riboseqc_dir, output_dir,
  bedgraph_plus_pattern: "*_P_sites_plus.bedgraph"
  bedgraph_minus_pattern: "*_P_sites_minus.bedgraph"
filtering:                               # 二阶段阈值（_functions.R 有同款兜底默认）
  prelim_reads_min: 9        # Stage1 每样本最小 reads
  prelim_pN_min: 0.5         # Stage1 每样本最小周期性
  prelim_cross_sample: 1     # Stage1 ≥N 样本同时满足
  p_site_gse_min: 9          # Stage2 每样本最小 P-site reads
  p_site_pct_min: 0.5        # Stage2 每样本最小 P-site 比例
  p_site_pos_min: 2          # Stage2 最小 P-site 位置偏移
  final_cross_sample: 2      # Stage2 ≥N 样本可重现
ggribo:
  max_per_biotype: 50
  top_n_samples: 10          # 每 ORF 展示 top N 样本
  sort_metric: "p_site_GSE"  # reads | p_site_GSE | p_site_pct | pN
  extend_bp: 200
  workers: 8
  biotypes: [uORF, dORF, doORF, intORF, uoORF, lncRNA, intergenic]  # null=全部非 CDS
```

## 四步流程

1. **Step1 初筛**：`quarto render 01_prelim_analysis.qmd -P config:$CONFIG`
   → `prelim_orfs_for_psite.bed`（BED12，col4=orf_id）+ `prelim_orfs.tsv`。
   逻辑（`apply_prelim_filters()`）：`orf_biotype != CDS` 且每样本
   `reads > 9 且 pN > 0.5`，≥1 样本满足。
2. **Step2 P-site 纯度**：
   ```bash
   python3 compute_psite_fast.py --bed prelim_orfs_for_psite.bed \
     --riboseqc-dir $RIBOSEQC_DIR --output psite_purity.tsv \
     --bedtools /usr/bin/bedtools --workers 8
   ```
   bedtools map 版：`sort -k1,1 -k2,2n bg | bedtools map -a BED -b - -c 4 -o sum -null 0`
   （strand-aware，plus bedgraph 只映射 plus 链 BED）。**28 秒完成 341K ORF × 23 样本**。
   输出列：`orf_id chrom start end strand` + 每样本
   `{s}_p_site_GSE/_reads_GSE/_p_site_pct/_p_site_pos/_not_p_site_GSE` + 全局
   `total_psites total_reads global_p_site_pct global_p_site_pos global_p_site_pos_sd n_samples_with_psites`
3. **Step3 P-site 过滤**：`02_psite_filtering.qmd` 合并 prelim_orfs.tsv +
   psite_purity.tsv → 每样本 `p_site_GSE>9 且 p_site_pct>0.5`，≥final_cross_sample
   样本 → 叠加 `length_aa ≤ 50` → `final_orfs_for_ggribo.tsv`（列：orf_id/chrom/start/
   end/strand/orf_biotype/gencode_biotype/total_reads_expr/max_pN/total_psites_real/
   mean_p_site_pct/tier/ocs/n_detecting/length_aa/is_classified）
4. **Step4 ggRibo 批量图**：`03_generate_ggribo.qmd` → `ggribo_plots/{biotype}/{orf_id}.png`
   + 临时 `tmp_gtf/`（可删）

## p_site_pct 口径（重要）

`02_psite_filtering.qmd` / `_functions.R::apply_psite_filters()` 里：
`pct = psite_purity 的 p_site_GSE / expression_summary 的 {s}_reads`（两整数同单位）。
**不用** compute_psite_*.py 自带 `reads_GSE`/`pct` 列（coverage bedgraph 求和 =
RPM 单位，做分母无意义——曾因此踩坑后弃用）。

## fast vs purity 两版 compute_psite

| | compute_psite_fast.py | compute_psite_purity.py |
|---|---|---|
| 依赖 | bedtools | 纯 Python（numpy+bisect+ProcessPoolExecutor） |
| 速度 | 28s / 341K×23 | 3.1h / 207K×23（16 workers） |
| p_site_pos | **恒 0.00**（不追踪位置） | 真实值（Σvalue×pos/Σvalue，负链从 ORF end 反向） |
| BED12 外显子块 | 整个区间 | 按 blocks 精确计算 |
| 适用 | run_all.sh 默认（Stage2 不用 pos 时可接受） | 需要 pos/位置加权时 |

## rice 实测数字（供预期管理）

- 输入 394,053 unified ORFs（含 210,761 intergenic）
- output_with_orfquant_strandfix（当前最新，final_cross_sample=1）：
  207,413 prelim → **10,255 final** → 10,255 张 ggRibo PNG
- 旧 README 数字（86,375 → 35,073 → 5,843 → 1,590）对应更早的 output_with_orfquant/
  运行，已过时；output/（旧单报告版）final = 6,188 行
- 文件体量：psite_purity.tsv 209MB（121 列）、prelim_orfs.tsv 18.4MB、
  step2 日志记录了 341,026 / 207,413 两次不同规模的运行

## 四级评分框架（与二阶段过滤互补，~/riboseq/post_analysis/ 旧体系）

`design_integrated_scoring.md` + `config_thresholds.yaml` + `scripts/integrate_orf_scores.R`：
综合分 = 0.25·S_purity + 0.30·S_expression + 0.25·S_confidence + 0.20·S_biology；
Tier：Gold≥0.70 / Silver≥0.50 / Bronze≥0.30 / <0.30=Weak。
S_purity 的 Yuanliang 规则：total_psites≥10、pct≥0.5、pos≥2.0、sd≤5.0。
legacy strict_sorf_pass 的 `max_pn≥5` 与二阶段 `prelim_pN_min:0.5` 是**同一列
不同量纲阈值**，勿混用。
