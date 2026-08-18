# ggRibo 绘图

包：hsinyenwu/ggRibo 的 clone（`~/riboseq/ggRibo`，分支 v2026.05.21 = 0.3.9.1；
系统已装 0.3.9；R ≥ 4.0 可用）。**没有 pipeline module，没有 Singularity 定义**
——评估结论是 ggRibo 保持 post-hoc 独立脚本（docs/integration_expression_ggribo_assessment.md）。

## 输入数据（二选一）

1. **`{sample}_ggribo.tsv`**（推荐）：RiboseQC 模块自动生成（optional channel），
   4 列无 header：`count \t chrom \t position \t strand`。position 是 **1-based**
   （bedgraph 0-based start → 取第 3 列）。
   生成逻辑（modules/local/riboseqc/analysis/main.nf L166-192）：
   ```bash
   awk -v OFS='\t' '{print $4, $1, $3, "+"}' "${prefix}_P_sites_plus.bedgraph"  > "${prefix}_ggribo.tsv"
   awk -v OFS='\t' '{print $4, $1, $3, "-"}' "${prefix}_P_sites_minus.bedgraph" >> "${prefix}_ggribo.tsv"
   ```
2. **bedgraph 直接给**：`list(plus="..._P_sites_plus.bedgraph", minus="..._P_sites_minus.bedgraph")`
   （create_seq_input 支持 bam/bigwig/bedgraph/tabular 四种扩展名）

RNA-seq 一律不画（所有用户脚本 `include_rna=FALSE`）。

## 两条调用路径

### A. 批量（post_analysis 标准流程）

- run_all.sh Step4 自动跑（见 post_analysis_workflow.md）；rice strandfix 产出
  10,255 张 PNG
- maize 脚本集（run/maize/scripts/）：`run_ggribo_batch.sh`（每批 50 ORF、top 10
  样本）、`run_full_ggribo.sh`（P-site 计算→两级过滤→绘图全流程）、
  `finish_psite_ggribo.sh`、`run_ggribo_parallel.sh`——maize 实测 **13,959 张 PNG**
  （~200KB/张）
- rice 并行版：`run/rice/scripts/run_ggribo_parallel.sh`（GNU parallel，12 workers，
  失败重试 1 次）；`run/rice/post_analysis/run_ggribo.R`（mclapply 16 workers，
  每 ORF 手写 4 行最小 GTF）

### B. 单基因 / 指定清单

`bin/plot_orf_ggribo.R`（optparse CLI，与 maize/rice 版同源）：

```bash
Rscript bin/plot_orf_ggribo.R \
  --orf-meta unified_orfs.metadata.tsv \
  --expression unified_orfs_expression_summary.tsv \
  --psites-dir result/riboseqc \
  --gtf unified_orfs.gtf \
  --orf-ids ORF_1_Os01g0100100,ORF_2_Os01g0100200 \   # 或 --orf-ids-file / --n-top-orfs
  --output-dir plots --extend 200 --n-samples-per-orf 3
```

demo 入口：`scripts/test_ggRibo.sh` / `test_ggRibo2.sh`（后者带
`--backend auto|ggribo|manual`，manual 是纯 ggplot2/patchwork 自绘无需 ggRibo）。

## 主模式实现要点（post_analysis/scripts/plot_orf_ggribo.R，借鉴 Yuanliang）

```r
# 1. 自定义 gtf_import（不调 ggRibo::gtf_import，避免 BiocGenerics 命名空间冲突）
gtf_import_custom(gtf_path, format="gtf", dataSource="AMP", organism=org_species)
#    → 写全局 Txome_Range（最小 Range_info R6 类：exonsByTx/txByGene/cdsByTx/
#      fiveUTR/threeUTR/tx_to_gene，用 txdbmaker::makeTxDbFromGFF 构建）
# 2. 信号输入 + 全局变量（ggRibo 函数隐式依赖）
seq_result <- create_seq_input(ribo_files=bg_file_list, sample_names=sample_name_list,
                               include_rna=FALSE)
assign("inputs_full", seq_result, envir=.GlobalEnv)
# 3. 画图
p <- ggRibo(gene_id=ids$gene_id, tx_id=ids$tx_id, Extend=200, NAME=title_str,
            Riboseq=seq_result$Riboseq, SampleNames=sample_name_list,
            GRangeInfo=Txome_Range, data_types=rep("Ribo-seq", n_samples),
            Y_scale="each", plot_genomic_direction=TRUE, show_seq=FALSE,
            ribo_linewidth=0.6)
ggsave(out_file, p, width=14, height=max(5, n_samples*2.5), dpi=150, limitsize=FALSE)
```

`ggRibo_tx` 是转录本（exon-spliced）坐标版；`ggRibo_decom` 是三框分解版
（额外参数 `plot_unassigned_reads, frame_logic, nth_sample`）。

## 必踩的坑

1. **隐式全局变量**：`ggRibo()` 默认实参引用 `inputs_full$RNAseq/Riboseq`、
   `Samples`、`Txome_Range`、`RNAseqBamPairorSingle` —— 必须先
   `assign("Txome_Range", ..., .GlobalEnv)` + `assign("inputs_full", ..., .GlobalEnv)`
   （plot_orf_ggribo.R 头部注释明确写了这条）
2. **GTF 缺行**：unified_orfs.gtf 没有 gene/transcript 行，makeTxDbFromGFF 解析失败
   → 每 ORF 生成最小 GTF：1 gene + 1 transcript + CDS/exon（多外显子按
   metadata `exon_blocks` 展开、含 CDS frame；`_functions.R::create_orf_gtf` /
   demo 的 `augment_gtf_for_txdb()` 都做了这件事）
3. **命名空间冲突**：ggRibo 与 ORFquant/RiboseQC 有同源 BiocGenerics 冲突
   （assessment 文档 247-248 行列为已知问题）→ 不要调用 ggRibo 自己的
   `gtf_import`；用自定义 gtf_import_custom，或 demo 的
   `get("Range_info", envir=asNamespace("ggRibo"))` 取类
4. **版本不一致**：clone 0.3.9.1 vs 已装 0.3.9，`library(ggRibo)` 加载的是已装版
   （跑前确认行为一致，必要时 `devtools::load_all("~/riboseq/ggRibo")`）
5. **bedgraph 坐标制**：ggRibo 用 1-based position（取 bedgraph 第 3 列 end）
6. **长文件名**：maize 脚本对超长 ORF ID 截断 80 字符 + digest 哈希
7. **ORF 未找到**：metadata/GTF 里查不到 orf_id 时报
   `"ORF %s not found in %s"`；样本 bedgraph 缺失时跳过或报错
8. **Y_scale/data_types 校验**：`Y_scale` 只能 "all"/"each"；`data_types` 长度必须
   等于样本数；`GRangeInfo` 必须提供

## 输出

- 批量：`ggribo_plots/{biotype}/{orf_id}.png`（width 12-14、height 随样本数
  2.5/样本、dpi 150、`limitsize=FALSE`）；`plot_unified_orfs_ggribo.R` 还支持合并
  `all_selected_orfs.pdf`
- 历史样例：~/riboseq/post_analysis/coverage_plots/（32 张 rice + 30 张 maize）、
  sorf_strict_ggribo/（1,113 张 maize strict sORF）、
  run/maize/result/ggribo_plots/（13,959 张）
