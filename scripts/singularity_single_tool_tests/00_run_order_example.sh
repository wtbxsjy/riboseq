#!/usr/bin/env bash
#
# 顺序示例（从 BAM 开始，单工具脚本串联）
#
# 目的：把 scripts/singularity_single_tool_tests/ 下的每个单工具脚本按“常见依赖关系”串起来，
#      作为查阅/拷贝用的参考模板。
#
# 重要说明：
# - 这是“示例顺序 + 占位路径”，默认不保证在你当前环境能直接跑通。
# - 你需要把下面的输入路径变量改成你自己的：BAM / FASTA / GTF / FAI / rRNA FASTA。
# - 各工具的真实参数/镜像版本以本目录内的 01..11 脚本为准（它们 mirror pipeline modules）。
#
# 建议用法：
# - 仅查阅：直接打开此文件。
# - 想真的跑：把 DRY_RUN=0，并填写输入路径。

set -euo pipefail

###############################################################################
# 0) 全局配置
###############################################################################

# DRY_RUN=1 时只打印命令，不执行。
DRY_RUN=${DRY_RUN:-1}

run() {
  echo "+ $*"
  if [[ "$DRY_RUN" == "0" ]]; then
    eval "$@"
  fi
}

# 如果你的数据在 /mnt/c/...（WSL 常见），通常需要 bind。
# 不影响纯查阅；真正跑时建议打开。
export BIND_EXTRA=${BIND_EXTRA:-"/mnt:/mnt"}

# 目录约定：输出都写到本目录下的 out_* 子目录。
BASE_OUT=${BASE_OUT:-"$(pwd)"}

###############################################################################
# 1) 输入路径（你需要修改这里）
###############################################################################

SAMPLE=${SAMPLE:-"S1"}

# 必需输入
BAM=${BAM:-"/path/to/S1.bam"}                 # 输入 BAM（建议已排序）
FASTA=${FASTA:-"/path/to/genome.fa"}          # genome fasta
GTF=${GTF:-"/path/to/annotation.gtf"}         # annotation gtf
FAI=${FAI:-"/path/to/genome.fa.fai"}          # genome fasta index（samtools faidx 产生）
RRNA=${RRNA:-"/path/to/rrna.fa"}              # rp-bp 需要的 rRNA fasta

# 可选：如果你想在示例里也跑 TI（translation initiation）相关 BAM，可自行加变量。

###############################################################################
# 2) Step A: sORF predictors 前的 BAM 过滤
###############################################################################

# 过滤规则与 pipeline 一致（见 01_sorf_bam_filter.sh 以及 pipeline 的 sorf_* 参数）：
# - 只保留 primary mapped reads，排除：unmapped(0x4) secondary(0x100) duplicate(0x400) supplementary(0x800)
# - 按 SEQ 长度过滤（默认 28-30）
# - unique：auto（优先 NH==1，否则 MAPQ>=阈值），MAPQ 默认 60
# - 排除指定 contig（regex）

OUT_01="$BASE_OUT/out_01_filter"

# 默认参数：与 nextflow.config 的 params.sorf_* 保持一致
SORF_FILTER=${SORF_FILTER:-1}
SORF_UNIQUE_MODE=${SORF_UNIQUE_MODE:-"auto"}
SORF_UNIQUE_MAPQ=${SORF_UNIQUE_MAPQ:-60}
SORF_READ_LEN_MIN=${SORF_READ_LEN_MIN:-28}
SORF_READ_LEN_MAX=${SORF_READ_LEN_MAX:-30}
EXCLUDE_REGEX=${EXCLUDE_REGEX:-'^(chr)?(M|MT|Mt|chrM|chrMT|chrMt|ChrM|ChrMT|ChrMt)$|^(chr)?(C|CP|Pt|chrC|chrCP|chrPt|ChrC|ChrCP|ChrPt)$|^chrUn_.*|.*_random$|.*_alt$|.*_fix$'}

if [[ "${SORF_FILTER}" == "1" || "${SORF_FILTER}" == "true" || "${SORF_FILTER}" == "TRUE" ]]; then
  run "./01_sorf_bam_filter.sh \
    --sample '${SAMPLE}' \
    --bam '${BAM}' \
    --fai '${FAI}' \
    --unique-mode '${SORF_UNIQUE_MODE}' \
    --mapq '${SORF_UNIQUE_MAPQ}' \
    --len-min '${SORF_READ_LEN_MIN}' \
    --len-max '${SORF_READ_LEN_MAX}' \
    --exclude-regex '${EXCLUDE_REGEX}' \
    --outdir '${OUT_01}'"
else
  echo "[INFO] SORF_FILTER is disabled; using unfiltered BAM for predictors (example only)."
  mkdir -p "${OUT_01}"
  BAM_FILT="${BAM}"
fi

if [[ -z "${BAM_FILT:-}" ]]; then
  BAM_FILT="$OUT_01/${SAMPLE}.sorf.filtered.bam"
fi

###############################################################################
# 3) Step B: RiboseQC（准备 annotation，一次即可）
###############################################################################

OUT_02="$BASE_OUT/out_02_riboseqc_annot"

run "./02_riboseqc_prepareannotation.sh \
  --gtf '${GTF}' \
  --fasta '${FASTA}' \
  --outdir '${OUT_02}'"

# 注意：RiboseQC 的 *_Rannot 文件名取决于你的 GTF basename。
# 真实跑时请在 OUT_02 目录里找到那个 *_Rannot，并把变量指向它。
RANNOT=${RANNOT:-"$OUT_02/REPLACE_ME_Rannot"}

###############################################################################
# 4) Step C: RiboseQC analysis（建议 pre/post 各跑一次，方便对照）
###############################################################################

OUT_03="$BASE_OUT/out_03_riboseqc_prefilter"
OUT_04="$BASE_OUT/out_04_riboseqc_postfilter"

run "./03_riboseqc_analysis.sh \
  --sample '${SAMPLE}_prefilter' \
  --bam '${BAM}' \
  --annotation '${RANNOT}' \
  --fasta '${FASTA}' \
  --outdir '${OUT_03}' \
  --fast-mode TRUE"

run "./03_riboseqc_analysis.sh \
  --sample '${SAMPLE}_postfilter' \
  --bam '${BAM_FILT}' \
  --annotation '${RANNOT}' \
  --fasta '${FASTA}' \
  --outdir '${OUT_04}' \
  --fast-mode TRUE"

# ORFquant 需要 RiboseQC 的 *_for_ORFquant 文件。
FOR_ORFQUANT=${FOR_ORFQUANT:-"$OUT_03/${SAMPLE}_prefilter_for_ORFquant"}

###############################################################################
# 5) Step D: ORFquant（读取 RiboseQC 的 *_for_ORFquant + *_Rannot）
###############################################################################

OUT_05="$BASE_OUT/out_05_orfquant"

run "./04_orfquant_run.sh \
  --sample '${SAMPLE}' \
  --for-orfquant '${FOR_ORFQUANT}' \
  --annotation '${RANNOT}' \
  --fasta '${FASTA}' \
  --cpus 4 \
  --outdir '${OUT_05}'"

# 如果你所在环境无法访问 GitHub 下载 ORFquant 源码包，可在真实跑时加：
#   --orfquant-pkg /path/to/ORFquant-1.02.tar.gz

###############################################################################
# 6) Step E: RiboTISH（quality -> predict；predict 使用过滤后的 BAM）
###############################################################################

OUT_06="$BASE_OUT/out_06_ribotish_quality"
OUT_07="$BASE_OUT/out_07_ribotish_predict"

run "./05_ribotish_quality.sh \
  --sample '${SAMPLE}' \
  --bam '${BAM_FILT}' \
  --gtf '${GTF}' \
  --cpus 4 \
  --outdir '${OUT_06}'"

RIBOPARA=${RIBOPARA:-"$OUT_06/${SAMPLE}.para.py"}

run "./06_ribotish_predict.sh \
  --sample '${SAMPLE}' \
  --bam '${BAM_FILT}' \
  --gtf '${GTF}' \
  --fasta '${FASTA}' \
  --ribopara '${RIBOPARA}' \
  --cpus 4 \
  --outdir '${OUT_07}'"

###############################################################################
# 7) Step F: Ribotricer（prepare-orfs -> detect-orfs；detect 使用过滤后的 BAM）
###############################################################################

OUT_08="$BASE_OUT/out_08_ribotricer_index"
OUT_09="$BASE_OUT/out_09_ribotricer_detect"

# prepare-orfs 的 prefix 是“索引前缀”，与 sample 无关。
# 这里用 genome 只是示例。
RIBOTRICER_PREFIX=${RIBOTRICER_PREFIX:-"genome"}

run "./07_ribotricer_prepareorfs.sh \
  --prefix '${RIBOTRICER_PREFIX}' \
  --gtf '${GTF}' \
  --fasta '${FASTA}' \
  --outdir '${OUT_08}'"

RIBOTRICER_INDEX=${RIBOTRICER_INDEX:-"$OUT_08/${RIBOTRICER_PREFIX}_candidate_orfs.tsv"}

# strandedness 示例：forward/reverse/unstranded（unstranded 在 ribotricer 里曾有已知问题，脚本默认不传）。
STRANDEDNESS=${STRANDEDNESS:-"forward"}

run "./08_ribotricer_detectorfs.sh \
  --sample '${SAMPLE}' \
  --bam '${BAM_FILT}' \
  --index '${RIBOTRICER_INDEX}' \
  --stranded '${STRANDEDNESS}' \
  --outdir '${OUT_09}'"

###############################################################################
# 8) Step G: rp-bp（prepare-genome 一次 -> predict 每个 sample；predict 使用过滤后的 BAM）
###############################################################################

OUT_10="$BASE_OUT/out_10_rpbp_genome"
OUT_11="$BASE_OUT/out_11_rpbp_predict"

run "./09_rpbp_prepare_genome.sh \
  --fasta '${FASTA}' \
  --gtf '${GTF}' \
  --rrna '${RRNA}' \
  --cpus 8 \
  --outdir '${OUT_10}'"

ORFS_GEN=${ORFS_GEN:-"$OUT_10/transcript-index/genome.orfs-genomic.bed.gz"}
ORFS_EX=${ORFS_EX:-"$OUT_10/transcript-index/genome.orfs-exons.bed.gz"}

run "./10_rpbp_predict.sh \
  --sample '${SAMPLE}' \
  --bam '${BAM_FILT}' \
  --orfs-genomic '${ORFS_GEN}' \
  --orfs-exons '${ORFS_EX}' \
  --cpus 8 \
  --outdir '${OUT_11}'"

###############################################################################
# 9) Step H: RiboCode（可选：低深度数据失败很常见）
###############################################################################

OUT_12="$BASE_OUT/out_12_ribocode"

run "./11_ribocode_detect.sh \
  --sample '${SAMPLE}' \
  --bam '${BAM}' \
  --gtf '${GTF}' \
  --fasta '${FASTA}' \
  --stranded '${STRANDEDNESS}' \
  --outdir '${OUT_12}'"

###############################################################################
# 10) 完成提示
###############################################################################

echo ""
echo "[DONE] 这是顺序示例脚本。"
echo "- 当前 DRY_RUN=${DRY_RUN}（1=只打印；0=实际执行）"
echo "- 真正跑时请务必设置：BAM/FASTA/GTF/FAI/RRNA/RANNOT/FOR_ORFQUANT 等变量为真实路径。"
