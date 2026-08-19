# Ensembl/GENCODE 分类目录手工构建（gffread 路线）

当物种不在 prepare_workflow.py 的 ENSEMBL_SPECIES_MAP 中、无网络、或已有本地
filtered GTF 时，用 gffread 从本地 GTF + 基因组 FASTA 构建 GENCODE mapper 所需的
分类目录。真实案例：maize（`run/maize/scripts/prepare_maize_ensembl_dir.sh`）。

## 目录要求

`--orf_classify_ensembl_dir` 指向的目录必须含这 5 个名字（mapper 硬编码）：

| 名称 | 内容 |
|---|---|
| `TRANSCRIPTOME_FASTA` | 转录组 FASTA（gffread -w 输出） |
| `SORTED_TRANSCRIPTOME_GTF` | 按 chrom + start 排序的 GTF |
| `PROTEOME_FASTA` | 蛋白组 FASTA（gffread -y 输出） |
| `TRANSCRIPT_SUPPORT` | transcript_id + TSL + APPRIS |
| `PSITES_BED` | 每条转录本 1 个 P-site 位置 |

（标准 Ensembl 下载目录里这 5 个是 symlink；手工构建时直接产出同名文件即可。）

## 五步构建

```bash
GTF=/path/to/species.genome.filtered.gtf     # 建议用 pipeline 产出的 filtered GTF
GENOME=/path/to/species.genome.fa
OUTDIR=/path/to/Ens58_<species>              # 命名跟随 Ensembl 风格
mkdir -p "$OUTDIR" && cd "$OUTDIR"

# 1. 转录组 FASTA（maize 实测 77,341 条）
gffread -w TRANSCRIPTOME_FASTA -g "$GENOME" "$GTF"

# 2. 蛋白组 FASTA（maize 实测 72,539 条）
gffread -y PROTEOME_FASTA -g "$GENOME" "$GTF"

# 3. 排序 GTF（自然顺序用 -k4,4n；maize 实测 130 万行）
grep -v "^#" "$GTF" | sort -k1,1 -k4,4n > SORTED_TRANSCRIPTOME_GTF

# 4. PSITES_BED：从 start_codon 行提取，+链 start-1、-链 end-1（BED 0-based），
#    每条转录本取 1 个位置（maize 实测 71,695 条）
# 5. TRANSCRIPT_SUPPORT：从 transcript 行生成，TSL=NA（BioMart 不可用时的兜底）
```

第 4、5 步的 Python 实现（按 maize 案例口径）：

```python
# generate_psites_and_support.py
import sys

gtf, out_psites, out_support = sys.argv[1:4]
support_rows = {}   # transcript_id -> tsl/appris
psite_rows = []     # (chrom, pos0, pos1, tx_id, strand)
with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        chrom, src, feat, start, end, _, strand, _, attrs = f
        attrs = {k: v.strip('"') for k, v in
                 (a.split(" ", 1) for a in attrs.rstrip(";").split("; ")) if k}
        tx = attrs.get("transcript_id")
        if feat == "transcript" and tx:
            support_rows[tx] = (tx, "NA", "NA")   # TSL, APPRIS
        elif feat == "start_codon" and tx:
            s, e = int(start), int(end)
            if strand == "+":
                psite_rows.append((chrom, s - 1, s, tx, strand))
            else:
                psite_rows.append((chrom, e - 1, e, tx, strand))
with open(out_support, "w") as fh:
    for t in sorted(support_rows):
        fh.write("\t".join(t) + "\n")
with open(out_psites, "w") as fh:
    for chrom, s, e, tx, strand in sorted(psite_rows):
        fh.write(f"{chrom}\t{s}\t{e}\t{tx}\t0\t{strand}\n")
```

## 三个必踩的坑

### 坑 1：gffread 蛋白 FASTA header 与 GTF protein_id 不匹配（最严重）

gffread `-y` 输出的 header 是 **transcript_id**（`>Zm00001eb000010_T001`），
而 GTF 里 ORF 的标识是 **protein_id**（`Zm00001eb000010_P001`）。
GENCODE mapper 按 protein_id 建立 protein_seq_map，匹配不上 → 所有
protein_coding 转录本被跳过。maize 实测后果：660K ORF 只分类出 2190 个，
其中 2187 个是 lncRNA。

**修复**：把 PROTEOME_FASTA 的 header 重命名为 GTF 的 protein_id：

```python
# rename_protein_headers.py  <proteome.fa> <gtf> <out.fa>
# 从 GTF 建 transcript_id -> protein_id 映射，然后重写 FASTA header
```

重命名后务必验证映射覆盖率（理想 100%）。

### 坑 2：BED6 未排序

交给 mapper 的 BED（含 PSITES_BED 及后续 classify 的输入）必须排序
（`sort -k1,1 -k2,2n`），否则分类结果错误——maize 案例 sort 后重跑才正确。

### 坑 3：TRANSCRIPT_SUPPORT 的 TSL 列

TSL/APPRIS 从 BioMart 下载失败时统一填 `NA` 是 mapper 可接受的兜底
（prepare_workflow.py 自动路线在 BioMart 失败时也这样做，日志会有 WARN）。
若 mapper 对 TSL 有额外期望（如过滤），先确认其解析代码再决定是否填值。

## 验证清单

- 5 个文件都存在且非空；GTF 行数、转录本数、蛋白数数量级合理
- `grep -c ">" PROTEOME_FASTA` 与 GTF `protein_id` 唯一数一致
- PSITES_BED 行数 ≈ 转录本数（每条转录本 1 行）
- 分类跑完检查 biotype 分布（⚠️ orf_biotype 是第 10 列且文件有 header，
  `tail -n +2 | cut -f10`）——lncRNA 占比异常高（>90%）说明坑 1 复发了；只剩
  lncRNA/CDS 两类则是版本后缀不匹配的 2 类塌缩（2026-08-19 已修，见
  riboseq-orf-analysis classifiers.md）
