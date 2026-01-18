#!/bin/bash
set -euo pipefail

#############################################################################
# GENCODE M35 (Ensembl 112) 小鼠 Ensembl 目录准备脚本
#
# 用途：为 gencode-riboseqORFs 创建兼容的 Ensembl 注释目录
# 适用于：GENCODE M35 (GRCm39, Ensembl 112)
#
# 使用方法：
#   bash prepare_gencode_m35_ensembl.sh /path/to/output_dir
#############################################################################

# 颜色输出
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}GENCODE M35 Ensembl 目录准备脚本${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# 输出目录
OUTPUT_DIR=${1:-"./Ens_GENCODE_M35"}

if [[ -d "$OUTPUT_DIR" ]]; then
    echo -e "${YELLOW}[警告] 目录已存在: $OUTPUT_DIR${NC}"
    read -p "是否覆盖？(y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "取消操作"
        exit 0
    fi
    rm -rf "$OUTPUT_DIR"
fi

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo -e "${GREEN}[步骤 1/6] 下载 GENCODE M35 蛋白质序列...${NC}"
wget -q --show-progress \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.pc_translations.fa.gz \
    -O gencode.vM35.pc_translations.fa.gz

gunzip -f gencode.vM35.pc_translations.fa.gz
ln -sf gencode.vM35.pc_translations.fa PROTEOME_FASTA
echo -e "${GREEN}✓ 蛋白质序列: $(grep -c '^>' PROTEOME_FASTA) 条${NC}"

echo ""
echo -e "${GREEN}[步骤 2/6] 下载 GENCODE M35 转录本序列...${NC}"
wget -q --show-progress \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.transcripts.fa.gz \
    -O gencode.vM35.transcripts.fa.gz

gunzip -f gencode.vM35.transcripts.fa.gz
ln -sf gencode.vM35.transcripts.fa TRANSCRIPTOME_FASTA
echo -e "${GREEN}✓ 转录本序列: $(grep -c '^>' TRANSCRIPTOME_FASTA) 条${NC}"

echo ""
echo -e "${GREEN}[步骤 3/6] 下载 GENCODE M35 GTF 注释...${NC}"
wget -q --show-progress \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz \
    -O gencode.vM35.annotation.gtf.gz

gunzip -f gencode.vM35.annotation.gtf.gz

echo -e "${YELLOW}[步骤 3/6] 排序 GTF 文件...${NC}"
grep -v "^#" gencode.vM35.annotation.gtf | \
    sort -k1,1 -k4,4n > SORTED_TRANSCRIPTOME_GTF
echo -e "${GREEN}✓ 排序 GTF: $(wc -l < SORTED_TRANSCRIPTOME_GTF) 行${NC}"

echo ""
echo -e "${GREEN}[步骤 4/6] 提取转录本支持信息...${NC}"
cat > extract_transcript_support.py << 'PYTHON_SCRIPT'
#!/usr/bin/env python3
"""
从 GENCODE GTF 提取转录本支持信息
输出格式：Transcript stable ID | TSL | APPRIS
"""
import re
import sys

# 输出 header
print("Transcript stable ID\tTranscript Support Level (TSL)\tAPPRIS")

transcripts_seen = set()

with open("gencode.vM35.annotation.gtf") as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        # 只处理 transcript 行
        if fields[2] != "transcript":
            continue

        attrs = fields[8]

        # 提取 transcript_id（移除版本号）
        tid_match = re.search(r'transcript_id "([^"]+)"', attrs)
        if not tid_match:
            continue
        tid_full = tid_match.group(1)
        tid = tid_full.split('.')[0]  # 移除版本号

        # 避免重复
        if tid in transcripts_seen:
            continue
        transcripts_seen.add(tid)

        # 提取 TSL
        tsl_match = re.search(r'transcript_support_level "([^"]+)"', attrs)
        if tsl_match:
            tsl = tsl_match.group(1)
            # TSL 格式转换：NA, 1, 2, 3, 4, 5
            if tsl == "NA":
                tsl = "tslNA"
            else:
                tsl = f"tsl{tsl}"
        else:
            tsl = "tslNA"

        # GENCODE 小鼠注释通常没有 APPRIS，使用 NA
        # 但可以尝试从 tag 中提取
        if 'tag "APPRIS_principal"' in attrs:
            appris = "appris_principal"
        elif 'tag "APPRIS_candidate"' in attrs:
            appris = "appris_candidate"
        else:
            appris = "apprisNA"

        print(f"{tid}\t{tsl}\t{appris}")

print(f"# Total transcripts: {len(transcripts_seen)}", file=sys.stderr)
PYTHON_SCRIPT

python3 extract_transcript_support.py > TRANSCRIPT_SUPPORT 2> support.log
TRANSCRIPT_COUNT=$(tail -1 support.log | grep -o '[0-9]*')
echo -e "${GREEN}✓ 转录本支持: ${TRANSCRIPT_COUNT} 条${NC}"
rm support.log

echo ""
echo -e "${GREEN}[步骤 5/6] 生成 P-site BED 文件...${NC}"
cat > generate_psites.py << 'PYTHON_SCRIPT'
#!/usr/bin/env python3
"""
从 GTF 生成 P-site BED 文件
基于 start_codon 位置生成 P-site
"""
import re
import sys

psites = {}  # {transcript_id: (chrom, psite, strand)}

with open("SORTED_TRANSCRIPTOME_GTF") as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        feature = fields[2]

        # 只处理 start_codon 特征
        if feature == "start_codon":
            chrom = fields[0]
            start = int(fields[3])  # 1-based
            end = int(fields[4])
            strand = fields[6]

            # 提取 transcript_id
            tid_match = re.search(r'transcript_id "([^"]+)"', fields[8])
            if not tid_match:
                continue
            tid = tid_match.group(1)

            # P-site 位置：起始密码子的第一个核苷酸
            # BED 格式需要 0-based
            if strand == '+':
                psite_0based = start - 1
                psite_1based = start
            else:  # strand == '-'
                psite_0based = end - 1
                psite_1based = end

            # 保存（避免重复）
            if tid not in psites:
                psites[tid] = (chrom, psite_0based, psite_1based, strand)

# 输出 BED 格式（0-based，半开区间）
for tid, (chrom, psite_0, psite_1, strand) in sorted(psites.items()):
    print(f"{chrom}\t{psite_0}\t{psite_1}\t{tid}\t.\t{strand}")

print(f"# Generated {len(psites)} P-sites", file=sys.stderr)
PYTHON_SCRIPT

python3 generate_psites.py > PSITES_BED 2> psites.log
PSITE_COUNT=$(tail -1 psites.log | grep -o '[0-9]*')
echo -e "${GREEN}✓ P-site 位置: ${PSITE_COUNT} 个${NC}"
rm psites.log

echo ""
echo -e "${GREEN}[步骤 6/6] 验证文件完整性...${NC}"

# 验证所有必需文件
REQUIRED_FILES=(
    "PROTEOME_FASTA"
    "TRANSCRIPTOME_FASTA"
    "SORTED_TRANSCRIPTOME_GTF"
    "TRANSCRIPT_SUPPORT"
    "PSITES_BED"
)

ALL_OK=true
for file in "${REQUIRED_FILES[@]}"; do
    if [[ -f "$file" ]]; then
        SIZE=$(du -h "$file" | cut -f1)
        LINES=$(wc -l < "$file" 2>/dev/null || echo "N/A")
        echo -e "${GREEN}  ✓ $file${NC} ($SIZE, $LINES lines)"
    else
        echo -e "${RED}  ✗ $file (缺失)${NC}"
        ALL_OK=false
    fi
done

echo ""
echo -e "${GREEN}========================================${NC}"
if [[ "$ALL_OK" == "true" ]]; then
    echo -e "${GREEN}✅ 成功！Ensembl 目录已准备完成${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo "目录位置: $(pwd)"
    echo ""
    echo "文件清单："
    ls -lh PROTEOME_FASTA TRANSCRIPTOME_FASTA SORTED_TRANSCRIPTOME_GTF TRANSCRIPT_SUPPORT PSITES_BED
    echo ""
    echo "现在可以运行 GENCODE ORF mapper："
    echo -e "${YELLOW}  ./16_gencode_orf_mapper.sh --ensembl-dir $(pwd) ...${NC}"
else
    echo -e "${RED}❌ 出现错误，请检查上述缺失文件${NC}"
    exit 1
fi
