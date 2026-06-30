#!/bin/bash
set -euo pipefail

###############################################################################
# GENCODE ORF Mapper 容器构建脚本
#
# 用途：构建 gencode-riboseqORFs 的 Singularity 容器镜像
# 输出：gencode-orf-mapper.sif (可分发到其他服务器)
#
# 使用方法：
#   bash build_container.sh [--output /path/to/output.sif] [--build-method METHOD]
#
# 参数：
#   --output PATH      输出镜像路径 [默认: ./gencode-orf-mapper.sif]
#   --build-method     构建方法: sudo/fakeroot/remote [默认: 自动检测]
#   --force           强制重新构建（覆盖已存在的镜像）
#   --test            构建后运行测试
#
###############################################################################

# 颜色输出
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 默认参数
OUTPUT_SIF="${OUTPUT_SIF:-gencode-orf-mapper.sif}"
BUILD_METHOD="${BUILD_METHOD:-auto}"
FORCE_BUILD=false
RUN_TEST=false
DEF_FILE="Singularity.def"

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --output)
            OUTPUT_SIF="$2"
            shift 2
            ;;
        --build-method)
            BUILD_METHOD="$2"
            shift 2
            ;;
        --force)
            FORCE_BUILD=true
            shift
            ;;
        --test)
            RUN_TEST=true
            shift
            ;;
        -h|--help)
            grep "^#" "$0" | grep -v "#!/bin/bash" | sed 's/^# //' | sed 's/^#//'
            exit 0
            ;;
        *)
            echo -e "${RED}错误: 未知参数 $1${NC}"
            echo "使用 --help 查看帮助"
            exit 1
            ;;
    esac
done

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}GENCODE ORF Mapper 容器构建${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# 检查定义文件
if [[ ! -f "$DEF_FILE" ]]; then
    echo -e "${RED}错误: 找不到 Singularity 定义文件: $DEF_FILE${NC}"
    exit 1
fi

echo -e "${BLUE}[配置]${NC}"
echo "  定义文件: $DEF_FILE"
echo "  输出镜像: $OUTPUT_SIF"
echo "  构建方法: $BUILD_METHOD"
echo ""

# 检查输出文件是否存在
if [[ -f "$OUTPUT_SIF" ]] && [[ "$FORCE_BUILD" == "false" ]]; then
    echo -e "${YELLOW}警告: 输出文件已存在: $OUTPUT_SIF${NC}"
    read -p "是否覆盖？(y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "取消构建"
        exit 0
    fi
    rm -f "$OUTPUT_SIF"
fi

# 检测可用的构建方法
detect_build_method() {
    if [[ "$BUILD_METHOD" != "auto" ]]; then
        echo "$BUILD_METHOD"
        return
    fi

    echo -e "${BLUE}[检测] 自动检测可用的构建方法...${NC}" >&2

    # 检查是否有 sudo 权限
    if sudo -n true 2>/dev/null; then
        echo -e "${GREEN}  ✓ 检测到 sudo 权限${NC}" >&2
        echo "sudo"
        return
    fi

    # 检查是否支持 fakeroot
    if singularity build --help 2>&1 | grep -q "fakeroot"; then
        echo -e "${GREEN}  ✓ 检测到 fakeroot 支持${NC}" >&2
        echo "fakeroot"
        return
    fi

    # 检查是否支持 remote build
    if singularity remote list 2>/dev/null | grep -q "SylabsCloud"; then
        echo -e "${YELLOW}  ! 使用远程构建（需要 Sylabs Cloud 账户）${NC}" >&2
        echo "remote"
        return
    fi

    echo -e "${RED}  ✗ 未检测到可用的构建方法${NC}" >&2
    echo "none"
}

# 执行构建
build_container() {
    local method=$1
    local output=$2
    local def=$3

    case "$method" in
        sudo)
            echo -e "${BLUE}[构建] 使用 sudo 方法构建...${NC}"
            sudo singularity build "$output" "$def"
            ;;
        fakeroot)
            echo -e "${BLUE}[构建] 使用 fakeroot 方法构建...${NC}"
            singularity build --fakeroot "$output" "$def"
            ;;
        remote)
            echo -e "${BLUE}[构建] 使用远程构建...${NC}"
            echo -e "${YELLOW}  注意: 需要登录 Sylabs Cloud${NC}"
            singularity build --remote "$output" "$def"
            ;;
        none)
            echo -e "${RED}错误: 没有可用的构建方法${NC}"
            echo ""
            echo "请尝试以下方法之一："
            echo "  1. 使用 sudo: sudo bash build_container.sh"
            echo "  2. 使用 fakeroot: bash build_container.sh --build-method fakeroot"
            echo "  3. 使用远程构建: bash build_container.sh --build-method remote"
            echo ""
            echo "或在有 sudo 权限的服务器上构建后传输镜像文件"
            exit 1
            ;;
        *)
            echo -e "${RED}错误: 未知的构建方法: $method${NC}"
            exit 1
            ;;
    esac
}

# 检测构建方法
DETECTED_METHOD=$(detect_build_method)

if [[ "$DETECTED_METHOD" == "none" ]]; then
    exit 1
fi

echo ""
echo -e "${BLUE}[开始构建]${NC}"
echo "  时间: $(date)"
echo "  方法: $DETECTED_METHOD"
echo ""

# 执行构建
START_TIME=$(date +%s)

if build_container "$DETECTED_METHOD" "$OUTPUT_SIF" "$DEF_FILE"; then
    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))

    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}✅ 构建成功！${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo "  镜像文件: $OUTPUT_SIF"
    echo "  文件大小: $(du -h "$OUTPUT_SIF" | cut -f1)"
    echo "  构建时间: ${DURATION} 秒"
    echo ""

    # 运行测试
    if [[ "$RUN_TEST" == "true" ]]; then
        echo -e "${BLUE}[测试] 验证容器...${NC}"

        echo "  测试 1: 检查 Python 版本"
        singularity exec "$OUTPUT_SIF" python3 --version

        echo "  测试 2: 检查 bedtools"
        singularity exec "$OUTPUT_SIF" bedtools --version

        echo "  测试 3: 检查 gffread"
        singularity exec "$OUTPUT_SIF" gffread --version

        echo "  测试 4: 检查 Biopython"
        singularity exec "$OUTPUT_SIF" python3 -c "import Bio; print('Biopython:', Bio.__version__)"

        echo "  测试 5: 检查主脚本"
        singularity exec "$OUTPUT_SIF" python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py --help | head -5

        echo ""
        echo -e "${GREEN}  ✓ 所有测试通过${NC}"
    fi

    echo -e "${YELLOW}下一步：${NC}"
    echo "  1. 测试镜像: bash build_container.sh --test"
    echo "  2. 分发镜像: scp $OUTPUT_SIF user@server:/path/to/containers/"
    echo "  3. 使用镜像: singularity exec $OUTPUT_SIF python3 /opt/gencode-riboseqORFs/ORF_mapper_to_GENCODE_v1.1.py --help"
    echo ""
else
    echo ""
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}❌ 构建失败${NC}"
    echo -e "${RED}========================================${NC}"
    echo ""
    echo "请检查上述错误信息"
    exit 1
fi
