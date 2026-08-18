# ORFquant 非模式生物修复史（saga 摘要）

目的：rice/maize 等非模式生物没有预构建 BSgenome 包，ORFquant 报
`library(NULL)` 崩溃。以下按时间线摘要，细节只需在再次踩坑时回对话记录查
（ba515e20 会话 L463-L856、L1846-L2015、L4971-L7342；最终文档
`docs/MODEL_VS_NONMODEL_ORGANISMS.md`）。

## 第一战：load_annotation monkey-patch（5 轮迭代）

**症状**：`forge_BSgenome=FALSE` → `GTF_annotation$genome_package` 为 NULL →
`library(NULL)` 报 `'package' must be of length 1`。

**修复链（每轮踩一个新坑）**：

1. 第一版 patch 用 `load(envir=environment())` + `return(annotation)`——
   丢失原版 `<<-` 全局赋值 → `object 'GTF_annotation' not found`
2. 修正为复刻原逻辑：`GTF_annotation <- get(load(path))` + `GTF_annotation <<-`
   + `genome_seq <<-`（从 `singularity exec ... Rscript -e
   'print(ORFquant:::load_annotation)'` 读原码得出）
3. Groovy 里 R 的 `$` 被当插值 → 写成 `\$genome`；`cat()` 的 `\n` 也要 `\\n`
4. bash heredoc `<<RSCRIPTEOF` 未加引号 + `set -u` → `genome: unbound
   variable` → 改 `<<'RSCRIPTEOF'`（还误写过 `<'RSCRIPTEOF'`）
5. **两个入口都要 patch**：`PREPARE_FOR_ORFQUANT_CORRECTED`（prepare_for_ORFquant）
   和 `ORFQUANT_RUN`（run_ORFquant 内部也调 load_annotation）
6. `R_LIBS_USER` 要**追加**（`${task.workDir}/Rlibs:${R_LIBS_USER}`）不能覆盖，
   否则容器预装 ORFquant 不可见
7. 修复在远端同步时丢过两次（dev-rice_run 重建），需重新应用

## 第二战：BSgenome saga（~15 commits）

- monkey-patch 使 `genome_seq=NULL` 后，ORFquant 内部 `seqinfo(NULL)` →
  **99.94% 基因失败**（9579/9585）——不是 bug，是固有限制
- 尝试 `forge_BSgenome=TRUE`：BSgenomeForge 不支持自定义基因组 → revert
  （commit 1bafcff）；`.copySeqFile` 因 `inst/extdata` 不存在失败 → 补 patch
  自动建目录；2bit 构建用 `export(replaceAmbiguities(genome,new="N"), "genome.2bit")`
- `is(genome, "FaFile")` 对 `FaFile_Circ` 子类返回 FALSE（容器内 S4 注册不完整）
  → patch 改"先查 genome_package、再用 genome 字段兜底"
- Bioc 3.20 移除 `disjointExons`：ORFquant dev (v1.3.2) 还在调，master (v1.02.0)
  已改 `exonicParts(orfann, linked.to.single.gene.only=F)` → 补丁替换两处
- 容器构建失败史：R 3.6.2 装不了 dev（需 R≥4.0）→ Docker Hub Bioc 3.20 拉取
  超时 → rocker/r-base:4.4.2 失败 → 本地 bioc_3_20.sif 缺 RiboseQC 依赖 →
  sed 补丁误删函数头 → apptainer 缓存旧源码（`apptainer cache clean --force`）
  → 最终 `containers/Singularity.orfquant.bioc320.def`
  （`Bootstrap: localimage From bioc_3_20.sif`）

## 最终方案

- revert forge_BSgenome；FaFile monkey-patch
- 直接用 `~/riboseq/ORFquant` dev 源码（用户已 clone）
- 运行时从 bind-mount 的 `/opt/ORFquant` 安装（`R CMD INSTALL -l /tmp/rlibs`，
  配 `NXF_SINGULARITY_BINDPATH` 免重建容器，见 riboseq-pipeline-run skill）
- dev-rice_run → dev fast-forward 合并（7 文件 +384/−108）
- 文档：`docs/MODEL_VS_NONMODEL_ORGANISMS.md`

## 运行形态演化（OOM 攻坚）

| 方案 | 结果 |
|---|---|
| SIF + mirai 16 daemon 单样本 | exit 137（FUSE squashfuse 空闲超时） |
| sandbox + mirai | exit 137 OOM（sandbox 可见内存 ~54GB vs 宿主 188GB；16 daemon × 2.5GB） |
| sandbox + mclapply fork | R6 对象 fork 子进程析构崩溃 |
| **最终可行** | sandbox + mclapply n_cores=1 + MAX_PARALLEL=8 跨样本并行 |
| 后续优化 | mirai daemon 改 FaFile 磁盘句柄（2.5GB→300MB/daemon）；主进程 100+GB→6.3GB |

## 提示

- 现成补丁容器：orfquant_patched.sif（NAMESPACE 冲突修复）与
  orfquant_mirai.sif（v1.3.2 mirai 并行后端）
- 每样本正常耗时 5-11 分钟、峰值内存 17-18GB（请求 111GB）
