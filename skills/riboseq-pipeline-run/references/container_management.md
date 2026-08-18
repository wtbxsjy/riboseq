# 容器构建与管理

机器环境：singularity-ce 4.3.0 与 apptainer 1.5.2 并存；生产运行全部
`-profile singularity`；docker 仅测试/CI 用；另有 conda 双模式（见下）。

## 自定义容器清单（每项目一份，放 run/{project}/containers/）

| 参数 | 镜像 | 说明 |
|---|---|---|
| `--orfquant_container` | orfquant_patched.sif (~1.8GB) | 补丁版：NAMESPACE 冲突修复；另有 orfquant_mirai.sif（v1.3.2 mirai 并行） |
| `--unify_orf_container` | unify_orf.sif (~190MB) | ORF 统一（**不含 bedtools**） |
| `--gencode_orf_mapper_container` | gencode_orf_mapper.sif (~339MB) | 需 bedtools + BioPython |
| `--ribowaltz_container` | ribowaltz.sif (~430MB) | Bioc 3.20，缺 txdbmaker 时运行时安装 |
| `--price_container` | gedi.sif | PRICE/GEDI，从 Google Drive 挂载盘读 |
| `--rpbp_container` | rpbp.sif | 人类项目用 |

## 构建

```bash
# 标准构建（定义文件在 containers/Singularity.xxx.def）
apptainer build --fakeroot -F run/rice/containers/orfquant_patched.sif \
    containers/Singularity.orfquant.patched.def

# --fakeroot 卡住时：去掉重试；再不行就后台跑
nohup apptainer build -F out.sif containers/Singularity.xxx.def > /tmp/container_build.log 2>&1 &

# 构建结果不对时：旧源码被 apptainer 缓存
apptainer cache clean --force

# 直接从 Docker Hub pull
singularity pull --force docker://bioconductor/bioconductor_docker:RELEASE_3_20 /tmp/bioc_320.sif
singularity pull out.sif docker://bioconductor/bioconductor_docker:RELEASE_3_20
```

构建 Bioc 3.20 类镜像的已知失败史（不要重蹈）：Docker Hub 拉取超时（>2h）→
本地 bioc_3_20.sif 缺 RiboseQC 依赖（rmarkdown/DT/ggpubr/viridis）→ sed 补丁误删
函数头 → 缓存旧源码。最终方案：`Bootstrap: localimage From bioc_3_20.sif` 的
`Singularity.orfquant.bioc320.def`。

## 免重建技巧：运行时 bind 开发源码

```bash
export NXF_SINGULARITY_BINDPATH="/home/25119231r/riboseq/ORFquant:/opt/ORFquant"
```

把宿主机 ORFquant 开发源码 bind 进容器路径，容器内 `R CMD INSTALL -l /tmp/rlibs`
运行时安装即可生效——改 R 代码不用重建 sif。注意 `R_LIBS_USER` 要**追加**
（`${task.workDir}/Rlibs:${R_LIBS_USER}`）不能覆盖，否则容器预装包不可见。

## sandbox 与 OOM 限制（ORFquant 攻坚史）

- SIF 模式：Apptainer FUSE squashfuse 空闲超时杀容器（`SINGULARITY_FUSE_TIMEOUT=-1`
  无效）；exit 137
- sandbox：`apptainer build --sandbox out_sandbox in.sif` 解包绕开 FUSE，但
  **sandbox 内可见内存 ≈54GB（宿主 188GB）**——多进程方案要按 54GB 预算
- **最终可行方案**：sandbox + mclapply n_cores=1（每样本串行，无 fork/daemon）+ 
  MAX_PARALLEL=8 跨样本进程级并行（`run/maize/scripts/run_orfquant_parallel.sh`）
- 后续优化：mirai daemon 用 FaFile 磁盘句柄（每 daemon 2.5GB → 300MB）；
  主进程内存 100+GB → 6.3GB
- 判定 OOM：`dmesg | grep -i "oom\|killed"`、`cat /sys/fs/cgroup/memory.max`

## conda 双模式（绕开 Apptainer 全部问题）

动机："conda bypasses all Apptainer issues"。nextflow.config 已定义 `conda`
profile（`conda.enabled=true`，docker/singularity 关闭）：

```bash
nextflow run . -profile test,conda
```

ORFquant 相关参数：`orfquant_conda_env`、`orfquant_mirai_container`、
`orfquant_mirai_daemons`。SIF 模式与 conda 模式同套参数互为镜像。
