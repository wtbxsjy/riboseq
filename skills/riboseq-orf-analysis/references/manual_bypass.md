# 手动 bypass：跳过 pipeline 跑 unify / 分类

场景：pipeline 已跑完各工具 per-sample 预测（或部分），但 unify/分类失败、
想换参数重跑、或根本没跑 pipeline（human_GSE157490 62 样本的真实做法）。

## 入口：scripts/run_orf.py

自动检测 `orfont/`（项目根下优化包）：可用则走 DuckDB 优化路径，
否则退回原脚本子进程。子命令：`unify` / `classify-gencode` /
`classify-orftype` / `classify-orfquant`。

```bash
# unify（orfont 路径参数；原脚本路径同义）
python3 scripts/run_orf.py unify \
    --gtf human_SARS2.genome.filtered.gtf \
    --fasta human_SARS2.genome.fa \
    --output unified_orfs \
    --min_len 6 --threads 2 --frame-merge-min-overlap 0.9 \
    --ribotish <ribotish/postfilter>/*_pred.txt \
    --ribotricer <ribotricer/postfilter>/*_translating_ORFs.tsv \
    --ribocode <ribocode>/*_collapsed.gtf.gz \
    --orfquant <orfquant>/*_Detected_ORFs.gtf.gz \
    --price <price>/*.orfs.tsv

# 分类（注意：wrapper CLI 是 --input/--output_dir，不是 --bed/--metadata）
python3 scripts/run_orf.py classify-gencode \
    --input unified_orfs --output_dir out_gencode \
    --ensembl_dir Ens58_oryza_sativa --cpus 16
python3 scripts/run_orf.py classify-orftype \
    --input unified_orfs --output_dir out_orftype --gtf ref.gtf
```

**orfont unify 额外参数**：`--duckdb-db`（持久化 DB 文件）、
`--duckdb-memory-limit`（默认 32GB）、`--bedgraph-dir`/`--sample-list`
（stream 统计）、`--per-tool-output`。未移植的功能（frame-merge/seq-cluster/
bedgraph 统计/per-tool 输出）会自动 fallback 到原脚本。

## 容器内执行（unify_orf.sif 模式）

容器缺 duckdb（orfont 依赖），需先 `pip install --user duckdb`。标准姿势
（human_GSE157490 实际命令的简化）：

```bash
singularity exec --no-home --pid -B /home/25119231r/riboseq \
    run/{project}/containers/unify_orf.sif bash -c '
export HOME="$RUNTIME_DIR"; export PYTHONUSERBASE="$RUNTIME_DIR/.pylibs"
export PATH="$PYTHONUSERBASE/bin:$PATH"; export PIP_NO_CACHE_DIR=1
pip install --user --no-cache-dir duckdb
python3 -u /home/25119231r/riboseq/riboseq/scripts/run_orf.py unify ...
'
```

要点：`--no-home --pid` 隔离容器家目录；`-B` bind 仓库根；`PYTHONUSERBASE`
指到宿主机可写目录（容器 FS 只读，不能装到默认位置）。

**用自带包装脚本**（自动处理上述细节）：

```bash
bash skills/riboseq-orf-analysis/scripts/run_orf_in_container.sh \
  ~/riboseq/run/rice/containers/unify_orf.sif \
  -- unify --gtf ... --fasta ... --output unified_orfs --ribotish ...
# 包装脚本选项：--runtime-dir DIR（默认 <cwd>/orf_runtime）、--bind PATH、
# --no-duckdb、--dry-run、--log FILE
```

## orfont 实测与坑（human_GSE157490）

- 老版 DuckDB 全量内存 insert 在 4-6 个 Ribo-TISH 文件后必崩（任何配置）
  → 修复为流式分批 insert（500K/批）→ 8 个 RPF 20.5min、6GB RSS（15x 加速）
- RiboCode `.gtf.gz` 的 `endswith('.gtf')` 漏匹配 → 62 文件全部
  "missing required columns"（已修 `_open()` fallback）
- GENCODE 分类 982K ORFs 时 mapper 崩溃 → 分批 + `csv.field_size_limit(sys.maxsize)`
- 后台跑法：`nohup singularity exec ... > unify.log 2>&1 &`

## 各工具输入文件速查

| 工具 | 文件模式 | 位置 |
|---|---|---|
| Ribo-TISH | `{sample}_pred.txt` | result/orf_predictions/ribotish/postfilter/ |
| Ribotricer | `{sample}_translating_ORFs.tsv` | result/orf_predictions/ribotricer/postfilter/ |
| ORFquant | `{sample}_Detected_ORFs.gtf.gz` | result/orf_predictions/orfquant/ |
| PRICE | `{sample}.orfs.tsv` | result/orf_predictions/price/ |
| RiboCode | `{sample}_collapsed.gtf.gz`（优先）等 7 类 | result/orf_predictions/ribocode/ |
