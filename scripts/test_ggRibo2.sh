Rscript ./R/plot_unified_orfs_ggribo.R \
	--metadata /home/25119231r/MyDrive/sORF_Discovery_Project/results/mouse_Mucosal_Immunity/orf_unification/unified_orfs.metadata.tsv \
  --gtf /home/25119231r/MyDrive/sORF_Discovery_Project/results/mouse_Mucosal_Immunity/orf_unification/unified_orfs.gtf \
  --riboseqc-dir /home/25119231r/MyDrive/sORF_Discovery_Project/results/mouse_Mucosal_Immunity/riboseqc \
  --outdir ../test_results/ggribo_backend_optimized_demo \
  --orf-ids ORF_8_ENSMUSG00000033793.13,ORF_6_ENSMUSG00000025903.15 \
  --backend manual \
  --signal unique \
  --format png

