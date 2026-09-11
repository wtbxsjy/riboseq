#!/bin/bash
# ─── Per-read-length frame preference on an UNFILTERED BAM ──────────────────
# Usage: bash frame_preference_by_length.sh <bam> <gtf> <outdir> <sample> [subsample] [regex] [unique]
#
# Why: the pipeline's sORF filter (--sorf_read_len_min/max, default 28-30) runs
# BEFORE RiboseQC and riboWaltz, so both QC tools only ever see 28-30 nt reads
# and can never tell you whether the discarded lengths carried real footprint
# signal. This script measures, for every read length, how concentrated the
# read 5' ends are in one of the three CDS frames.
#
# Frame definition (important): the frame of a read is its position in the
# TRANSCRIPT's coding sequence, NOT its genomic offset from the CDS start —
# introns are not multiples of 3, so a genomic distance gives the wrong frame
# for every read past the first intron (measured: apparent periodicity 87%
# before the first intron vs 34% after it, purely an artefact). Each CDS exon
# therefore carries its cumulative coding offset in transcript order.
#
# Annotation scope: one transcript per gene (longest CDS). Otherwise a P-site
# covered by N transcripts of the same gene is counted N times.
#
# Output: per read length, the % of 5' ends in each frame + the max frame %.
# A constant per-length P-site offset shifts WHICH frame dominates but not the
# concentration, so max-frame % is a valid length-independent periodicity score.
set -euo pipefail

BAM="$1"; GTF="$2"; OUT="$3"; SAMPLE="$4"; FRAC="${5:-0.1}"
REGEX="${6:-}"; UNIQUE="${7:-no}"
mkdir -p "$OUT"
PREFIX="$OUT/$SAMPLE"

echo "[$SAMPLE] 1/4 CDS exons (with coding offsets) + origins from GTF ..." >&2
SEL="$PREFIX.canonical_tx.txt"
awk -F'\t' -v sel="$SEL" '
  $3 == "CDS" {
    if (!match($9, /transcript_id "[^"]+"/)) next
    tid = substr($9, RSTART+15, RLENGTH-16)
    gid = "."
    if (match($9, /gene_id "[^"]+"/)) gid = substr($9, RSTART+9, RLENGTH-10)
    len[tid] += $5 - $4 + 1
    gene[tid] = gid
  }
  END {
    for (t in len) {
      g = gene[t]
      if (!(g in best) || len[t] > len[best[g]]) best[g] = t
    }
    for (g in best) print best[g] > sel
  }
' "$GTF"
echo "  canonical transcripts (longest CDS per gene): $(wc -l < "$SEL")" >&2

awk -F'\t' -v OFS='\t' -v bed="$PREFIX.cds_exons.bed" -v org="$PREFIX.tx_origin.tsv" -v sel="$SEL" '
  BEGIN { while ((getline l < sel) > 0) keep[l] = 1; close(sel) }
  $3 == "CDS" {
    if (!match($9, /transcript_id "[^"]+"/)) next
    tid = substr($9, RSTART+15, RLENGTH-16)
    if (!(tid in keep)) next
    i = n[tid]++
    xs[tid, i] = $4; xe[tid, i] = $5      # 1-based inclusive GTF coords
    chrom[tid] = $1; strd[tid] = $7
  }
  END {
    for (t in n) {
      m = n[t]
      for (i = 0; i < m; i++) idx[i] = i
      for (i = 1; i < m; i++) {            # insertion sort by exon start
        v = idx[i]; j = i - 1
        while (j >= 0 && xs[t, idx[j]] > xs[t, v]) { idx[j+1] = idx[j]; j-- }
        idx[j+1] = v
      }
      cum = 0; origin = ""
      if (strd[t] == "+") {
        for (i = 0; i < m; i++) {
          e = idx[i]
          if (origin == "") origin = xs[t, e]
          print chrom[t], xs[t,e]-1, xe[t,e], t, "+", cum > bed
          cum += xe[t,e] - xs[t,e] + 1
        }
      } else {
        for (i = m - 1; i >= 0; i--) {
          e = idx[i]
          if (origin == "") origin = xe[t, e]
          print chrom[t], xs[t,e]-1, xe[t,e], t, "-", cum > bed
          cum += xe[t,e] - xs[t,e] + 1
        }
      }
      print t, origin, strd[t], cum > org
    }
  }
' "$GTF"
sort -k1,1 -k2,2n "$PREFIX.cds_exons.bed" -o "$PREFIX.cds_exons.bed"
echo "  CDS exons: $(wc -l < "$PREFIX.cds_exons.bed")" >&2

echo "[$SAMPLE] 2/4 read 5-prime ends (subsample $FRAC, regex=${REGEX:-none}, unique=$UNIQUE) ..." >&2
samtools view -F 0xD04 -s "$FRAC" "$BAM" | awk -F'\t' -v OFS='\t' \
    -v re="$REGEX" -v uniq="$UNIQUE" '
  {
    # same non-length criteria as SORF_BAM_FILTER (modules/local/sorf_bam_filter)
    if (re != "" && $3 ~ re) next
    if (uniq == "nh") {
      nh = ""
      for (i = 12; i <= NF; i++) if ($i ~ /^NH:i:/) { split($i, a, ":"); nh = a[3]; break }
      if (nh != 1) next
    }
    flag = $2; chrom = $3; pos = $4; c = $6
    ref = 0
    while (match(c, /[0-9]+[MDN=X]/)) {
      n = substr(c, RSTART, RLENGTH-1) + 0
      op = substr(c, RSTART+RLENGTH-1, 1)
      if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X") ref += n
      c = substr(c, RSTART+RLENGTH)
    }
    len = length($10)
    if (and(flag, 16)) { five = pos + ref - 1; strand = "-" }
    else               { five = pos;           strand = "+" }
    print chrom, five-1, five, len, strand
  }
' | sort -k1,1 -k2,2n > "$PREFIX.read5.bed"
echo "  reads: $(wc -l < "$PREFIX.read5.bed")" >&2

echo "[$SAMPLE] 3/4 intersect with CDS ..." >&2
bedtools intersect -a "$PREFIX.read5.bed" -b "$PREFIX.cds_exons.bed" -wa -wb \
  > "$PREFIX.read5_in_cds.tsv"

echo "[$SAMPLE] 4/4 frame tally (transcript coding offset) ..." >&2
awk -F'\t' -v OFS='\t' -v out="$PREFIX.frame_by_length.tsv" '
  {
    # A = read: 1 chrom 2 start0 3 five 4 len 5 strand
    # B = CDS:  6 chrom 7 start0 8 end 9 tid 10 strand 11 cumulative offset
    rlen = $4; rstrand = $5
    if (rstrand != $10) next              # same-strand CDS only
    if (rstrand == "+") within = $3 - ($7 + 1)
    else                within = $8 - $3
    if (within < 0) next
    f = ($11 + within) % 3
    cnt[rlen SUBSEP f]++; tot[rlen]++
  }
  END {
    print "read_length", "n_reads_in_CDS", "frame0_pct", "frame1_pct", "frame2_pct", "max_frame_pct" > out
    for (l in tot) {
      c0 = cnt[l SUBSEP 0] + 0; c1 = cnt[l SUBSEP 1] + 0; c2 = cnt[l SUBSEP 2] + 0
      p0 = 100*c0/tot[l]; p1 = 100*c1/tot[l]; p2 = 100*c2/tot[l]
      m = p0; if (p1 > m) m = p1; if (p2 > m) m = p2
      printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", l, tot[l], p0, p1, p2, m >> out
    }
  }
' "$PREFIX.read5_in_cds.tsv"
sort -k1,1n "$PREFIX.frame_by_length.tsv" -o "$PREFIX.frame_by_length.tsv"
echo "[$SAMPLE] DONE → $PREFIX.frame_by_length.tsv" >&2
