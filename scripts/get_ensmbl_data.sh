RELEASE=58
   ASSEMBLY=IRGSP-1.0
   SPECIES=oryza_sativa
   SPECIES_CAPS=Oryza_sativa
   OUTDIR=Ens${RELEASE}_rice
   BASE_URL=https://ftp.ensemblgenomes.org/pub/plants/release-${RELEASE}

   mkdir -p "$OUTDIR" && cd "$OUTDIR"

   GTF_GZ=${SPECIES_CAPS}.${ASSEMBLY}.${RELEASE}.gtf.gz
   CDNA_GZ=${SPECIES_CAPS}.${ASSEMBLY}.cdna.all.fa.gz
   NCRNA_GZ=${SPECIES_CAPS}.${ASSEMBLY}.ncrna.fa.gz
   PEP_GZ=${SPECIES_CAPS}.${ASSEMBLY}.pep.all.fa.gz

   wget -c ${BASE_URL}/gtf/${SPECIES}/${GTF_GZ}
   wget -c ${BASE_URL}/fasta/${SPECIES}/cdna/${CDNA_GZ}
   wget -c ${BASE_URL}/fasta/${SPECIES}/ncrna/${NCRNA_GZ} || true
   wget -c ${BASE_URL}/fasta/${SPECIES}/pep/${PEP_GZ}

   GTF=${GTF_GZ%.gz}
   CDNA=${CDNA_GZ%.gz}
   NCRNA=${NCRNA_GZ%.gz}
   PEP=${PEP_GZ%.gz}

   gunzip -c "$GTF_GZ" > "$GTF"
   gunzip -c "$CDNA_GZ" > "$CDNA"
   if [[ -f "$NCRNA_GZ" ]]; then gunzip -c "$NCRNA_GZ" > "$NCRNA"; fi
   gunzip -c "$PEP_GZ" > "$PEP"

   GTF_SORTED=${SPECIES_CAPS}.${ASSEMBLY}.${RELEASE}.sorted.gtf
   ( grep "^#" "$GTF" || true; grep -v "^#" "$GTF" | sort -k1,1 -k4,4n -k5,5n ) > "$GTF_SORTED"

   TRANSCRIPTOME=${SPECIES_CAPS}.${ASSEMBLY}.transcriptome.fa
   if [[ -f "$NCRNA" ]]; then
     awk 'BEGIN{FS="."} /^>/{print $1; next} {print}' "$CDNA" "$NCRNA" > "$TRANSCRIPTOME"
   else
     awk 'BEGIN{FS="."} /^>/{print $1; next} {print}' "$CDNA" > "$TRANSCRIPTOME"
   fi

   awk 'BEGIN{FS="."} /^>/{print $1; next} {print}' "$PEP" > "${PEP}.tmp" && mv "${PEP}.tmp" "$PEP"

   PSITES_BED=psites.bed
   python3 - <<'PY' "$GTF_SORTED" "$PSITES_BED"
   import sys
   gtf_file, output_bed = sys.argv[1], sys.argv[2]
   with open(gtf_file, 'r') as gtf, open(output_bed, 'w') as bed:
       for line in gtf:
           if line.startswith('#'):
               continue
           fields = line.rstrip().split('\t')
           if len(fields) < 9 or fields[2] != 'start_codon':
               continue
           chrom = fields[0]
           start = int(fields[3]) - 1
           end = int(fields[4])
           strand = fields[6]
           attrs = fields[8]
           transcript_id = ""
           for attr in attrs.split(';'):
               attr = attr.strip()
               if attr.startswith('transcript_id'):
                   transcript_id = attr.split('"')[1]
                   break
           if transcript_id:
               bed.write(f"{chrom}\t{start}\t{end}\t{transcript_id}\t.\t{strand}\n")
   PY

   TSL_FILE=transcript_support_level.txt
   python3 - <<'PY' "$GTF_SORTED" "$TSL_FILE"
   import re, sys
   gtf_file, output_file = sys.argv[1], sys.argv[2]
   transcripts = {}
   with open(gtf_file, 'r') as gtf:
       for line in gtf:
           if line.startswith('#'):
               continue
           fields = line.rstrip().split('\t')
           if len(fields) < 9 or fields[2] != 'transcript':
               continue
           attrs = fields[8]
           transcript_id = ''
           tsl = 'NA'
           appris = 'NA'
           for attr in attrs.split(';'):
               attr = attr.strip()
               if attr.startswith('transcript_id'):
                   transcript_id = attr.split('"')[1]
               elif 'transcript_support_level' in attr or 'tsl' in attr.lower():
                   match = re.search(r'(\d+)', attr)
                   if match:
                       tsl = match.group(1)
           if transcript_id:
               transcripts[transcript_id] = {'tsl': tsl, 'appris': appris}
   with open(output_file, 'w') as out:
       out.write("transcript_id\tTSL\tAPPRIS\n")
       for tid, data in transcripts.items():
           out.write(f"{tid}\t{data['tsl']}\t{data['appris']}\n")
   PY

   ln -sf "$GTF_SORTED" SORTED_TRANSCRIPTOME_GTF
   ln -sf "$TRANSCRIPTOME" TRANSCRIPTOME_FASTA
   ln -sf "$PEP" PROTEOME_FASTA
   ln -sf "$TSL_FILE" TRANSCRIPT_SUPPORT
   ln -sf "$PSITES_BED" PSITES_BED
