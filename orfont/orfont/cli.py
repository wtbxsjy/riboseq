#!/usr/bin/env python3
"""Unified CLI for ORF unification and classification.

Usage:
    orfont run --pipeline riboseq --ribotish a.txt --gtf ref.gtf --fasta ref.fa -o results
    orfont unify --ribotish a.txt --ribotricer b.tsv --gtf ref.gtf --fasta ref.fa
    orfont classify-gencode --bed unified.bed --metadata unified.metadata.tsv --ensembl-dir ./Ens110
    orfont classify-orfquant --gtf unified.gtf --metadata unified.metadata.tsv --annotation ref.gtf
    orfont classify-orftype --metadata unified.metadata.tsv --gtf ref.gtf
    orfont convert ribotish --predict a.txt --fasta ref.fa --study-id S1
    orfont convert ribotricer --tsv b.tsv --fasta ref.fa --study-id S1
"""

import sys
import argparse
import logging


def main():
    parser = argparse.ArgumentParser(
        description="ORF ontology: unification and classification for Ribo-seq ORFs")
    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # --- run (full pipeline) ---
    run_parser = subparsers.add_parser("run", help="Run complete pipeline (unify + classify)")
    run_parser.add_argument("--pipeline", "-p", default="riboseq",
                            choices=["riboseq"], help="Target pipeline")
    run_parser.add_argument("--ribotish", nargs="*", default=[],
                            help="Ribo-TISH predict files")
    run_parser.add_argument("--ribotricer", nargs="*", default=[],
                            help="Ribotricer TSV files")
    run_parser.add_argument("--ribocode", nargs="*", default=[],
                            help="RiboCode GTF files")
    run_parser.add_argument("--orfquant", nargs="*", default=[],
                            help="ORFquant GTF files")
    run_parser.add_argument("--gtf", required=True, help="Reference GTF")
    run_parser.add_argument("--fasta", required=True, help="Genome FASTA")
    run_parser.add_argument("--ensembl-dir", help="Ensembl annotation dir (for GENCODE)")
    run_parser.add_argument("--output", "-o", default="./orf_results")
    run_parser.add_argument("--prefix", default="unified_orfs")
    run_parser.add_argument("--no-gencode", action="store_true")
    run_parser.add_argument("--no-orfquant", action="store_true")
    run_parser.add_argument("--no-orftype", action="store_true")
    run_parser.add_argument("--no-frame-merge", action="store_true")
    run_parser.add_argument("--seq-cluster", action="store_true")
    run_parser.add_argument("--gencode-impl", default="original",
                            choices=["original", "fast", "indexed_fast"])
    run_parser.add_argument("--cpus", type=int, default=1)
    run_parser.add_argument("--verbose", "-v", action="store_true")

    # --- unify ---
    unify_parser = subparsers.add_parser("unify", help="Run ORF unification only")
    unify_parser.add_argument("--ribotish", nargs="*", default=[])
    unify_parser.add_argument("--ribotricer", nargs="*", default=[])
    unify_parser.add_argument("--ribocode", nargs="*", default=[])
    unify_parser.add_argument("--orfquant", nargs="*", default=[])
    unify_parser.add_argument("--gtf", required=True)
    unify_parser.add_argument("--fasta", required=True)
    unify_parser.add_argument("--output", "-o", default=".")
    unify_parser.add_argument("--prefix", default="unified_orfs")
    unify_parser.add_argument("--no-frame-merge", action="store_true")
    unify_parser.add_argument("--seq-cluster", action="store_true")
    unify_parser.add_argument("--verbose", "-v", action="store_true")

    # --- classify-gencode ---
    gencode_parser = subparsers.add_parser("classify-gencode", help="GENCODE classification")
    gencode_parser.add_argument("--bed", required=True, help="Unified BED12")
    gencode_parser.add_argument("--metadata", required=True, help="Unified metadata TSV")
    gencode_parser.add_argument("--ensembl-dir", required=True, help="Ensembl annotation dir")
    gencode_parser.add_argument("--output", "-o", default=".")
    gencode_parser.add_argument("--impl", default="original",
                                choices=["original", "fast", "indexed_fast"])
    gencode_parser.add_argument("--cpus", type=int, default=1)
    gencode_parser.add_argument("--verbose", "-v", action="store_true")

    # --- classify-orfquant ---
    orfquant_parser = subparsers.add_parser("classify-orfquant", help="ORFquant classification")
    orfquant_parser.add_argument("--gtf", required=True, help="Unified GTF")
    orfquant_parser.add_argument("--metadata", required=True, help="Unified metadata TSV")
    orfquant_parser.add_argument("--annotation", required=True, help="Reference GTF")
    orfquant_parser.add_argument("--output", "-o", default=".")
    orfquant_parser.add_argument("--verbose", "-v", action="store_true")

    # --- classify-orftype ---
    orftype_parser = subparsers.add_parser("classify-orftype", help="ORF-type classification")
    orftype_parser.add_argument("--metadata", required=True, help="Unified metadata TSV")
    orftype_parser.add_argument("--gtf", required=True, help="Reference GTF")
    orftype_parser.add_argument("--output", "-o", default=".")
    orftype_parser.add_argument("--cpus", type=int, default=1)
    orftype_parser.add_argument("--verbose", "-v", action="store_true")

    # --- convert ---
    convert_parser = subparsers.add_parser("convert", help="Convert tool outputs to GENCODE format")
    convert_sub = convert_parser.add_subparsers(dest="convert_tool")

    ct_ribotish = convert_sub.add_parser("ribotish", help="Convert Ribo-TISH")
    ct_ribotish.add_argument("--predict", required=True)
    ct_ribotish.add_argument("--fasta", required=True)
    ct_ribotish.add_argument("--study-id", required=True)
    ct_ribotish.add_argument("--output", "-o", default=".")
    ct_ribotish.add_argument("--min-length", type=int, default=16)
    ct_ribotish.add_argument("--verbose", "-v", action="store_true")

    ct_ribotricer = convert_sub.add_parser("ribotricer", help="Convert Ribotricer")
    ct_ribotricer.add_argument("--tsv", required=True)
    ct_ribotricer.add_argument("--fasta", required=True)
    ct_ribotricer.add_argument("--study-id", required=True)
    ct_ribotricer.add_argument("--output", "-o", default=".")
    ct_ribotricer.add_argument("--min-length", type=int, default=16)
    ct_ribotricer.add_argument("--min-phase-score", type=float, default=0.5)
    ct_ribotricer.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    level = logging.DEBUG if getattr(args, 'verbose', False) else logging.INFO
    logging.basicConfig(
        level=level,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')

    if args.command == "run":
        from orfont.pipelines.riboseq import run
        run(
            ribotish_files=args.ribotish or None,
            ribotricer_files=args.ribotricer or None,
            ribocode_files=args.ribocode or None,
            orfquant_files=args.orfquant or None,
            gtf_path=args.gtf,
            fasta_path=args.fasta,
            ensembl_dir=args.ensembl_dir,
            output_dir=args.output,
            unify_prefix=args.prefix,
            run_gencode=not args.no_gencode,
            run_orfquant=not args.no_orftype if hasattr(args, 'no_orftype') else True,
            run_orftype=not args.no_orftype,
            frame_merge=not args.no_frame_merge,
            seq_cluster=args.seq_cluster,
            gencode_impl=args.gencode_impl,
            cpus=args.cpus,
        )

    elif args.command == "unify":
        from orfont.unification.builder import unify
        unify(
            ribotish_files=args.ribotish or None,
            ribotricer_files=args.ribotricer or None,
            ribocode_files=args.ribocode or None,
            orfquant_files=args.orfquant or None,
            gtf_path=args.gtf,
            fasta_path=args.fasta,
            output_dir=args.output,
            prefix=args.prefix,
            frame_merge=not args.no_frame_merge,
            seq_cluster=args.seq_cluster,
        )

    elif args.command == "classify-gencode":
        from orfont.classification.wrapper import classify_gencode
        classify_gencode(
            bed_path=args.bed,
            metadata_path=args.metadata,
            ensembl_dir=args.ensembl_dir,
            output_dir=args.output,
            gencode_impl=args.impl,
            cpus=args.cpus,
        )

    elif args.command == "classify-orfquant":
        from orfont.classification.wrapper import classify_orfquant
        classify_orfquant(
            gtf_path=args.gtf,
            metadata_path=args.metadata,
            ref_gtf=args.annotation,
            output_dir=args.output,
        )

    elif args.command == "classify-orftype":
        from orfont.classification.wrapper import classify_orftype
        classify_orftype(
            metadata_path=args.metadata,
            ref_gtf=args.gtf,
            output_dir=args.output,
            cpus=args.cpus,
        )

    elif args.command == "convert":
        if args.convert_tool == "ribotish":
            from orfont.converters.ribotish import convert as conv_ribotish
            conv_ribotish(
                predict_file=args.predict,
                fasta_file=args.fasta,
                study_id=getattr(args, 'study_id'),
                output_dir=args.output,
                min_length=args.min_length,
            )
        elif args.convert_tool == "ribotricer":
            from orfont.converters.ribotricer import convert as conv_ribotricer
            conv_ribotricer(
                tsv_file=args.tsv,
                fasta_file=args.fasta,
                study_id=getattr(args, 'study_id'),
                output_dir=args.output,
                min_length=args.min_length,
                min_phase_score=args.min_phase_score,
            )
        else:
            print(f"Unknown convert tool: {args.convert_tool}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
