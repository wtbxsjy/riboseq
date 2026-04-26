#!/usr/bin/env python3

import importlib.util
import tempfile
import textwrap
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("unify_orf_predictions.py")
SPEC = importlib.util.spec_from_file_location("unify_orf_predictions", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class InferSampleIdFromPredictionPathTests(unittest.TestCase):
    def test_preserves_dots_in_ribotish_sample_names(self):
        sample_id = MODULE.infer_sample_id_from_prediction_path(
            "/tmp/E15.5_Brain_Ribo_pred.txt",
            "_pred.txt",
        )
        self.assertEqual(sample_id, "E15.5_Brain_Ribo")

    def test_preserves_dots_in_merged_orfquant_sample_names(self):
        sample_id = MODULE.infer_sample_id_from_prediction_path(
            "/tmp/E15.5_Brain_Ribo_merged_Detected_ORFs.gtf",
            "_Detected_ORFs.gtf",
        )
        self.assertEqual(sample_id, "E15.5_Brain_Ribo_merged")

    def test_distinguishes_samples_that_only_share_prefix_before_dot(self):
        brain = MODULE.infer_sample_id_from_prediction_path(
            "/tmp/E15.5_Brain_Ribo_translating_ORFs.tsv",
            "_translating_ORFs.tsv",
        )
        liver = MODULE.infer_sample_id_from_prediction_path(
            "/tmp/E15.5_Liver_Ribo_translating_ORFs.tsv",
            "_translating_ORFs.tsv",
        )
        self.assertNotEqual(brain, liver)

    def test_preserves_dots_in_ribocode_gtf_sample_names(self):
        sample_id = MODULE.infer_sample_id_from_prediction_path(
            "/tmp/E15.5_Brain_Ribo.gtf",
            ".gtf",
        )
        self.assertEqual(sample_id, "E15.5_Brain_Ribo")


class ParseRiboCodeTests(unittest.TestCase):
    def test_parses_ribocode_gtf_with_sidecar_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            reference_gtf = tmpdir / "reference.gtf"
            ribocode_gtf = tmpdir / "sample.gtf"
            ribocode_txt = tmpdir / "sample.txt"

            reference_gtf.write_text(textwrap.dedent("""\
                chr1\ttest\tgene\t100\t220\t.\t+\t.\tgene_id "gene1"; gene_name "Gene1";
                chr1\ttest\texon\t100\t220\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; gene_name "Gene1";
                chr1\ttest\tCDS\t130\t210\t.\t+\t0\tgene_id "gene1"; transcript_id "tx1"; gene_name "Gene1";
            """))
            ribocode_gtf.write_text(textwrap.dedent("""\
                chr1\tRiboCode\tORF\t130\t210\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; orf_id "orf1";
                chr1\tRiboCode\texon\t130\t150\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; orf_id "orf1";
                chr1\tRiboCode\texon\t190\t210\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1"; orf_id "orf1";
            """))
            ribocode_txt.write_text(
                "ORF_ID\tadjusted_pval\tpval_combined\n"
                "orf1\t1e-04\t1e-03\n"
            )

            gtf_index = MODULE.GTFIndex(str(reference_gtf))
            candidates = MODULE.parse_ribocode(str(ribocode_gtf), gtf_index, "sample", min_len=0)

        self.assertEqual(len(candidates), 1)
        candidate = candidates[0]
        self.assertEqual(candidate.chrom, "chr1")
        self.assertEqual(candidate.blocks, ((130, 150), (190, 210)))
        self.assertEqual(candidate.tid, "tx1")
        self.assertEqual(candidate.gid, "gene1")
        self.assertIn(("RiboCode", "sample"), candidate.sources)
        self.assertEqual(candidate.tool_pvalues["RiboCode"], 1e-04)
        self.assertAlmostEqual(candidate.tool_scores["RiboCode"], 4.0)


if __name__ == "__main__":
    unittest.main()
