#!/usr/bin/env python3

import importlib.util
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


if __name__ == "__main__":
    unittest.main()
