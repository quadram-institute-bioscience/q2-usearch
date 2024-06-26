from qiime2.plugin.testing import TestPluginBase
from unittest.mock import patch
import os
from qiime2 import Artifact
from q2_types.feature_data import DNAFASTAFormat, FeatureData, Sequence
from q2_usearch._unoise import unoise3
from q2_usearch._format import USEARCHTextFile
from q2_usearch._utils import USearchError


class TestUnoise3(TestPluginBase):
    package = "q2_usearch"

    def setUp(self):
        super().setUp()

        # Load the real test file
        test_dir = self.get_data_path(".")
        sequences_fp = os.path.join(test_dir, "uniques_5k.fa")

        # Import the file as a QIIME 2 Artifact
        sequences_artifact = Artifact.import_data(FeatureData[Sequence], sequences_fp)

        # Get the DNAFASTAFormat representation
        self.sequences = sequences_artifact.view(DNAFASTAFormat)

    @patch("q2_usearch._unoise.run_command")
    def test_unoise3_default_params(self, mock_run_command):
        # Test with default parameters
        zotus_seqs, tabbed_out = unoise3(self.sequences)

        # Check if run_command was called with the correct command
        expected_cmd = [
            "usearch",
            "-unoise3",
            str(self.sequences),
            "-zotus",
            str(zotus_seqs),
            "-tabbedout",
            str(tabbed_out),
            "-minsize",
            "8",
            "-unoise_alpha",
            "2.0",
        ]
        mock_run_command.assert_called_once_with(expected_cmd)

        # Check if the returned objects are of the correct type
        self.assertIsInstance(zotus_seqs, DNAFASTAFormat)
        self.assertIsInstance(tabbed_out, USEARCHTextFile)

    @patch("q2_usearch._unoise.run_command")
    def test_unoise3_custom_params(self, mock_run_command):
        # Test with custom parameters
        zotus_seqs, tabbed_out = unoise3(
            self.sequences, minsize=10, unoise_alpha=1.5
        )

        # Check if run_command was called with the correct command
        expected_cmd = [
            "usearch",
            "-unoise3",
            str(self.sequences),
            "-zotus",
            str(zotus_seqs),
            "-tabbedout",
            str(tabbed_out),
            "-minsize",
            "10",
            "-unoise_alpha",
            "1.5",
        ]
        mock_run_command.assert_called_once_with(expected_cmd)

    @patch("q2_usearch._unoise.run_command")
    def test_unoise3_error_handling(self, mock_run_command):
        # Simulate an error in the USEARCH command
        mock_run_command.side_effect = USearchError("USEARCH command failed")

        # Check if the function raises the USearchError
        with self.assertRaises(USearchError):
            unoise3(self.sequences)

    def test_unoise3_invalid_minsize(self):
        # Test with invalid minsize (negative value)
        with self.assertRaises(ValueError):
            unoise3(self.sequences, minsize=-1)

    def test_unoise3_invalid_unoise_alpha(self):
        # Test with invalid unoise_alpha (negative value)
        with self.assertRaises(ValueError):
            unoise3(self.sequences, unoise_alpha=-0.5)
