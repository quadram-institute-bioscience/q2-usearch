from qiime2.plugin.testing import TestPluginBase
from unittest.mock import patch
import os
from qiime2 import Artifact
from q2_types.feature_data import DNAFASTAFormat, FeatureData, Sequence
from q2_usearch._cluster import cluster_otus
from q2_usearch._utils import USearchError
from q2_usearch._format import USEARCHTextFile


class TestClusterOtus(TestPluginBase):
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

    @patch("q2_usearch._cluster.run_command")
    def test_cluster_otus_default_params(self, mock_run_command):
        # Test with default parameters
        otus_seqs, uparse_out = cluster_otus(self.sequences)

        # Check if run_command was called with the correct command
        expected_cmd = [
            "usearch",
            "-cluster_otus",
            str(self.sequences),
            "-otus",
            str(otus_seqs),
            "-uparseout",
            str(uparse_out),
            "-relabel",
            "Otu",
            "-minsize",
            "2",
            "-threads",
            "1",
        ]
        mock_run_command.assert_called_once_with(expected_cmd)

        # Check if the returned objects are of the correct type
        self.assertIsInstance(otus_seqs, DNAFASTAFormat)
        self.assertIsInstance(uparse_out, USEARCHTextFile)

    @patch("q2_usearch._cluster.run_command")
    def test_cluster_otus_custom_params(self, mock_run_command):
        # Test with custom parameters
        otus_seqs, uparse_out = cluster_otus(
            self.sequences, minsize=5, relabel="CustomOtu", threads=4
        )

        # Check if run_command was called with the correct command
        expected_cmd = [
            "usearch",
            "-cluster_otus",
            str(self.sequences),
            "-otus",
            str(otus_seqs),
            "-uparseout",
            str(uparse_out),
            "-relabel",
            "CustomOtu",
            "-minsize",
            "5",
            "-threads",
            "4",
        ]
        mock_run_command.assert_called_once_with(expected_cmd)

    @patch("q2_usearch._cluster.run_command")
    def test_cluster_otus_error_handling(self, mock_run_command):
        # Simulate an error in the USEARCH command
        mock_run_command.side_effect = USearchError("USEARCH command failed")

        # Check if the function raises the USearchError
        with self.assertRaises(USearchError):
            cluster_otus(self.sequences)