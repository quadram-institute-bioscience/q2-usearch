# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
# Copyright (c) 2024, Thanh Le Viet, Core Bioinformatics, Quadram Institute Bioscience
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import subprocess
import filecmp
import unittest
from unittest.mock import patch
from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import QIIME1DemuxDirFmt
from q2_usearch._fastqx import fastx_uniques
from q2_usearch._utils import USearchError

class TestFastxUniques(TestPluginBase):
    package = "q2_usearch.tests"

    def setUp(self):
        super().setUp()

        # Load the real test file
        sequences_fp = self.get_data_path("merged_sequences.qza")

        # Import the file as a QIIME 2 Artifact
        sequences_artifact = Artifact.load(sequences_fp)

        # Get the QIIME1DemuxDirFmt representation
        self.sequences = sequences_artifact.view(QIIME1DemuxDirFmt)

    @patch("q2_usearch._fastqx.run_command")
    def test_fastx_uniques_default_params(self, mock_run_command):
        # Test with default parameters
        unique_sequences = fastx_uniques(self.sequences)

        # Check if run_command was called with the correct command
        mock_run_command.assert_called_once()
        cmd = mock_run_command.call_args[0][0]
        print("Command:", " ".join(cmd))
        self.assertIn("usearch", cmd[0])
        self.assertIn("-fastx_uniques", cmd)
        self.assertIn("-sizeout", cmd)
        self.assertIn("-relabel", cmd)

        # Check if the returned object is of the correct type
        self.assertIsInstance(unique_sequences, DNAFASTAFormat)

    @patch("q2_usearch._fastqx.run_command")
    def test_fastx_uniques_custom_params(self, mock_run_command):
        # Test with custom parameters
        unique_sequences = fastx_uniques(self.sequences, sizeout=False, relabel=False)

        # Check if run_command was called with the correct command
        mock_run_command.assert_called_once()
        cmd = mock_run_command.call_args[0][0]
        print("Command:", " ".join(cmd))
        self.assertIn("usearch", cmd[0])
        self.assertIn("-fastx_uniques", cmd)
        self.assertNotIn("-sizeout", cmd)
        self.assertNotIn("-relabel", cmd)

        # Check if the returned object is of the correct type
        self.assertIsInstance(unique_sequences, DNAFASTAFormat)

    @patch("q2_usearch._fastqx.run_command")
    def test_fastx_uniques_error_handling(self, mock_run_command):
        # Simulate an error in the USEARCH command
        mock_run_command.side_effect = USearchError("USEARCH command failed")

        # Check if the function raises the USearchError
        with self.assertRaises(USearchError):
            fastx_uniques(self.sequences)

    def test_fastx_uniques_output_content(self):
        # Run the function with real data (not mocked)
        unique_sequences = fastx_uniques(self.sequences)

        # Check if the output file exists and is not empty
        self.assertTrue(os.path.exists(str(unique_sequences)))
        self.assertGreater(os.path.getsize(str(unique_sequences)), 0, "Output file is empty")

        # Check the content of the output file
        with open(str(unique_sequences), 'r') as f:
            content = f.read()
            self.assertIn(">", content, "Output file does not contain FASTA format sequences")

        # Count the number of sequences
        sequence_count = content.count(">")
        self.assertGreater(sequence_count, 0, "No sequences found in the output file")
        print(f"Number of unique sequences: {sequence_count}")

    def test_fastx_uniques_vs_usearch(self):
        # Run fastx_uniques
        unique_sequences = fastx_uniques(self.sequences)

        # Run USEARCH directly
        with tempfile.TemporaryDirectory() as temp_dir:
            usearch_output = os.path.join(temp_dir, 'usearch_uniques.fasta')
            input_file = os.path.join(str(self.sequences), 'seqs.fna')

            # Run USEARCH command
            usearch_cmd = [
                "usearch",
                "-fastx_uniques", input_file,
                "-fastaout", usearch_output,
                "-sizeout",
                "-relabel", "@"
            ]
            subprocess.run(usearch_cmd, check=True)

            # Compare the outputs
            self.assertTrue(filecmp.cmp(str(unique_sequences), usearch_output),
                            "Outputs differ between fastx_uniques and direct USEARCH command")

if __name__ == "__main__":
    unittest.main()