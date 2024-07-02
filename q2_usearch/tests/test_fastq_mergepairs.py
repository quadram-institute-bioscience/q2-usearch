import os
import pandas as pd
import unittest
from unittest.mock import patch
from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
)
from q2_usearch._fastqx import fastq_mergepairs
from q2_usearch._utils import USearchError


class TestFastqMergepairs(TestPluginBase):
    package = "q2_usearch.tests"

    def setUp(self):
        super().setUp()

        # Load the real test file
        paired_end_sequences_fp = self.get_data_path("sub-sample.qza")

        # Import the file as a QIIME 2 Artifact
        paired_end_sequences_artifact = Artifact.load(paired_end_sequences_fp)

        # Get the SingleLanePerSamplePairedEndFastqDirFmt representation
        self.paired_end_sequences = paired_end_sequences_artifact.view(
            SingleLanePerSamplePairedEndFastqDirFmt
        )

    @patch("q2_usearch._fastqx.run_command")
    def test_fastq_mergepairs_default_params(self, mock_run_command):
        # Test with default parameters
        merged_seqs, unmerged_seqs = fastq_mergepairs(self.paired_end_sequences)

        # Check if run_command was called with the correct command
        for call in mock_run_command.call_args_list:
            cmd = call[0][0]
            if "usearch" in cmd[0]:
                print("Command:", " ".join(cmd))
                # This is the USEARCH command
                self.assertIn("-fastq_mergepairs", cmd)

        # Check if the returned objects are of the correct type
        self.assertIsInstance(merged_seqs, SingleLanePerSampleSingleEndFastqDirFmt)
        self.assertIsInstance(unmerged_seqs, SingleLanePerSamplePairedEndFastqDirFmt)

    @patch("q2_usearch._fastqx.run_command")
    def test_fastq_mergepairs_custom_params(self, mock_run_command):
        # Test with custom parameters
        merged_seqs, unmerged_seqs = fastq_mergepairs(
            self.paired_end_sequences,
            maxdiffs=10,
            pctid=85,
            minmergelen=100,
            maxmergelen=300,
            minqual=10,
            minovlen=20,
            threads=4,
        )

        # Check if run_command was called with the correct command
        for call in mock_run_command.call_args_list:
            cmd = call[0][0]
            if "usearch" in cmd[0]:
                print("Command:", " ".join(cmd))
                # This is the USEARCH command
                self.assertIn("-fastq_mergepairs", cmd)
                self.assertIn("-fastq_maxdiffs", cmd)
                self.assertIn("10", cmd)
                self.assertIn("-fastq_pctid", cmd)
                self.assertIn("85", cmd)
                self.assertIn("-fastq_minmergelen", cmd)
                self.assertIn("100", cmd)
                self.assertIn("-fastq_maxmergelen", cmd)
                self.assertIn("300", cmd)
                self.assertIn("-fastq_minqual", cmd)
                self.assertIn("10", cmd)
                self.assertIn("-fastq_minovlen", cmd)
                self.assertIn("20", cmd)
                self.assertIn("-threads", cmd)
                self.assertIn("4", cmd)

        # Check if the returned objects are of the correct type
        self.assertIsInstance(merged_seqs, SingleLanePerSampleSingleEndFastqDirFmt)
        self.assertIsInstance(unmerged_seqs, SingleLanePerSamplePairedEndFastqDirFmt)

    @patch("q2_usearch._fastqx.run_command")
    def test_fastq_mergepairs_error_handling(self, mock_run_command):
        # Simulate an error in the USEARCH command
        mock_run_command.side_effect = USearchError("USEARCH command failed")

        # Check if the function raises the USearchError
        with self.assertRaises(USearchError):
            fastq_mergepairs(self.paired_end_sequences)

    def test_fastq_mergepairs_output_content(self):
        # Run the function with real data (not mocked)
        merged_seqs, unmerged_seqs = fastq_mergepairs(self.paired_end_sequences)

        # Check the content of the merged sequences manifest
        merged_manifest = pd.read_csv(
            os.path.join(str(merged_seqs), merged_seqs.manifest.pathspec),
            header=0,
            comment="#",
        )
        self.assertGreater(
            len(merged_manifest), 0, "Merged sequences manifest is empty"
        )

        # Check the content of the unmerged sequences manifest
        unmerged_manifest = pd.read_csv(
            os.path.join(str(unmerged_seqs), unmerged_seqs.manifest.pathspec),
            header=0,
            comment="#",
        )
        self.assertGreater(
            len(unmerged_manifest), 0, "Unmerged sequences manifest is empty"
        )

        # Check if there are exactly 3 FASTQ files for each sample
        sample_ids = merged_manifest["sample-id"].unique()
        for sample_id in sample_ids:
            merged_files = merged_manifest[merged_manifest["sample-id"] == sample_id][
                "filename"
            ].tolist()
            unmerged_files = unmerged_manifest[
                unmerged_manifest["sample-id"] == sample_id
            ]["filename"].tolist()

            self.assertEqual(
                len(merged_files),
                1,
                f"Expected 1 merged file for sample {sample_id}, but found {len(merged_files)}",
            )
            self.assertEqual(
                len(unmerged_files),
                2,
                f"Expected 2 unmerged files for sample {sample_id}, but found {len(unmerged_files)}",
            )

            # Check if files exist and are not empty
            for file_path in merged_files + unmerged_files:
                full_path = os.path.join(
                    str(merged_seqs if file_path in merged_files else unmerged_seqs),
                    file_path,
                )
                self.assertTrue(
                    os.path.exists(full_path), f"File does not exist: {full_path}"
                )
                self.assertGreater(
                    os.path.getsize(full_path), 0, f"File is empty: {full_path}"
                )

        # Additional checks for file format
        for file_path in (
            merged_manifest["filename"].tolist()
            + unmerged_manifest["filename"].tolist()
        ):
            full_path = os.path.join(
                str(
                    merged_seqs
                    if file_path in merged_manifest["filename"].tolist()
                    else unmerged_seqs
                ),
                file_path,
            )
            self.assertTrue(
                full_path.endswith(".fastq.gz"),
                f"File is not in .fastq.gz format: {full_path}",
            )

        # Check total number of files
        total_files = len(merged_manifest) + len(unmerged_manifest)
        expected_total = len(sample_ids) * 3
        self.assertEqual(
            total_files,
            expected_total,
            f"Expected {expected_total} total files, but found {total_files}",
        )


if __name__ == "__main__":
    unittest.main()
