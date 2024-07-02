import os
import pandas as pd
import unittest
from unittest.mock import patch
from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt
from q2_usearch._fastqx import fastx_truncate
from q2_usearch._utils import USearchError

import logging as logger
import gzip
logger.basicConfig(level=logger.DEBUG)

class TestFastxTruncate(TestPluginBase):
    package = "q2_usearch.tests"

    def setUp(self):
        super().setUp()

        # Load the real test file
        sequences_fp = self.get_data_path("merged_sequences.qza")

        # Import the file as a QIIME 2 Artifact
        sequences_artifact = Artifact.load(sequences_fp)

        # Get the SingleLanePerSampleSingleEndFastqDirFmt representation
        self.sequences = sequences_artifact.view(
            SingleLanePerSampleSingleEndFastqDirFmt
        )

    # @patch("q2_usearch._fastqx.run_command")
    def test_fastx_truncate_default_params(self):
        # Test with default parameters
        truncated_seqs = fastx_truncate(self.sequences)
        # Check if the returned object is of the correct type
        self.assertIsInstance(truncated_seqs, SingleLanePerSampleSingleEndFastqDirFmt)
    
    def test_fastx_truncate_custom_params(self):
        # Test with custom parameters
        truncated_seqs = fastx_truncate(
            self.sequences,
            trunclen=100,
            stripleft=10,
            stripright=5,
            padlen=150,
            relabel=True,
        )

        # Check if the returned object is of the correct type
        self.assertIsInstance(truncated_seqs, SingleLanePerSampleSingleEndFastqDirFmt)

        # Check the content of the truncated sequences manifest
        manifest_path = os.path.join(
            str(truncated_seqs), truncated_seqs.manifest.pathspec
        )
        logger.debug(f"Manifest path: {manifest_path}")

        manifest = pd.read_csv(manifest_path, header=0, comment="#")
        self.assertGreater(len(manifest), 0, "Truncated sequences manifest is empty")
        logger.debug(f"Manifest content:\n{manifest}")
        # Check if at least one truncated sequence file exists and is not empty
        truncated_file = os.path.join(str(truncated_seqs), manifest.iloc[0]["filename"])
        logger.debug(f"Truncated file path: {truncated_file}")

        self.assertTrue(
            os.path.exists(truncated_file),
            f"Truncated sequence file does not exist: {truncated_file}",
        )
        self.assertGreater(
            os.path.getsize(truncated_file),
            0,
            f"Truncated sequence file is empty: {truncated_file}",
        )
        sample_id = manifest.iloc[0]["sample-id"]
        # Check if the sequences are actually truncated and processed according to the parameters
        try:
            with gzip.open(truncated_file, "rt") as f:
                for line in f:
                    if line.startswith("@"):  # FASTQ header
                        self.assertTrue(
                            sample_id in line,
                            f"Relabeling not applied: {line}",
                        )
                        continue
                    if line.startswith("+"):  # Quality score identifier
                        break
                    self.assertLessEqual(
                        len(line.strip()),
                        100,
                        f"Sequence is longer than specified truncation length: {line}",
                    )
        except gzip.BadGzipFile:
            logger.error(f"File is not a valid gzip file: {truncated_file}")
            raise
        except Exception as e:
            logger.error(f"Error reading file {truncated_file}: {str(e)}")
            raise

        logger.info("Test completed successfully")

    @patch("q2_usearch._fastqx.run_command")
    def test_fastx_truncate_error_handling(self, mock_run_command):
        # Simulate an error in the USEARCH command
        mock_run_command.side_effect = USearchError("USEARCH command failed")
        # Test with invalid parameters to trigger an error
        with self.assertRaises(USearchError):
            fastx_truncate(self.sequences)  # Invalid truncation length

    def test_fastx_truncate_output_content(self):
        # Run the function with real data
        truncated_seqs = fastx_truncate(self.sequences, trunclen=100)

        # Check the content of the truncated sequences manifest
        manifest = pd.read_csv(
            os.path.join(str(truncated_seqs), truncated_seqs.manifest.pathspec),
            header=0,
            comment="#",
        )
        self.assertGreater(len(manifest), 0, "Truncated sequences manifest is empty")

        # Check if at least one truncated sequence file exists and is not empty
        truncated_file = os.path.join(str(truncated_seqs), manifest.iloc[0]["filename"])
        self.assertTrue(
            os.path.exists(truncated_file), "Truncated sequence file does not exist"
        )
        self.assertGreater(
            os.path.getsize(truncated_file), 0, "Truncated sequence file is empty"
        )

        # Check if the sequences are actually truncated
        with gzip.open(truncated_file, "r") as f:
            for line in f:
                if line.startswith(b"@"):  # FASTQ header
                    continue
                if line.startswith(b"+"):  # Quality score identifier
                    break
                self.assertLessEqual(
                    len(line.strip()),
                    100,
                    "Sequence is longer than specified truncation length",
                )
    

if __name__ == "__main__":
    unittest.main()