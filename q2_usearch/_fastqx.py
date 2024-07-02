# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
# Copyright (c) 2024, Thanh Le Viet, Core Bioinformatics, Quadram Institute Bioscience
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import gzip
import shutil
import tempfile
import pandas as pd
import shlex
import yaml
from typing import List

from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import (
    QIIME1DemuxDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqManifestFormat,
    YamlFormat,
)


from ._utils import run_command, validate_params
import logging as logger

logger.basicConfig(level=logger.DEBUG, format="%(asctime)s - %(message)s")

_mp_defaults = {
    "maxdiffs": 5,
    "pctid": 90,
    "nostagger": False,
    "minmergelen": 50,
    "maxmergelen": 270,
    "minqual": 0,
    "minovlen": 16,
    "trunctail": 2,
    "minlen": 64,
    "relabel": "@",
    "threads": 1,
}


def fastq_join():
    pass


def fastq_mergepairs(
    demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
    maxdiffs: int = _mp_defaults["maxdiffs"],
    pctid: int = _mp_defaults["pctid"],
    nostagger: bool = _mp_defaults["nostagger"],
    minmergelen: int = _mp_defaults["minmergelen"],
    maxmergelen: int = _mp_defaults["maxmergelen"],
    minqual: int = _mp_defaults["minqual"],
    minovlen: int = _mp_defaults["minovlen"],
    relabel: str = _mp_defaults["relabel"],
    threads: int = _mp_defaults["threads"],
) -> (SingleLanePerSampleSingleEndFastqDirFmt, SingleLanePerSamplePairedEndFastqDirFmt):  # type: ignore
    _, merged, unmerged = _merge_pairs_w_command_output(
        demultiplexed_seqs,
        maxdiffs,
        pctid,
        nostagger,
        minmergelen,
        maxmergelen,
        minqual,
        minovlen,
        threads,
    )
    return merged, unmerged


def fastx_uniques(
    sequences: QIIME1DemuxDirFmt, sizeout: bool = True, relabel: bool = True
) -> DNAFASTAFormat:  # type: ignore
    # TODO: the software crashes when uc/tableout is defined.
    # The magic converting from fastq to fasta happens here https://github.com/qiime2/q2-types/blob/70b511c9657e3b464d1b6c0ed18673a3f0990a48/q2_types/per_sample_sequences/_transformer.py#L188
    unique_sequences = DNAFASTAFormat()
    seqs_fp = f"{sequences}/seqs.fna"
    logger.debug(f"Seqs file: {seqs_fp}")
    _relabel = (
        "-relabel @" if relabel else ""
    )  # relabels the sequences with the sample name using the special symbol @
    _sizeout = "-sizeout" if sizeout else ""
    _cmd = f"usearch -fastx_uniques {seqs_fp} -fastaout {unique_sequences} {_sizeout} {_relabel}".strip()
    cmd = shlex.split(_cmd)
    logger.debug(f"Command: {cmd}")
    run_command(cmd)
    return unique_sequences


def fastx_truncate(
    unique_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
    trunclen: int = 200,
    stripleft: int = 0,
    stripright: int = 0,
    padlen: int = 200,
    relabel: bool = False,
) -> SingleLanePerSampleSingleEndFastqDirFmt:  # type: ignore
    validate_params([trunclen, stripleft, stripright, padlen])
    truncated_seqs = SingleLanePerSampleSingleEndFastqDirFmt()

    # Read the manifest file
    manifest = pd.read_csv(
        os.path.join(str(unique_seqs), unique_seqs.manifest.pathspec),
        header=0,
        comment="#",
    )
    # logger.info(f"MANIFEST: {manifest}")

    # Update file paths in manifest
    manifest["filename"] = manifest["filename"].apply(
        lambda x: os.path.join(str(unique_seqs), x)
    )

    # Create a new manifest for truncated sequences
    truncated_manifest = FastqManifestFormat()
    truncated_manifest_fh = truncated_manifest.open()
    _write_manifest_header(truncated_manifest_fh)

    # Initialize a list to store information about missing files
    missing_files = []

    # Process each sample
    for i, row in manifest.iterrows():
        sample_id = row["sample-id"]
        input_fp = row["filename"]

        # Create a temporary directory for this sample
        with tempfile.TemporaryDirectory() as temp_dir:
            # Check if input file is compressed
            if input_fp.endswith(".gz"):
                # Decompress the file
                uncompressed_fp = os.path.join(
                    temp_dir, f"{sample_id}_uncompressed.fastq"
                )
                with gzip.open(input_fp, "rb") as f_in:
                    with open(uncompressed_fp, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
            else:
                uncompressed_fp = input_fp

            gz_truncated_path, fq_truncated_path = _get_output_paths(
                truncated_seqs, sample_id, i, 1
            )

            # Generate output path
            # _output_file_name = f"{sample_id}_L001_R1_001.fastq.gz"
            _output_file_name = os.path.basename(fq_truncated_path)
            output_fp = os.path.join(temp_dir, _output_file_name)

            # Prepare USEARCH command arguments
            cmd = [
                "usearch",
                "-fastx_truncate",
                uncompressed_fp,
                "-fastqout",
                fq_truncated_path,
            ]

            if trunclen is not None:
                cmd.extend(["-trunclen", str(trunclen)])
            if stripleft is not None:
                cmd.extend(["-stripleft", str(stripleft)])
            if stripright is not None:
                cmd.extend(["-stripright", str(stripright)])
            if padlen is not None:
                cmd.extend(["-padlen", str(padlen)])
            if relabel:
                cmd.extend(["-relabel", f"{sample_id}_"])

            # Process reads
            run_command(cmd)
            # Check if the output file exists and is not empty
            if (
                os.path.exists(fq_truncated_path)
                and os.path.getsize(fq_truncated_path) > 0
            ):
                # Compress output file
                final_output_fp = os.path.join(str(truncated_seqs), _output_file_name)
                with open(fq_truncated_path, "rb") as f_in:
                    with gzip.open(gz_truncated_path, "wb") as f_out:
                        logger.info(f"Copying {output_fp} to {final_output_fp}")
                        shutil.copyfileobj(f_in, f_out)
                os.remove(fq_truncated_path)
                logger.info(f"Output file created: {final_output_fp}")
                logger.info(
                    f"{sample_id},{os.path.basename(final_output_fp)},forward\n"
                )
                # Update manifest
                # Bear in mind the fasqt file is copied into a gzipped file
                truncated_manifest_fh.write(
                    f"{sample_id},{os.path.basename(final_output_fp)}.gz,forward\n"
                )
            else:
                logger.info("No sample passed!")
                # Record missing file
                missing_files.append(
                    {
                        "sample_id": sample_id,
                        "input_file": input_fp,
                        "reason": "Output file not created or empty",
                    }
                )

    truncated_manifest_fh.close()
    truncated_seqs.manifest.write_data(truncated_manifest, FastqManifestFormat)
    logger.info(f"MANIFEST: {truncated_manifest.path}")
    with open(truncated_manifest.path, "r") as f:
        logger.info(f.read())
    # Copy metadata
    metadata = YamlFormat()
    with open(os.path.join(str(unique_seqs), unique_seqs.metadata.pathspec), "r") as f:
        metadata.path.write_text(f.read())
    truncated_seqs.metadata.write_data(metadata, YamlFormat)

    return truncated_seqs


def fastq_filter(
    input_seqs: SingleLanePerSampleSingleEndFastqDirFmt,
    fastq_truncqual: int = None,
    fastq_maxee: float = None,
    fastq_trunclen: int = None,
    fastq_minlen: int = None,
    fastq_stripleft: int = None,
    fastq_maxee_rate: float = None,
    fastq_maxns: int = None,
    relabel: bool = False,
    fastq_eeout: bool = False,
    sample: str = None,
    threads: int = 1,
) -> SingleLanePerSampleSingleEndFastqDirFmt:
    filtered_seqs = SingleLanePerSampleSingleEndFastqDirFmt()

    # Read the manifest file
    manifest = pd.read_csv(
        os.path.join(str(input_seqs), input_seqs.manifest.pathspec),
        header=0,
        comment="#",
    )

    # Update file paths in manifest
    manifest["filename"] = manifest["filename"].apply(
        lambda x: os.path.join(str(input_seqs), x)
    )

    # Create a new manifest for filtered sequences
    filtered_manifest = FastqManifestFormat()
    filtered_manifest_fh = filtered_manifest.open()
    _write_manifest_header(filtered_manifest_fh)

    # Process each sample
    for _, row in manifest.iterrows():
        sample_id = row["sample-id"]
        input_fp = row["filename"]

        # Generate output paths
        output_gz, output_fq = _get_output_paths(filtered_seqs, sample_id, 0, 1)
        discarded_gz, discarded_fq = _get_output_paths(
            filtered_seqs, f"{sample_id}_discarded", 0, 1
        )

        # Prepare USEARCH command arguments
        cmd_args = [
            "usearch",
            "-fastq_filter",
            input_fp,
            "-fastqout",
            output_fq,
            "-fastqout_discarded",
            discarded_fq,
            "-threads",
            str(threads),
        ]

        if fastq_truncqual is not None:
            cmd_args.extend(["-fastq_truncqual", str(fastq_truncqual)])
        if fastq_maxee is not None:
            cmd_args.extend(["-fastq_maxee", str(fastq_maxee)])
        if fastq_trunclen is not None:
            cmd_args.extend(["-fastq_trunclen", str(fastq_trunclen)])
        if fastq_minlen is not None:
            cmd_args.extend(["-fastq_minlen", str(fastq_minlen)])
        if fastq_stripleft is not None:
            cmd_args.extend(["-fastq_stripleft", str(fastq_stripleft)])
        if fastq_maxee_rate is not None:
            cmd_args.extend(["-fastq_maxee_rate", str(fastq_maxee_rate)])
        if fastq_maxns is not None:
            cmd_args.extend(["-fastq_maxns", str(fastq_maxns)])
        if relabel:
            cmd_args.extend(["-relabel", f"{sample_id}_"])
        if fastq_eeout:
            cmd_args.append("-fastq_eeout")
        if sample:
            cmd_args.extend(["-sample", sample])

        # Run USEARCH command
        run_command(cmd_args)

        # Compress output files
        run_command(["gzip", output_fq])
        run_command(["gzip", discarded_fq])

        # Update manifest
        filtered_manifest_fh.write(f"{sample_id},{output_gz.name},forward\n")

    filtered_manifest_fh.close()
    filtered_seqs.manifest.write_data(filtered_manifest, FastqManifestFormat)

    # Copy metadata
    metadata = YamlFormat()
    with open(os.path.join(str(input_seqs), input_seqs.metadata.pathspec), "r") as f:
        metadata.path.write_text(f.read())
    filtered_seqs.metadata.write_data(metadata, YamlFormat)

    return filtered_seqs


def fastx_get_sample_names():
    pass


def _merge_pairs_w_command_output(
    demultiplexed_seqs: SingleLanePerSamplePairedEndFastqDirFmt,
    maxdiffs: int = _mp_defaults["maxdiffs"],
    pctid: int = _mp_defaults["pctid"],
    nostagger: bool = _mp_defaults["nostagger"],
    minmergelen: int = _mp_defaults["minmergelen"],
    maxmergelen: int = _mp_defaults["maxmergelen"],
    minqual: int = _mp_defaults["minqual"],
    minovlen: int = _mp_defaults["minovlen"],
    relabel: str = _mp_defaults["relabel"],
    threads: int = _mp_defaults["threads"],
) -> (
    List[str],
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
):  # type: ignore
    # create formats
    """
    Merges paired-end reads from demultiplexed sequences using USEARCH.

    Parameters:
    - demultiplexed_seqs (SingleLanePerSamplePairedEndFastqDirFmt): Input sequences.
    - maxdiffs, pctid, minmergelen, maxmergelen, minqual, minovlen (int): Parameters controlling merging.
    - nostagger (bool): If True, excludes staggered reads from merging.
    - relabel (str): Character for relabeling merged reads.
    - threads (int): Number of threads for merging.

    Returns:
    - Tuple with the USEARCH command, directory format objects for merged and unmerged reads.

    Initializes output formats, merges reads per sample with USEARCH, handles outputs, and writes manifest files.
    """
    merged = SingleLanePerSampleSingleEndFastqDirFmt()
    unmerged = SingleLanePerSamplePairedEndFastqDirFmt()

    # create manifests
    merged_manifest = FastqManifestFormat()
    merged_manifest_fh = merged_manifest.open()
    unmerged_manifest = FastqManifestFormat()
    unmerged_manifest_fh = unmerged_manifest.open()

    # write manifest headers
    _write_manifest_header(merged_manifest_fh, add_warning=True)
    _write_manifest_header(unmerged_manifest_fh)

    logger.debug(
        f"Manifests demultiplexed_seqs: {demultiplexed_seqs.manifest.pathspec}"
    )
    # generate input reads iterable
    manifest = pd.read_csv(
        os.path.join(str(demultiplexed_seqs), demultiplexed_seqs.manifest.pathspec),
        header=0,
        comment="#",
    )

    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(demultiplexed_seqs), x)
    )

    phred_offset = yaml.load(
        open(
            os.path.join(str(demultiplexed_seqs), demultiplexed_seqs.metadata.pathspec)
        ),
        Loader=yaml.SafeLoader,
    )["phred-offset"]

    id_to_fps = manifest.pivot(
        index="sample-id", columns="direction", values="filename"
    )

    for i, (sample_id, (fwd_fp, rev_fp)) in enumerate(id_to_fps.iterrows()):
        # The barcode id and lane number are not relevant for either format.
        # We might ultimately want to use a dir format other than these which
        # doesn't care about this information.
        # The read number (direction) is only relevant for the unmerged reads.

        gz_merged_path, fq_merged_path = _get_output_paths(merged, sample_id, i, 1)
        logger.debug(f"fq merged path: {fq_merged_path}")
        gz_unmerged_fwd_path, fq_unmerged_fwd_path = _get_output_paths(
            unmerged, sample_id, i, 1
        )
        gz_unmerged_rev_path, fq_unmerged_rev_path = _get_output_paths(
            unmerged, sample_id, i, 2
        )
        # Input
        _cmd = f"usearch -fastq_mergepairs {fwd_fp} -reverse {rev_fp} -fastqout {fq_merged_path}"
        # Output
        _cmd += f" -fastqout_notmerged_fwd {fq_unmerged_fwd_path} -fastqout_notmerged_rev {fq_unmerged_rev_path}"
        # Options
        _relabel = f"-relabel {relabel}" if relabel else ""
        _cmd += f" -fastq_maxdiffs {maxdiffs} \
                   -fastq_pctid {pctid} \
                    {'-fastq_nostagger' if nostagger else ''} \
                    -fastq_minmergelen {minmergelen} \
                    -fastq_maxmergelen {maxmergelen} \
                    -fastq_minqual {minqual} \
                    -fastq_minovlen {minovlen} \
                    {_relabel} \
                    -threads {threads}"
        logger.debug(f"Command: {_cmd}")
        # TODO: add report option
        cmd = shlex.split(_cmd)

        run_command(cmd)

        run_command(
            ["gzip", fq_merged_path, fq_unmerged_fwd_path, fq_unmerged_rev_path]
        )

        merged_manifest_fh.write(f"{sample_id},{gz_merged_path.name},forward\n")
        unmerged_manifest_fh.write(f"{sample_id},{gz_unmerged_fwd_path.name},forward\n")
        unmerged_manifest_fh.write(f"{sample_id},{gz_unmerged_rev_path.name},reverse\n")

    merged_manifest_fh.close()
    unmerged_manifest_fh.close()
    merged.manifest.write_data(merged_manifest, FastqManifestFormat)
    unmerged.manifest.write_data(unmerged_manifest, FastqManifestFormat)

    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({"phred-offset": phred_offset}))
    merged.metadata.write_data(metadata, YamlFormat)
    unmerged.metadata.write_data(metadata, YamlFormat)

    return cmd, merged, unmerged


def _get_output_paths(
    format_: SingleLanePerSampleSingleEndFastqDirFmt,
    sample_id: str,
    barcode_id: int,
    direction: int,
) -> tuple[str, str]:
    """
    Generate output paths for the given format, sample ID, barcode ID, and direction.

    Parameters:
    - format_ (SingleLanePerSampleSingleEndFastqDirFmt): The format object to generate paths for.
    - sample_id (str): The sample ID for the output paths.
    - barcode_id (int): The barcode ID for the output paths.
    - direction (int): The direction for the output paths.

    Returns:
    - tuple: A tuple containing the path and the path stripped of '.gz'.

    """
    path = format_.sequences.path_maker(
        sample_id=sample_id, barcode_id=barcode_id, lane_number=1, read_number=direction
    )
    return path, str(path).strip(".gz")


def _write_manifest_header(manifest_fh: str, add_warning: bool = False) -> None:
    manifest_fh.write("sample-id,filename,direction\n")
    if add_warning:
        manifest_fh.write("")
        manifest_fh.write("# direction is not meaningful for joined reads\n")
