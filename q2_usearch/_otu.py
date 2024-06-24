# ----------------------------------------------------------------------------
# Copyright (c) 2024, Thanh Le Viet, Core Bioinformatics, Quadram Institute Bioscience
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import sqlite3

import biom
import skbio
import pandas as pd

from qiime2 import Metadata
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_table import BIOMV210Format
from q2_types.per_sample_sequences import (
    QIIME1DemuxDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqManifestFormat,
    YamlFormat,
)

from ._utils import (
    run_command,
    _fasta_with_sizes,
    _error_on_nonoverlapping_ids,
    USearchError,
    _uc_to_sqlite,
    _collapse_f_from_sqlite,
    _fasta_from_sqlite,
)
import shlex

from ._format import USEARCHTextFile


def otutab(
    reads: SingleLanePerSampleSingleEndFastqDirFmt,
    otudb: DNAFASTAFormat,
    pct_identity: float = 0.97,
    out_denoised: bool = False,
    threads: int = 1,
) -> (
    USEARCHTextFile,
    BIOMV210Format,
    USEARCHTextFile,
    SingleLanePerSampleSingleEndFastqDirFmt,
):  # type: ignore
    tabbed_out = USEARCHTextFile()
    tabbed_biomv210 = BIOMV210Format()
    mapout = USEARCHTextFile()
    unmapped = SingleLanePerSampleSingleEndFastqDirFmt()

    _otu = "-zotus" if out_denoised else "-otus"
    _cmd = f"usearch -otutab {reads} {_otu} {otudb} -otutabout {tabbed_out} -biomout otutab.biom -mapout {mapout} -id {pct_identity} -nomatchedfq {unmapped} -sizeout -threads {threads}".strip()
    cmd = shlex.split(_cmd)
    run_command(cmd)
    with biom.util.biom_open("otutab.biom") as f:
        tabbed_biomv210 = BIOMV210Format(f.name)
    return mapout, tabbed_biomv210, tabbed_out, unmapped

# Not available in v12
# def otu_norm(otutable: USEARCHTextFile, sample_size: int = 1000) -> USEARCHTextFile:
#     tabbed_out = USEARCHTextFile()
#     _cmd = f"usearch -otutab_norm {otutable} -sample_size {sample_size} -output {tabbed_out}".strip()
#     run_command(shlex.split(_cmd))
#     return tabbed_out
