# ----------------------------------------------------------------------------
# Copyright (c) 2024, Thanh Le Viet, Core Bioinformatics, Quadram Institute Bioscience
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import DNAFASTAFormat
from ._utils import run_command, validate_params

import shlex

from ._format import USEARCHTextFile


def unoise3(
    sequences: DNAFASTAFormat, minsize: int = 8, unoise_alpha: float = 2.0
) -> (DNAFASTAFormat, USEARCHTextFile):  # type: ignore

    validate_params([minsize, unoise_alpha])

    zotus_seqs = DNAFASTAFormat()
    tabbed_out = USEARCHTextFile()

    _minsize = f"-minsize {minsize}"
    _unoise_alpha = f"-unoise_alpha {unoise_alpha}"
    _cmd = f"usearch -unoise3 {sequences} -zotus {zotus_seqs} -tabbedout {tabbed_out} {_minsize} {_unoise_alpha}".strip()
    cmd = shlex.split(_cmd)
    run_command(cmd)
    return zotus_seqs, tabbed_out
