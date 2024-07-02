# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
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


def cluster_otus(
    sequences: DNAFASTAFormat,
    minsize: int = 2,
    relabel: str = "Otu",
    threads: int = 1,
) -> (DNAFASTAFormat, USEARCHTextFile):  # type: ignore
    otus_seqs = DNAFASTAFormat()
    uparse_out = USEARCHTextFile()

    validate_params([minsize, threads])

    _relabel = f"-relabel {relabel}"
    _minsize = f"-minsize {minsize}"
    _cmd = f"usearch -cluster_otus {sequences} -otus {otus_seqs} -uparseout {uparse_out} {_relabel} {_minsize} -threads {threads}".strip()
    cmd = shlex.split(_cmd)
    run_command(cmd)
    return otus_seqs, uparse_out


def cluster_fast():
    pass


def cluster_smallmem():
    pass


def cluster_mt():
    pass
