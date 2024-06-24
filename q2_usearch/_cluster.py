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
from ._utils import run_command, _fasta_with_sizes, _error_on_nonoverlapping_ids, USearchError, _uc_to_sqlite, _collapse_f_from_sqlite, _fasta_from_sqlite
import shlex

from ._format import USEARCHTextFile

def cluster_fast(sequences: DNAFASTAFormat, 
                 table: biom.Table,
                 perc_identity: float = 0.9,
                 strand: str = 'plus',
                 threads: int = 1
                ) -> (biom.Table, DNAFASTAFormat): # type: ignore
    clustered_sequences = DNAFASTAFormat()
    with tempfile.NamedTemporaryFile() as fasta_with_sizes:
        with tempfile.NamedTemporaryFile() as out_uc:
            _fasta_with_sizes(str(sequences), fasta_with_sizes.name, table)
            cmd = ['usearch',
                   '--cluster_size', fasta_with_sizes.name,
                   '--id', str(perc_identity),
                   '--centroids', str(clustered_sequences),
                   '--uc', out_uc.name,
                   '--strand', str(strand),
                   '--qmask', 'none',  # ensures no lowercase DNA chars
                   '--xsize',
                   '--threads', str(threads),
                   '--minseqlength', '1',
                   '--fasta_width', '0']
            run_command(cmd)
            out_uc.seek(0)

            conn = _uc_to_sqlite(out_uc)
            collapse_f = _collapse_f_from_sqlite(conn)

    table = table.collapse(collapse_f, norm=False, min_group_size=1,
                           axis='observation',
                           include_collapsed_metadata=False)

    return table, clustered_sequences


def cluster_otus(
    sequences: DNAFASTAFormat,
    minsize: int = 2,
    relabel: str = "Otu",
    threads: int = 1,
) -> (DNAFASTAFormat, USEARCHTextFile):  # type: ignore
    otus_seqs = DNAFASTAFormat()
    uparse_out = USEARCHTextFile()

    _relabel = f"-relabel {relabel}"
    _minsize = f"-minsize {minsize}"
    _cmd = f"usearch -cluster_otus {sequences} -otus {otus_seqs} -uparseout {uparse_out} {_relabel} {_minsize} -threads {threads}".strip()
    cmd = shlex.split(_cmd)
    run_command(cmd)
    return otus_seqs, uparse_out


def cluster_smallmem(
    ctx, sequences, table, reference_sequences, perc_identity, strand="plus", threads=1
):

    cluster_features_closed_reference = ctx.get_action(
        'vsearch', 'cluster_features_closed_reference')
    filter_features = ctx.get_action('feature_table', 'filter_features')
    cluster_features_de_novo = ctx.get_action(
        'vsearch', 'cluster_features_de_novo')
    merge = ctx.get_action('feature_table', 'merge')
    merge_seqs = ctx.get_action('feature_table', 'merge_seqs')

    skipped_closed_ref = True
    try:
        closed_ref_table, rep_seqs, unmatched_seqs = \
            cluster_features_closed_reference(
                sequences=sequences, table=table,
                reference_sequences=reference_sequences,
                perc_identity=perc_identity,
                strand=strand, threads=threads)
        skipped_closed_ref = False
    except USearchError:  # No matches
        pass

    # If cluster_features_closed_reference fails to match, we need to
    # pass the source data into cluster_features_de_novo wholesale.
    if skipped_closed_ref:
        unmatched_seqs, closed_ref_table = sequences, table

    # It is possible that all of the sequences matched the reference database,
    # if that is the case, don't worry about running cluster_features_de_novo.
    if unmatched_seqs.view(pd.Series).size > 0:
        unmatched_seqs_md = unmatched_seqs.view(Metadata)
        unmatched_table, = filter_features(table=table,
                                           metadata=unmatched_seqs_md)

        de_novo_table, de_novo_seqs = cluster_features_de_novo(
            sequences=unmatched_seqs, table=unmatched_table,
            perc_identity=perc_identity, strand=strand, threads=threads)

        if skipped_closed_ref:
            merged_reference_seqs, = merge_seqs(data=[reference_sequences,
                                                      de_novo_seqs])
            outputs = (de_novo_table, de_novo_seqs, merged_reference_seqs)
        else:
            merged_table, = merge(
                tables=[closed_ref_table, de_novo_table],
                overlap_method='error_on_overlapping_feature')

            merged_rep_seqs, = merge_seqs(data=[rep_seqs, de_novo_seqs])

            merged_reference_seqs, = merge_seqs(data=[reference_sequences,
                                                      de_novo_seqs])
            outputs = (merged_table, merged_rep_seqs, merged_reference_seqs)
    else:  # skipped de novo
        outputs = (closed_ref_table, rep_seqs, reference_sequences)
    return outputs


def cluster_mt():
    pass
