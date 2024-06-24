# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
# Copyright (c) 2024, Thanh Le Viet, Core Bioinformatics, Quadram Institute Bioscience
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin
from qiime2.plugin import Citations, Plugin, Metadata, SemanticType
from q2_usearch import __version__
from q2_usearch._methods import duplicate_table
import q2_usearch._fastqx
import q2_usearch._cluster

# from q2_usearch._type import UchimeStats
from q2_usearch._format import UchimeStatsFmt, UchimeStatsDirFmt, USEARCHTextFile, USEARCHDirFmt
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    Sequences,
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
    JoinedSequencesWithQuality,
)


citations = Citations.load("citations.bib", package="q2_usearch")

plugin = Plugin(
    name="usearch",
    version=__version__,
    website="https://github.com/quadram-institute-bioscience/q2-usearch",
    package="q2_usearch",
    description="Plugin for clustering and dereplicating with usearch(12). Analysing steps with usearch should be following this recommendation https://drive5.com/usearch/manual/uparse_pipeline.html.",
    short_description="Plugin for clustering and dereplicating with usearch(12).",
    # Please retain the plugin-level citation of 'Caporaso-Bolyen-2024'
    # as attribution of the use of this template, in addition to any citations
    # you add.
    citations=[citations["Edgar2018"], citations["Caporaso-Bolyen-2024"]],
)

UsearchText = SemanticType("UsearchText")
plugin.register_formats(USEARCHTextFile, USEARCHDirFmt)
plugin.register_semantic_types(UsearchText)
plugin.register_semantic_type_to_format(UsearchText, artifact_format=USEARCHDirFmt)


plugin.methods.register_function(
    function=q2_usearch._fastqx.fastq_mergepairs,
    inputs={"demultiplexed_seqs": SampleData[PairedEndSequencesWithQuality]},
    parameters={
        "maxdiffs": qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        "pctid": qiime2.plugin.Int % qiime2.plugin.Range(0, 100),
        "nostagger": qiime2.plugin.Bool,
        "minmergelen": qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        "maxmergelen": qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        "minqual": qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        "minovlen": qiime2.plugin.Int % qiime2.plugin.Range(5, None),
        "relabel": qiime2.plugin.Str,
        "threads": qiime2.plugin.Int,
    },
    outputs=[
        ("merged_sequences", SampleData[JoinedSequencesWithQuality]),
        ("unmerged_sequences", SampleData[PairedEndSequencesWithQuality]),
    ],
    input_descriptions={
        "demultiplexed_seqs": (
            "The demultiplexed paired-end sequences to " "be merged."
        ),
    },
    parameter_descriptions={
        "maxdiffs": (
            "Maximum number of mismatches in the area of overlap " "during merging."
        ),
        "pctid": (
            "Minimum percentage of identity for matches in the area of overlap "
            "during merging."
        ),
        "nostagger": ("If set to True, disallow merging of staggered read pairs."),
        "minmergelen": ("Minimum length of the merged read to be retained."),
        "maxmergelen": ("Maximum length of the merged read to be retained."),
        "minqual": ("Minimum quality score to keep a base during merging."),
        "minovlen": (
            "Minimum length of the area of overlap between reads " "during merging."
        ),
        "threads": (
            "The number of threads to use for computation. Does "
            "not scale much past 4 threads."
        ),
    },
    output_descriptions={
        "merged_sequences": "The merged sequences.",
        "unmerged_sequences": "The unmerged paired-end reads.",
    },
    name="Merge paired-end reads.",
    description=(
        "Merge paired-end sequence reads using usearch's "
        "fastq_mergepairs function. See the usearch documentation for "
        "details on how paired-end merging is performed, and for "
        "more information on the parameters to this method."
    ),
)

plugin.methods.register_function(
    function=q2_usearch._fastqx.fastx_truncate,
    inputs={"demultiplexed_seqs": SequencesWithQuality},
    parameters={
        "trunclen": qiime2.plugin.Int % qiime2.plugin.Range(1, None),
        "stripleft": qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        "stripright": qiime2.plugin.Int % qiime2.plugin.Range(0, None),
        "padlen": qiime2.plugin.Int % qiime2.plugin.Range(1, None),
        "relabel": qiime2.plugin.Bool,
        "threads": qiime2.plugin.Int % qiime2.plugin.Range(1, None),
    },
    outputs=[("truncated_seqs", SequencesWithQuality)],
    input_descriptions={
        "demultiplexed_seqs": "The single-end demultiplexed sequences to be truncated. Input sequences should be the output of mergepairs step."
    },
    parameter_descriptions={
        "trunclen": "Truncate sequences to the specified length.",
        "stripleft": "Strip the specified number of bases from the start of each sequence.",
        "stripright": "Strip the specified number of bases from the end of each sequence.",
        "padlen": "Pad sequences with N's to the specified length.",
        "relabel": "Relabel sequences with their sample ID.",
        "threads": "Number of threads to use for processing.",
    },
    output_descriptions={"truncated_seqs": "The resulting truncated sequences."},
    name="Truncate FASTQ sequences",
    description="This method truncates FASTQ sequences using the USEARCH fastx_truncate command.",
)

plugin.methods.register_function(
    function=q2_usearch._fastqx.fastx_uniques,
    inputs={
        "sequences": (
            SampleData[Sequences]
            | SampleData[SequencesWithQuality]
            | SampleData[JoinedSequencesWithQuality]
        )
    },
    parameters={
        "sizeout": qiime2.plugin.Bool,
        "relabel": qiime2.plugin.Bool,
    },
    outputs=[
        ("unique_sequences", FeatureData[Sequence]),
    ],
    input_descriptions={
        "sequences": "The sequences to be dereplicated.",
    },
    parameter_descriptions={
        "sizeout": ("Size annotations should be added to the output sequence labels."),
        "relabel": ("Re-label the dereplicated sequences."),
    },
    output_descriptions={
        "unique_sequences": "The dereplicated sequences.",
    },
    name="Dereplicate sequences.",
    description=(
        "Find the set of unique sequences in an input file, also called dereplication. "
        "Input is a FASTA or FASTQ file. Sequences are compared letter-by-letter and "
        "must be identical over the full length of both sequences (substrings do not "
        "match). Case is ignored, so an upper-case letter matches a lower-case letter."
    ),
)

plugin.methods.register_function(
    function=q2_usearch._cluster.cluster_otus,
    inputs={"sequences": FeatureData[Sequence]},
    parameters={
        "minsize": qiime2.plugin.Int,
        "relabel": qiime2.plugin.Str,
        "threads": qiime2.plugin.Int,
    },
    outputs=[("otus", FeatureData[Sequence]), ("uparse_out", UsearchText)],
    input_descriptions={"sequences": "The sequences to be clustered into OTUs."},
    parameter_descriptions={
        "minsize": "The minimum size for OTUs.",
        "relabel": "The prefix for OTU labels.",
        "threads": "The number of threads to use for clustering.",
    },
    output_descriptions={
        "otus": "The resulting OTU sequences.",
        "uparse_out": "UPARSE output information as metadata.",
    },
    name="Cluster OTUs using UPARSE",
    description="Clusters sequences into OTUs using the UPARSE algorithm.",
)