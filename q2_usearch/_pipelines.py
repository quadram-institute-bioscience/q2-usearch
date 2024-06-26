import qiime2.plugin
from q2_types.feature_data import FeatureData, Sequence
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import SampleData
from ._format import USEARCHDirFmt


# Define the pipeline function
def denoise_pipeline(
    ctx,
    sequences,
    minsize=2,
    relabel="Otu",
    unoise_alpha=2.0,
    pct_identity=0.97,
    out_denoised=True,
    threads=1,
):
    # Step 1: fastx_uniques
    unique_seqs = ctx.get_action("usearch", "fastx_uniques")
    unique_output = unique_seqs(sequences, sizeout=True, relabel=True)

    # Step 2: cluster_otus
    cluster = ctx.get_action("usearch", "cluster_otus")
    cluster_output = cluster(
        unique_output.unique_seqs, minsize=minsize, relabel=relabel, threads=threads
    )

    # Step 3: unoise3
    unoise = ctx.get_action("usearch", "unoise3")
    unoise_output = unoise(
        cluster_output.clustered_seqs,
        minsize=minsize,
        unoise_alpha=unoise_alpha,
        threads=threads,
    )

    # Step 4: otutab
    otutab = ctx.get_action("usearch", "otutab")
    otutab_output = otutab(
        sequences,
        unoise_output.zotus,
        pct_identity=pct_identity,
        out_denoised=out_denoised,
        threads=threads,
    )

    # Return all outputs
    return (
        otutab_output.feature_table,
        otutab_output.feature_sequences,
        otutab_output.otu_table,
        otutab_output.denoised_sequences if out_denoised else None
    )