from dataclasses import dataclass
from typing import List, Union

from dataclasses_json import dataclass_json
from latch import small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from .assembly import AssemblyOut, assembly_wf
from .binning import binning_wf
from .docs import megs_DOCS
from .kaiju import kaiju_wf
from .types import Sample, TaxonRank


@dataclass_json
@dataclass
class WfResults:
    assembly_results: List[LatchDir]
    binning_results: List[LatchDir]
    krona_plots: List[LatchFile]
    kaiju2table_outs: List[LatchFile]


@small_task
def organize_final_outputs(
    assembly_results: List[AssemblyOut],
    binning_results: List[LatchDir],
    krona_plots: List[LatchFile],
    kaiju2table_outs: List[LatchFile],
) -> WfResults:

    metaquast_results = [
        assembly_result.evaluation for assembly_result in assembly_results
    ]

    return WfResults(
        assembly_results=metaquast_results,
        binning_results=binning_results,
        krona_plots=krona_plots,
        kaiju2table_outs=kaiju2table_outs,
    )


@workflow(megs_DOCS)
def megs(
    samples: List[Sample],
    kaiju_ref_db: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    taxon_rank: TaxonRank = TaxonRank.species,
    min_count: int = 2,
    k_min: int = 21,
    k_max: int = 141,
    k_step: int = 12,
    min_contig_len: int = 200,
) -> WfResults:
    """Metagenomic assembly with MEGAHit
    megs
    ----------
    megs assembles metagenomic reads with MEGAHit.
    """

    assembly_dirs = assembly_wf(
        samples=samples,
        min_count=min_count,
        k_min=k_min,
        k_max=k_max,
        k_step=k_step,
        min_contig_len=min_contig_len,
    )

    # Binning
    binning_results = binning_wf(samples=samples, megahit_out=assembly_dirs)

    kaiju2table_outs, krona_plots = kaiju_wf(
        samples=samples,
        kaiju_ref_db=kaiju_ref_db,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
        taxon_rank=taxon_rank,
    )

    organized_outputs = organize_final_outputs(
        assembly_results=assembly_dirs,
        binning_results=binning_results,
        krona_plots=krona_plots,
        kaiju2table_outs=kaiju2table_outs,
    )

    return organized_outputs


LaunchPlan(
    megs,  # workflow name
    "Example Metagenome (Crohn's disease gut microbiome)",  # name of test data
    {
        "samples": [
            Sample(
                sample_name="SRR579291",
                read1=LatchFile("s3://latch-public/test-data/4318/SRR579291_1.fastq"),
                read2=LatchFile("s3://latch-public/test-data/4318/SRR579291_2.fastq"),
            ),
            Sample(
                sample_name="SRR579292",
                read1=LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
                read2=LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
            ),
        ],
        "min_count": 2,
        "k_min": 21,
        "k_max": 141,
        "k_step": 12,
        "min_contig_len": 200,
        "kaiju_ref_db": LatchFile(
            "s3://latch-public/test-data/4318/kaiju_db_viruses.fmi"
        ),
        "kaiju_ref_nodes": LatchFile(
            "s3://latch-public/test-data/4318/virus_nodes.dmp"
        ),
        "kaiju_ref_names": LatchFile(
            "s3://latch-public/test-data/4318/virus_names.dmp"
        ),
        "taxon_rank": TaxonRank.species,
    },
)
