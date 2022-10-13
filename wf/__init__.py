from dataclasses import dataclass
from typing import List, Union

from dataclasses_json import dataclass_json
from latch import small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from .assembly import assembly_wf
from .binning import binning_wf
from .docs import megs_DOCS
from .functional import FunctionalOutput, functional_wf
from .types import ProdigalOutput, Sample, fARGeneModel


@dataclass_json
@dataclass
class WfResults:
    binning_results: List[LatchDir]
    prodigal_results: List[LatchDir]
    macrel_results: List[LatchDir]
    fargene_results: List[LatchDir]
    gecco_results: List[LatchDir]


@small_task
def organize_final_outputs(
    functional_results: List[FunctionalOutput],
    binning_results: List[LatchDir],
) -> WfResults:

    return WfResults(
        binning_results=binning_results,
        prodigal_results=[func.prodigal_result for func in functional_results],
        macrel_results=[func.macrel_result for func in functional_results],
        fargene_results=[func.fargene_result for func in functional_results],
        gecco_results=[func.gecco_result for func in functional_results],
    )


@workflow(megs_DOCS)
def megs(
    samples: List[Sample],
    min_count: int = 2,
    k_min: int = 21,
    k_max: int = 141,
    k_step: int = 12,
    min_contig_len: int = 200,
    prodigal_output_format: ProdigalOutput = ProdigalOutput.gbk,
    fargene_hmm_model: fARGeneModel = fARGeneModel.class_a,
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

    # Functional
    functional_results = functional_wf(
        assembly_data=assembly_dirs,
        prodigal_output_format=prodigal_output_format,
        fargene_hmm_model=fargene_hmm_model,
    )

    organized_outputs = organize_final_outputs(
        functional_results=functional_results,
        binning_results=binning_results,
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
        "prodigal_output_format": ProdigalOutput.gff,
        "fargene_hmm_model": fARGeneModel.class_b_1_2,
    },
)
