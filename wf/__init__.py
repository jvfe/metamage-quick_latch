import itertools
from typing import List, Union

from latch import small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from .assembly import AssemblyOut, assembly_wf
from .binning import binning_wf
from .docs import maggie_DOCS
from .kaiju import kaiju_wf
from .types import Sample, TaxonRank


@small_task
def organize_final_outputs(
    kaiju_results: List[LatchFile],
    assembly_data: List[AssemblyOut],
    binning_results: List[LatchDir],
) -> List[Union[LatchFile, LatchDir]]:

    metassembly_results = [result.evaluation for result in assembly_data]
    all_outputs = [kaiju_results, metassembly_results, binning_results]

    return list(itertools.chain(*all_outputs))


@workflow(maggie_DOCS)
def maggie(
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
) -> List[Union[LatchFile, LatchDir]]:
    """Metagenomic assembly and binning

    maggie
    ----------
    """

    # Kaiju taxonomic classification
    krona_plots = kaiju_wf(
        samples=samples,
        kaiju_ref_db=kaiju_ref_db,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
        taxon_rank=taxon_rank,
    )

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

    organized_outputs = organize_final_outputs(
        kaiju_results=krona_plots,
        assembly_data=assembly_dirs,
        binning_results=binning_results,
    )

    return organized_outputs


LaunchPlan(
    maggie,  # workflow name
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
        "kaiju_ref_db": LatchFile(
            "s3://latch-public/test-data/4318/kaiju_db_nr_2021-02-24.fmi"
        ),
        "kaiju_ref_nodes": LatchFile("s3://latch-public/test-data/4318/nodes.dmp"),
        "kaiju_ref_names": LatchFile("s3://latch-public/test-data/4318/names.dmp"),
        "taxon_rank": TaxonRank.species,
        "min_count": 2,
        "k_min": 21,
        "k_max": 141,
        "k_step": 12,
        "min_contig_len": 200,
    },
)
