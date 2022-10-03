from typing import List, Union

from latch import workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from .assembly import assembly_wf
from .binning import binning_wf
from .docs import metamage_DOCS
from .kaiju import kaiju_wf
from .types import HostData, ProdigalOutput, Sample, TaxonRank, fARGeneModel


@workflow(metamage_DOCS)
def metamage(
    sample: Sample,
    host_data: HostData,
    kaiju_ref_db: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    sample_name: str = "metamage_sample",
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
    kaiju2table, krona_plot = kaiju_wf(
        read_dir=unaligned,
        kaiju_ref_db=kaiju_ref_db,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
        sample_name=sample_name,
        taxon_rank=taxon_rank,
    )

    assembly_dir, metassembly_results = assembly_wf(
        read_dir=unaligned,
        sample_name=sample_name,
        min_count=min_count,
        k_min=k_min,
        k_max=k_max,
        k_step=k_step,
        min_contig_len=min_contig_len,
    )

    # Binning
    binning_results = binning_wf(
        read_dir=unaligned, assembly_dir=assembly_dir, sample_name=sample_name
    )

    return [
        kaiju2table,
        krona_plot,
        metassembly_results,
        binning_results,
    ]


LaunchPlan(
    metamage,  # workflow name
    "Example Metagenome (Crohn's disease gut microbiome)",  # name of test data
    {
        "sample": Sample(
            read1=LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
            read2=LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
        ),
        "host_data": HostData(
            host_genome=LatchFile(
                "s3://latch-public/test-data/4318/Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz"
            ),
            host_name="homo_sapiens",
        ),
        "kaiju_ref_db": LatchFile(
            "s3://latch-public/test-data/4318/kaiju_db_plasmids.fmi"
        ),
        "kaiju_ref_nodes": LatchFile("s3://latch-public/test-data/4318/nodes.dmp"),
        "kaiju_ref_names": LatchFile("s3://latch-public/test-data/4318/names.dmp"),
        "sample_name": "SRR579292",
        "taxon_rank": TaxonRank.species,
        "min_count": 2,
        "k_min": 21,
        "k_max": 141,
        "k_step": 12,
        "min_contig_len": 200,
        "prodigal_output_format": ProdigalOutput.gff,
        "fargene_hmm_model": fARGeneModel.class_b_1_2,
    },
)
