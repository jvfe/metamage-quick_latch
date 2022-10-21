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
    """Metagenomic assembly, binning and taxonomic classification

    metamage-quick
    ----------

    metamage-quick is a workflow for taxonomic classification, assembly
    and binning of short-read host-associated metagenomics datasets.

    It's composed of:

    ## Assembly

    - [MEGAHIT](https://github.com/voutcn/megahit) for assembly [^1]
    - [MetaQuast](https://github.com/ablab/quast) for assembly evaluation

    ## Binning

    - BowTie2 and [Samtools](https://github.com/samtools/samtools)[^11] to
      building depth files for binning.
    - [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) for
      binning [^2]

    ## Taxonomic classification of reads

    - [Kaiju](https://github.com/bioinformatics-centre/kaiju) for
      taxonomic classification [^3]
    - [KronaTools](https://github.com/marbl/Krona/wiki/KronaTools) for
      visualizing taxonomic classification results

    # Output tree

    - |megs
      - |{sample_name}
        - |kaiju
        - |MEGAHIT
        - |MetaQuast - Assembly evaluation report
        - |{sample_name}_assembly_idx - BowTie Index from assembly data
        - |{sample_name}_assembly_sorted.bam - Reads aligned to assembly contigs
        - |METABAT

    # Where to get the data?

    - Kaiju indexes can be generated based on a reference database but
      you can also find some pre-built ones in the sidebar of the
      [Kaiju website](https://kaiju.binf.ku.dk/server).
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
