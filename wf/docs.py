from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchAuthor,
    LatchMetadata,
    LatchParameter,
    LatchRule,
    Params,
    Section,
    Spoiler,
    Text,
)

dmp_rule = LatchRule(regex=".dmp$", message="Must be a .dmp file")

PARAMS = {
    "samples": LatchParameter(
        display_name="Sample data",
        description="Paired-end FASTQ files",
        batch_table_column=True,
    ),
    "k_min": LatchParameter(
        display_name="Minimum kmer size",
        description="Must be odd and <=255",
    ),
    "k_max": LatchParameter(
        display_name="Maximum kmer size",
        description="Must be odd and <=255",
    ),
    "k_step": LatchParameter(
        display_name="Increment of kmer size of each iteration",
        description="Must be even and <=28",
    ),
    "min_count": LatchParameter(
        display_name="Minimum multiplicity for filtering (k_min+1)-mers",
    ),
    "min_contig_len": LatchParameter(
        display_name="Minimum length of contigs to output",
    ),
    "kaiju_ref_db": LatchParameter(
        display_name="Kaiju reference database (FM-index)",
        description="Kaiju reference database '.fmi' file.",
        detail="(.fmi)",
        rules=[LatchRule(regex=".fmi", message="Must be an .fmi index file")],
    ),
    "kaiju_ref_nodes": LatchParameter(
        display_name="Kaiju reference database nodes",
        description="Kaiju reference nodes, 'nodes.dmp' file.",
        detail="(.dmp)",
        rules=[dmp_rule],
    ),
    "kaiju_ref_names": LatchParameter(
        display_name="Kaiju reference database names",
        description="Kaiju reference taxon names, 'names.dmp' file.",
        detail="(.dmp)",
        rules=[dmp_rule],
    ),
    "taxon_rank": LatchParameter(
        display_name="Taxonomic rank (kaiju2table)",
        description="Taxonomic rank for summary table output (kaiju2table).",
    ),
    "prodigal_output_format": LatchParameter(
        display_name="Prodigal output file format",
        description="Specify main output file format (one of gbk, gff or sco).",
    ),
    "fargene_hmm_model": LatchParameter(
        display_name="fARGene's HMM model",
        description="The Hidden Markov Model that should be used to predict ARGs from the data",
    ),
}

FLOW = [
    Section(
        "Samples",
        Text(
            "Sample provided has to include an identifier for the sample (Sample name)"
            " and two files corresponding to the reads (paired-end)"
        ),
        Params("samples"),
    ),
    Section(
        "Assembly parameters",
        Text("Parameters for the assembly software MEGAHIT"),
        Params("k_min", "k_max", "k_step", "min_count", "min_contig_len"),
    ),
    Section(
        "Taxonomic classification",
        Text(
            "Parameters for the taxonomic classifier Kaiju"
            " choose which database to use and at which taxonomic"
            " level the final TSV report should be generated"
        ),
        Params("kaiju_ref_db", "kaiju_ref_nodes", "kaiju_ref_names", "taxon_rank"),
    ),
    Section(
        "Functional annotation parameters",
        Text("Options for the functional annotation subworkflow"),
        Params("prodigal_output_format", "fargene_hmm_model"),
    ),
]

metamage_DOCS = LatchMetadata(
    display_name="MetaMage-quick",
    documentation="https://github.com/jvfe/metamage-quick_latch/blob/main/README.md",
    author=LatchAuthor(
        name="jvfe",
        github="https://github.com/jvfe",
    ),
    repository="https://github.com/jvfe/metamage-quick_latch",
    license="MIT",
    parameters=PARAMS,
    tags=["NGS", "metagenomics", "MAG"],
    flow=FLOW,
)
