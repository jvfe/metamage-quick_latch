from latch.types import LatchAuthor, LatchMetadata, LatchParameter

megs_DOCS = LatchMetadata(
    display_name="megs",
    documentation="https://github.com/jvfe/megs_latch/blob/main/README.md",
    author=LatchAuthor(
        name="jvfe",
        github="https://github.com/jvfe",
    ),
    repository="https://github.com/jvfe/megs_latch",
    license="MIT",
    tags=["NGS", "metagenomics", "MAG"],
)

megs_DOCS.parameters = {
    "samples": LatchParameter(
        display_name="Sample data",
        description="Paired-end FASTQ files",
        batch_table_column=True,
    ),
    "k_min": LatchParameter(
        display_name="Minimum kmer size",
        description="Must be odd and <=255",
        section_title="MEGAHIT parameters",
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
}
