from latch.types import LatchDir, LatchFile

import wf.assembly
from wf.assembly import AssemblyOut
from wf.types import ProdigalOutput, Sample, TaxonRank, fARGeneModel

# wf.megs(
#     samples=[
#         Sample(
#             sample_name="SRR579291",
#             read1=LatchFile("s3://latch-public/test-data/4318/SRR579291_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/SRR579291_2.fastq"),
#         ),
#         Sample(
#             sample_name="SRR579292",
#             read1=LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
#         ),
#     ],
#     min_count=2,
#     k_min=21,
#     k_max=141,
#     k_step=12,
#     min_contig_len=200,
#     kaiju_ref_db=LatchFile("s3://latch-public/test-data/4318/kaiju_db_viruses.fmi"),
#     kaiju_ref_nodes=LatchFile("s3://latch-public/test-data/4318/virus_nodes.dmp"),
#     kaiju_ref_names=LatchFile("s3://latch-public/test-data/4318/virus_names.dmp"),
#     taxon_rank=TaxonRank.species,
#     prodigal_output_format=ProdigalOutput.gff,
#     fargene_hmm_model=fARGeneModel.class_b_1_2,
# )

wf.kaiju.kaiju_wf(
    samples=[
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
    kaiju_ref_db=LatchFile("s3://latch-public/test-data/4318/kaiju_db_viruses.fmi"),
    kaiju_ref_nodes=LatchFile("s3://latch-public/test-data/4318/virus_nodes.dmp"),
    kaiju_ref_names=LatchFile("s3://latch-public/test-data/4318/virus_names.dmp"),
    taxon_rank=TaxonRank.species,
)

# wf.binning.binning_wf(
#     samples=[
#         Sample(
#             sample_name="SRR579291",
#             read1=LatchFile("s3://latch-public/test-data/4318/SRR579291_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/SRR579291_2.fastq"),
#         ),
#         Sample(
#             sample_name="SRR579292",
#             read1=LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
#         ),
#     ],
#     megahit_out=[
#         AssemblyOut(
#             sample_name="SRR579291",
#             assembly_data=LatchFile("latch:///megs/SRR579291/SRR579291.contigs.fa"),
#             evaluation=LatchDir("latch:///megs/SRR579291/SRR579291_MetaQuast"),
#         ),
#         AssemblyOut(
#             sample_name="SRR579292",
#             assembly_data=LatchFile("latch:///megs/SRR579292/SRR579292.contigs.fa"),
#             evaluation=LatchDir("latch:///megs/SRR579292/SRR579292_MetaQuast"),
#         ),
#     ],
# )

# wf.assembly.assembly_wf(
# samples=[
#         Sample(
#             sample_name="SRR579291",
#             read1=LatchFile("s3://latch-public/test-data/4318/SRR579291_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/SRR579291_2.fastq"),
#         ),
#         Sample(
#             sample_name="SRR579292",
#             read1=LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
#             read2=LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
#         ),
#     ],
#     min_count= 2,
#     k_min= 21,
#     k_max= 141,
#     k_step= 12,
#     min_contig_len= 200,
# )
