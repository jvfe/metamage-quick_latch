import subprocess
from pathlib import Path

from latch import large_task, message, small_task, workflow
from latch.types import LatchDir, LatchFile

from .types import Sample


@large_task
def bowtie_assembly_build(assembly_dir: LatchDir, sample: Sample) -> LatchDir:

    assembly_name = f"{sample.sample_name}.contigs.fa"
    assembly_fasta = Path(assembly_dir.local_path, assembly_name)

    output_dir_name = f"{sample.sample_name}_assembly_idx"
    output_dir = Path(output_dir_name).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    _bt_idx_cmd = [
        "bowtie2/bowtie2-build",
        str(assembly_fasta),
        f"{str(output_dir)}/{sample.sample_name}",
        "--threads",
        "31",
    ]

    subprocess.run(_bt_idx_cmd)

    return LatchDir(
        str(output_dir), f"latch:///metamage/{sample.sample_name}/{output_dir_name}"
    )


@large_task
def bowtie_assembly_align(
    assembly_idx: LatchDir,
    sample: Sample,
) -> LatchFile:

    output_file_name = f"{sample.sample_name}_assembly_sorted.bam"

    output_file = Path(output_file_name).resolve()

    _bt_cmd = [
        "bowtie2/bowtie2",
        "-x",
        f"{assembly_idx.local_path}/{sample.sample_name}",
        "-1",
        sample.read1.local_path,
        "-2",
        sample.read2.local_path,
        "--threads",
        "31",
    ]

    bt_align_out = subprocess.Popen(
        _bt_cmd,
        stdout=subprocess.PIPE,
    )

    _sam_convert_cmd = [
        "samtools",
        "view",
        "-@",
        "31",
        "-bS",
    ]

    sam_convert_out = subprocess.Popen(
        _sam_convert_cmd, stdin=bt_align_out.stdout, stdout=subprocess.PIPE
    )

    _sam_sort_cmd = [
        "samtools",
        "sort",
        "-@",
        "31",
        "-o",
        output_file_name,
    ]

    subprocess.run(
        _sam_sort_cmd,
        stdin=sam_convert_out.stdout,
    )

    return LatchFile(
        str(output_file), f"latch:///metamage/{sample.sample_name}/{output_file_name}"
    )


@small_task
def summarize_contig_depths(assembly_bam: LatchFile, sample: Sample) -> LatchFile:

    output_file_name = f"{sample.sample_name}_depths.txt"
    output_file = Path(output_file_name).resolve()

    _jgi_cmd = [
        "jgi_summarize_bam_contig_depths",
        "--outputDepth",
        output_file_name,
        assembly_bam.local_path,
    ]

    subprocess.run(_jgi_cmd)

    return LatchFile(
        str(output_file), f"latch:///metamage/{sample.sample_name}/{output_file_name}"
    )


@large_task
def metabat2(
    assembly_dir: LatchDir,
    depth_file: LatchFile,
    sample: Sample,
) -> LatchDir:

    assembly_name = f"{sample.sample_name}.contigs.fa"
    assembly_fasta = Path(assembly_dir.local_path, assembly_name)

    output_dir_name = f"METABAT/{sample.sample_name}"
    output_dir = Path(output_dir_name).parent.resolve()

    _metabat_cmd = [
        "metabat2",
        "--saveCls",
        "-i",
        str(assembly_fasta),
        "-a",
        depth_file.local_path,
        "-o",
        output_dir_name,
    ]
    message(
        "info",
        {
            "title": "Binning contigs with MetaBat2",
            "body": f"Command: {' '.join(_metabat_cmd)}",
        },
    )
    subprocess.run(_metabat_cmd)

    return LatchDir(str(output_dir), f"latch:///metamage/{sample.sample_name}/METABAT/")


@workflow
def binning_wf(sample: Sample, assembly_dir: LatchDir) -> LatchDir:

    # Binning preparation
    built_assembly_idx = bowtie_assembly_build(assembly_dir=assembly_dir, sample=sample)
    aligned_to_assembly = bowtie_assembly_align(
        assembly_idx=built_assembly_idx, sample=sample
    )
    depth_file = summarize_contig_depths(
        assembly_bam=aligned_to_assembly, sample=sample
    )

    # Binning
    binning_results = metabat2(
        assembly_dir=assembly_dir, depth_file=depth_file, sample=sample
    )

    return binning_results
