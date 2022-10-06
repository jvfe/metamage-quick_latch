import subprocess
from dataclasses import dataclass
from pathlib import Path
from random import sample
from typing import List

from dataclasses_json import dataclass_json
from latch import large_task, map_task, message, small_task, workflow
from latch.types import LatchDir, LatchFile

from .assembly import AssemblyOut
from .types import Sample


@dataclass
@dataclass_json
class BwAlignInput:
    assembly_index: LatchDir
    read_data: Sample


@dataclass
@dataclass_json
class JgiInput:
    assembly_bam: LatchFile
    sample_name: str


@dataclass
@dataclass_json
class MetaBatInput:
    sample_name: str
    assembly_data: LatchDir
    depth_file: LatchFile


@large_task
def bowtie_assembly_build(megahit_out: AssemblyOut) -> LatchDir:

    sample_name = megahit_out.sample_name
    assembly_name = f"{sample_name}.contigs.fa"
    assembly_fasta = Path(megahit_out.assembly_data.local_path, assembly_name)

    output_dir_name = f"{sample_name}_assembly_idx"
    output_dir = Path(output_dir_name).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    _bt_idx_cmd = [
        "bowtie2/bowtie2-build",
        str(assembly_fasta),
        f"{str(output_dir)}/{sample_name}",
        "--threads",
        "31",
    ]

    subprocess.run(_bt_idx_cmd)

    return LatchDir(str(output_dir), f"latch:///maggie/{sample_name}/{output_dir_name}")


@small_task
def organize_bwalign_inputs(
    assembly_idxs: List[LatchDir], samples: List[Sample]
) -> List[BwAlignInput]:

    inputs = []
    for sample, assembly_idx in zip(samples, assembly_idxs):
        cur_input = BwAlignInput(assembly_index=assembly_idx, read_data=sample)

        inputs.append(cur_input)

    return inputs


@large_task
def bowtie_assembly_align(bwalign_input: BwAlignInput) -> JgiInput:

    sample_name = bwalign_input.read_data.sample_name
    output_file_name = f"{sample_name}_assembly_sorted.bam"

    output_file = Path(output_file_name).resolve()

    _bt_cmd = [
        "bowtie2/bowtie2",
        "-x",
        f"{bwalign_input.assembly_idx.local_path}/{sample_name}",
        "-1",
        bwalign_input.read_data.read1.local_path,
        "-2",
        bwalign_input.read_data.read2.local_path,
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

    return JgiInput(
        assembly_bam=LatchFile(
            str(output_file), f"latch:///maggie/{sample_name}/{output_file_name}"
        ),
        sample_name=sample_name,
    )


@small_task
def summarize_contig_depths(jgi_input: JgiInput) -> LatchFile:

    sample_name = jgi_input.sample_name
    output_file_name = f"{sample_name}_depths.txt"
    output_file = Path(output_file_name).resolve()

    _jgi_cmd = [
        "jgi_summarize_bam_contig_depths",
        "--outputDepth",
        output_file_name,
        jgi_input.assembly_bam.local_path,
    ]

    subprocess.run(_jgi_cmd)

    return LatchFile(
        str(output_file), f"latch:///maggie/{sample_name}/{output_file_name}"
    )


@small_task
def organize_metabat_inputs(
    assembly_data: List[AssemblyOut],
    depth_files: List[LatchFile],
) -> List[MetaBatInput]:

    inputs = []
    for assembly, depth_file in zip(assembly_data, depth_files):
        cur_input = MetaBatInput(
            sample_name=assembly.sample_name,
            assembly_data=assembly.assembly_data,
            depth_file=depth_file,
        )
        inputs.append(cur_input)

    return inputs


@large_task
def metabat2(metabat_input: MetaBatInput) -> LatchDir:

    sample_name = metabat_input.sample_name
    assembly_name = f"{sample_name}.contigs.fa"
    assembly_fasta = Path(metabat_input.assembly_data.local_path, assembly_name)

    output_dir_name = f"METABAT/{sample_name}"
    output_dir = Path(output_dir_name).parent.resolve()

    _metabat_cmd = [
        "metabat2",
        "--saveCls",
        "-i",
        str(assembly_fasta),
        "-a",
        metabat_input.depth_file.local_path,
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

    return LatchDir(str(output_dir), f"latch:///maggie/{sample_name}/METABAT/")


@workflow
def binning_wf(samples: List[Sample], megahit_out: List[AssemblyOut]) -> List[LatchDir]:

    # Binning preparation
    built_assembly_idxs = map_task(bowtie_assembly_build)(megahit_out=megahit_out)

    bwalign_inputs = organize_bwalign_inputs(
        assembly_idxs=built_assembly_idxs, samples=samples
    )

    aligned_to_assembly = map_task(bowtie_assembly_align)(bwalign_input=bwalign_inputs)

    depth_files = map_task(summarize_contig_depths)(jgi_input=aligned_to_assembly)

    metabat_inputs = organize_metabat_inputs(
        assembly_data=megahit_out, depth_files=depth_files
    )

    # Binning
    binning_results = map_task(metabat2)(metabat_input=metabat_inputs)

    return binning_results
