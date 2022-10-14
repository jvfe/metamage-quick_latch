"""
Read assembly and evaluation for metagenomics data
"""

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

from dataclasses_json import dataclass_json
from latch import large_task, map_task, message, small_task, workflow
from latch.types import LatchDir, LatchFile

from .types import Sample


@dataclass_json
@dataclass
class MegaHitInput:
    sample_name: str
    read1: LatchFile
    read2: LatchFile
    min_count: int
    k_min: int
    k_max: int
    k_step: int
    min_contig_len: int


@dataclass_json
@dataclass
class MegaHitOut:
    sample_name: str
    assembly_data: LatchDir


@dataclass_json
@dataclass
class AssemblyOut(MegaHitOut):
    evaluation: LatchDir


@small_task
def organize_megahit_inputs(
    samples: List[Sample],
    min_count: int,
    k_min: int,
    k_max: int,
    k_step: int,
    min_contig_len: int,
) -> List[MegaHitInput]:

    inputs = []
    for sample in samples:
        cur_input = MegaHitInput(
            sample_name=sample.sample_name,
            read1=sample.read1,
            read2=sample.read2,
            min_count=min_count,
            k_min=k_min,
            k_max=k_max,
            k_step=k_step,
            min_contig_len=min_contig_len,
        )

        inputs.append(cur_input)

    return inputs


@large_task
def megahit(megahit_input: MegaHitInput) -> LatchDir:

    sample_name = megahit_input.sample_name
    output_dir_name = "MEGAHIT"
    output_dir = Path(output_dir_name).resolve()

    _megahit_cmd = [
        "/root/megahit",
        "--min-count",
        str(megahit_input.min_count),
        "--k-min",
        str(megahit_input.k_min),
        "--k-max",
        str(megahit_input.k_max),
        "--k-step",
        str(megahit_input.k_step),
        "--out-dir",
        output_dir_name,
        "--out-prefix",
        sample_name,
        "--min-contig-len",
        str(megahit_input.min_contig_len),
        "-1",
        megahit_input.read1.local_path,
        "-2",
        megahit_input.read2.local_path,
    ]
    message(
        "info",
        {
            "title": "Assembling reads into contigs with MEGAHit",
            "body": f"Command: {' '.join(_megahit_cmd)}",
        },
    )
    subprocess.run(_megahit_cmd)

    return LatchDir(str(output_dir), f"latch:///megs/{sample_name}/{output_dir_name}")


@small_task
def organize_megahit_outs(
    samples: List[Sample], assembly_data: List[LatchDir]
) -> List[MegaHitOut]:

    outs = []
    for sample, assembly in zip(samples, assembly_data):
        cur_out = MegaHitOut(sample_name=sample.sample_name, assembly_data=assembly)
        outs.append(cur_out)

    return outs


@small_task
def metaquast(megahit_out: MegaHitOut) -> LatchDir:

    sample_name = megahit_out.sample_name
    assembly_name = f"{sample_name}.contigs.fa"
    assembly_fasta = Path(megahit_out.assembly_data.local_path, assembly_name)

    output_dir_name = "MetaQuast"
    output_dir = Path(output_dir_name).resolve()

    _metaquast_cmd = [
        "/root/metaquast.py",
        "--rna-finding",
        "--no-sv",
        "--max-ref-number",
        "0",
        "-l",
        sample_name,
        "-o",
        output_dir_name,
        str(assembly_fasta),
    ]
    message(
        "info",
        {
            "title": "Evaluating assembly with MetaQuast",
            "body": f"Command: {' '.join(_metaquast_cmd)}",
        },
    )
    subprocess.run(_metaquast_cmd)

    return LatchDir(str(output_dir), f"latch:///megs/{sample_name}/{output_dir_name}")


@small_task
def organize_assembly_outs(
    megahit_outs: List[MegaHitOut], metaquast_results: List[LatchDir]
) -> List[AssemblyOut]:

    outs = []

    for assembly, evaluation in zip(megahit_outs, metaquast_results):

        cur_out = AssemblyOut(
            sample_name=assembly.sample_name,
            assembly_data=assembly.assembly_data,
            evaluation=evaluation,
        )
        outs.append(cur_out)

    return outs


@workflow
def assembly_wf(
    samples: List[Sample],
    min_count: int,
    k_min: int,
    k_max: int,
    k_step: int,
    min_contig_len: int,
) -> List[AssemblyOut]:

    megahit_inputs = organize_megahit_inputs(
        samples=samples,
        min_count=min_count,
        k_min=k_min,
        k_max=k_max,
        k_step=k_step,
        min_contig_len=min_contig_len,
    )

    # Assembly
    assembly_data = map_task(megahit)(megahit_input=megahit_inputs)

    megahit_outs = organize_megahit_outs(samples=samples, assembly_data=assembly_data)

    metaquast_results = map_task(metaquast)(megahit_out=megahit_outs)

    return organize_assembly_outs(
        megahit_outs=megahit_outs, metaquast_results=metaquast_results
    )
