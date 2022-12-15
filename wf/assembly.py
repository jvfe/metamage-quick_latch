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
    read_data: Sample
    min_count: int
    k_min: int
    k_max: int
    k_step: int
    min_contig_len: int


@dataclass_json
@dataclass
class MegaHitOut:
    sample_name: str
    assembly_data: LatchFile


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
            read_data=sample,
            min_count=min_count,
            k_min=k_min,
            k_max=k_max,
            k_step=k_step,
            min_contig_len=min_contig_len,
        )

        inputs.append(cur_input)

    return inputs


@large_task
def megahit(megahit_input: MegaHitInput) -> MegaHitOut:

    sample_name = megahit_input.read_data.sample_name
    output_dir_name = f"{sample_name}_MEGAHIT"

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
        megahit_input.read_data.read1.local_path,
        "-2",
        megahit_input.read_data.read2.local_path,
    ]

    subprocess.run(_megahit_cmd)

    megahit_output = Path(output_dir_name, f"{sample_name}.contigs.fa").resolve()

    return MegaHitOut(
        sample_name=sample_name,
        assembly_data=LatchFile(
            str(megahit_output),
            f"latch:///metamage/{sample_name}/MEGAHIT/{sample_name}.contigs.fa",
        ),
    )


@small_task
def metaquast(megahit_out: MegaHitOut) -> LatchDir:

    print(megahit_out.assembly_data.local_path)
    sample_name = megahit_out.sample_name
    assembly_fasta = megahit_out.assembly_data.local_path

    output_dir_name = f"{sample_name}_MetaQuast"
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

    subprocess.run(_metaquast_cmd)

    return LatchDir(str(output_dir), f"latch:///metamage/{sample_name}/{output_dir_name}")


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

    metaquast_results = map_task(metaquast)(megahit_out=assembly_data)

    return organize_assembly_outs(
        megahit_outs=assembly_data, metaquast_results=metaquast_results
    )
