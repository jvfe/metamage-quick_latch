"""
Read assembly and evaluation for metagenomics data
"""

import subprocess
from pathlib import Path
from typing import Tuple

from latch import large_task, message, small_task, workflow
from latch.types import LatchDir

from .types import Sample


@large_task
def megahit(
    sample: Sample,
    min_count: int,
    k_min: int,
    k_max: int,
    k_step: int,
    min_contig_len: int,
) -> LatchDir:

    output_dir_name = "MEGAHIT"
    output_dir = Path(output_dir_name).resolve()

    _megahit_cmd = [
        "/root/megahit",
        "--min-count",
        str(min_count),
        "--k-min",
        str(k_min),
        "--k-max",
        str(k_max),
        "--k-step",
        str(k_step),
        "--out-dir",
        output_dir_name,
        "--out-prefix",
        sample.sample_name,
        "--min-contig-len",
        str(min_contig_len),
        "-1",
        sample.read1.local_path,
        "-2",
        sample.read2.local_path,
    ]
    message(
        "info",
        {
            "title": "Assembling reads into contigs with MEGAHit",
            "body": f"Command: {' '.join(_megahit_cmd)}",
        },
    )
    subprocess.run(_megahit_cmd)

    return LatchDir(
        str(output_dir), f"latch:///metamage/{sample.sample_name}/{output_dir_name}"
    )


@small_task
def metaquast(
    assembly_dir: LatchDir,
    sample: Sample,
) -> LatchDir:

    assembly_name = f"{sample.sample_name}.contigs.fa"
    assembly_fasta = Path(assembly_dir.local_path, assembly_name)

    output_dir_name = "MetaQuast"
    output_dir = Path(output_dir_name).resolve()

    _metaquast_cmd = [
        "/root/metaquast.py",
        "--rna-finding",
        "--no-sv",
        "--max-ref-number",
        "0",
        "-l",
        sample.sample_name,
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

    return LatchDir(
        str(output_dir), f"latch:///metamage/{sample.sample_name}/{output_dir_name}"
    )


@workflow
def assembly_wf(
    sample: Sample,
    min_count: int,
    k_min: int,
    k_max: int,
    k_step: int,
    min_contig_len: int,
) -> Tuple[LatchDir, LatchDir]:

    # Assembly
    assembly_dir = megahit(
        sample=sample,
        min_count=min_count,
        k_min=k_min,
        k_max=k_max,
        k_step=k_step,
        min_contig_len=min_contig_len,
    )
    metassembly_results = metaquast(assembly_dir=assembly_dir, sample=sample)

    return assembly_dir, metassembly_results
