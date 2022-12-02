import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List

from dataclasses_json import dataclass_json
from latch import map_task, medium_task, message, small_task, workflow
from latch.types import LatchDir, LatchFile

from .assembly import AssemblyOut
from .types import ProdigalOutput, fARGeneModel


@dataclass_json
@dataclass
class FunctionalInput:
    sample_name: str
    assembly_data: LatchFile
    prodigal_output_format: ProdigalOutput
    fargene_hmm_model: fARGeneModel


@dataclass_json
@dataclass
class FunctionalOutput:
    sample_name: str
    prodigal_result: LatchDir
    macrel_result: LatchDir
    fargene_result: LatchDir
    gecco_result: LatchDir


@small_task
def organize_functional_inputs(
    assembly_data: List[AssemblyOut],
    prodigal_output_format: ProdigalOutput,
    fargene_hmm_model: fARGeneModel,
) -> List[FunctionalInput]:

    ins = []
    for assembly in assembly_data:
        ins.append(
            FunctionalInput(
                sample_name=assembly.sample_name,
                assembly_data=assembly.assembly_data,
                prodigal_output_format=prodigal_output_format,
                fargene_hmm_model=fargene_hmm_model,
            )
        )

    return ins


@small_task
def macrel(functional_in: FunctionalInput) -> LatchDir:

    # Assembly data
    sample_name = functional_in.sample_name
    assembly_fasta = Path(functional_in.assembly_data.local_path)

    output_dir_name = "macrel_results"
    outdir = Path(output_dir_name).resolve()

    _macrel_cmd = [
        "macrel",
        "contigs",
        "--fasta",
        str(assembly_fasta),
        "--output",
        str(outdir),
        "--tag",
        sample_name,
        "--log-file",
        f"{str(outdir)}/{sample_name}_log.txt",
        "--threads",
        "8",
    ]

    subprocess.run(_macrel_cmd)

    return LatchDir(str(outdir), f"latch:///metamage/{sample_name}/{output_dir_name}")


@small_task
def fargene(functional_in: FunctionalInput) -> LatchDir:

    # Assembly data
    sample_name = functional_in.sample_name
    assembly_fasta = Path(functional_in.assembly_data.local_path)

    output_dir_name = "fargene_results"
    outdir = Path(output_dir_name).resolve()

    _fargene_cmd = [
        "fargene",
        "-i",
        str(assembly_fasta),
        "--hmm-model",
        functional_in.fargene_hmm_model.value,
        "-o",
        output_dir_name,
        "-p",
        "8",
    ]

    subprocess.run(_fargene_cmd)

    return LatchDir(str(outdir), f"latch:///metamage/{sample_name}/{output_dir_name}")


@small_task
def gecco(functional_in: FunctionalInput) -> LatchDir:

    # Assembly data
    sample_name = functional_in.sample_name
    assembly_fasta = Path(functional_in.assembly_data.local_path)

    output_dir_name = "gecco_results"
    outdir = Path(output_dir_name).resolve()

    _gecco_cmd = [
        "gecco",
        "run",
        "-g",
        str(assembly_fasta),
        "-o",
        output_dir_name,
        "-j",
        "4",
        "--force-tsv",
    ]

    subprocess.run(_gecco_cmd)

    return LatchDir(str(outdir), f"latch:///metamage/{sample_name}/{output_dir_name}")


@medium_task
def prodigal(functional_in: FunctionalInput) -> LatchDir:

    # Assembly data
    sample_name = functional_in.sample_name
    assembly_fasta = Path(functional_in.assembly_data.local_path)

    # A reference to our output.
    output_format = functional_in.prodigal_output_format
    output_dir_name = "prodigal_results"
    output_dir = Path(output_dir_name).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir.joinpath(f"{sample_name}.{output_format.value}")
    output_proteins = output_dir.joinpath(f"{sample_name}.faa")
    output_genes = output_dir.joinpath(f"{sample_name}.fna")
    output_scores = output_dir.joinpath(f"{sample_name}.cds")

    _prodigal_cmd = [
        "/root/prodigal",
        "-i",
        str(assembly_fasta),
        "-f",
        output_format.value,
        "-o",
        str(output_file),
        "-a",
        str(output_proteins),
        "-d",
        str(output_genes),
        "-s",
        str(output_scores),
    ]

    subprocess.run(_prodigal_cmd)

    return LatchDir(
        str(output_dir), f"latch:///metamage/{sample_name}/{output_dir_name}"
    )


@small_task
def organize_functional_outputs(
    inputs: List[FunctionalInput],
    prodigal_results: List[LatchDir],
    macrel_results: List[LatchDir],
    fargene_results: List[LatchDir],
    gecco_results: List[LatchDir],
) -> List[FunctionalOutput]:

    outs = []
    for sample, prod, macr, farg, gecc in zip(
        inputs, prodigal_results, macrel_results, fargene_results, gecco_results
    ):

        cur_out = FunctionalOutput(
            sample_name=sample.sample_name,
            prodigal_result=prod,
            macrel_result=macr,
            fargene_result=farg,
            gecco_result=gecc,
        )
        outs.append(cur_out)

    return outs


@workflow
def functional_wf(
    assembly_data: List[AssemblyOut],
    prodigal_output_format: ProdigalOutput,
    fargene_hmm_model: fARGeneModel,
) -> List[FunctionalOutput]:

    functional_ins = organize_functional_inputs(
        assembly_data=assembly_data,
        prodigal_output_format=prodigal_output_format,
        fargene_hmm_model=fargene_hmm_model,
    )

    # Functional annotation
    prodigal_results = map_task(prodigal)(functional_in=functional_ins)
    macrel_results = map_task(macrel)(functional_in=functional_ins)
    fargene_results = map_task(fargene)(functional_in=functional_ins)
    gecco_results = map_task(gecco)(functional_in=functional_ins)

    func_outs = organize_functional_outputs(
        inputs=functional_ins,
        prodigal_results=prodigal_results,
        macrel_results=macrel_results,
        fargene_results=fargene_results,
        gecco_results=gecco_results,
    )

    return func_outs
