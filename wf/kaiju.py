"""
Taxonomic classification of reads
"""

import subprocess
from pathlib import Path
from typing import Tuple

from latch import large_task, message, small_task, workflow
from latch.types import LatchDir, LatchFile

from .types import Sample, TaxonRank


@large_task
def taxonomy_classification_task(
    sample: Sample,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_db: LatchFile,
) -> LatchFile:
    """Classify metagenomic reads with Kaiju"""

    output_name = f"{sample.sample_name}_kaiju.out"
    kaiju_out = Path(output_name).resolve()

    _kaiju_cmd = [
        "kaiju",
        "-t",
        kaiju_ref_nodes.local_path,
        "-f",
        kaiju_ref_db.local_path,
        "-i",
        sample.read1.local_path,
        "-j",
        sample.read2.local_path,
        "-z",
        "2",
        "-o",
        str(kaiju_out),
    ]
    message(
        "info",
        {
            "title": "Taxonomically classifying reads with Kaiju",
            "body": f"Command: {' '.join(_kaiju_cmd)}",
        },
    )
    subprocess.run(_kaiju_cmd)

    return LatchFile(
        str(kaiju_out), f"latch:///metamage/{sample.sample_name}/kaiju/{output_name}"
    )


@small_task
def kaiju2table_task(
    kaiju_out: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    sample: Sample,
    taxon: TaxonRank,
) -> LatchFile:
    """Convert Kaiju output to TSV format"""

    output_name = f"{sample.sample_name}_kaiju.tsv"
    kaijutable_tsv = Path(output_name).resolve()

    _kaiju2table_cmd = [
        "kaiju2table",
        "-t",
        kaiju_ref_nodes.local_path,
        "-n",
        kaiju_ref_names.local_path,
        "-r",
        taxon.value,
        "-p",
        "-e",
        "-o",
        str(kaijutable_tsv),
        kaiju_out.local_path,
    ]

    subprocess.run(_kaiju2table_cmd)

    return LatchFile(
        str(kaijutable_tsv),
        f"latch:///metamage/{sample.sample_name}/kaiju/{output_name}",
    )


@small_task
def kaiju2krona_task(
    kaiju_out: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    sample: Sample,
) -> LatchFile:
    """Convert Kaiju output to Krona-readable format"""

    output_name = f"{sample.sample_name}_kaiju2krona.out"
    krona_txt = Path(output_name).resolve()

    _kaiju2krona_cmd = [
        "kaiju2krona",
        "-t",
        kaiju_ref_nodes.local_path,
        "-n",
        kaiju_ref_names.local_path,
        "-i",
        kaiju_out.local_path,
        "-o",
        str(krona_txt),
    ]

    subprocess.run(_kaiju2krona_cmd)

    return LatchFile(
        str(krona_txt), f"latch:///metamage/{sample.sample_name}/kaiju/{output_name}"
    )


@small_task
def plot_krona_task(krona_txt: LatchFile, sample: Sample) -> LatchFile:
    """Make Krona plot from Kaiju results"""
    output_name = f"{sample.sample_name}_krona.html"
    krona_html = Path(output_name).resolve()

    _kaiju2krona_cmd = ["ktImportText", "-o", str(krona_html), krona_txt.local_path]

    subprocess.run(_kaiju2krona_cmd)

    return LatchFile(
        str(krona_html), f"latch:///metamage/{sample.sample_name}/kaiju/{output_name}"
    )


@workflow
def kaiju_wf(
    sample: Sample,
    kaiju_ref_db: LatchFile,
    kaiju_ref_nodes: LatchFile,
    kaiju_ref_names: LatchFile,
    taxon_rank: TaxonRank,
) -> Tuple[LatchFile, LatchFile]:

    kaiju_out = taxonomy_classification_task(
        sample=sample,
        kaiju_ref_db=kaiju_ref_db,
        kaiju_ref_nodes=kaiju_ref_nodes,
    )
    kaiju2table_out = kaiju2table_task(
        kaiju_out=kaiju_out,
        sample=sample,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
        taxon=taxon_rank,
    )
    kaiju2krona_out = kaiju2krona_task(
        kaiju_out=kaiju_out,
        sample=sample,
        kaiju_ref_nodes=kaiju_ref_nodes,
        kaiju_ref_names=kaiju_ref_names,
    )
    krona_plot = plot_krona_task(krona_txt=kaiju2krona_out, sample=sample)

    return kaiju2table_out, krona_plot
