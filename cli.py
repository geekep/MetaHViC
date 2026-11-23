import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path

import rich_click as click
from rich.console import Console
from rich.padding import Padding
from rich.panel import Panel

import MetaHViC
from MetaHViC.utils import get_n_available_cpus

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.MAX_WIDTH = None
click.rich_click.STYLE_OPTIONS_TABLE_BOX = "SIMPLE"
click.rich_click.STYLE_COMMANDS_TABLE_SHOW_LINES = True
click.rich_click.STYLE_COMMANDS_TABLE_PAD_EDGE = True
click.rich_click.STYLE_COMMANDS_TABLE_BOX = "SIMPLE"

click.rich_click.COMMAND_GROUPS = {
    "MetaHViC": [
        {
            "name": "Modules",
            "commands": [
                # "annotate",
                # "find-proviruses",
                # "find-proviruses-plm",
                # "marker-classification",
                # "nn-classification",
                "host-classification",
                # "aggregated-classification",
                # "score-calibration",
                # "summary",
            ],
        },
    ]
}

click.rich_click.OPTION_GROUPS = {
    "MetaHViC host-classification": [
        {
            "name": "Basic options",
            "options": ["--cleanup", "--restart", "--threads", "--verbose"],
        },
        {
            "name": "Advanced options",
            "options": ["--single-window", "--batch-size"],
        },
        {
            "name": "Other",
            "options": ["--help", "--version"],
        },
    ],
}


def use_preset(ctx, param, value):
    if value is not None and any(
        ctx.get_parameter_source(i).name == "COMMANDLINE"
        for i in [
            "min_score",
            "max_fdr",
            "min_number_genes",
            "min_plasmid_marker_enrichment",
            "min_virus_marker_enrichment",
            "min_plasmid_hallmarks",
            "min_plasmid_hallmarks_short_seqs",
            "min_virus_hallmarks",
            "min_virus_hallmarks_short_seqs",
            "max_uscg",
        ]
    ):
        console = Console()
        console.print(
            Padding(
                "Try [blue]'MetaHViC end-to-end -h'[/] for help.",
                (0, 1, 0, 1),
            ),
            style="dim",
        )
        console.print(
            Padding(
                Panel(
                    "You cannot use filtering options ([bold cyan]--min-score[/], "
                    "[bold cyan]--max-fdr[/], etc.) together with a preset "
                    "([bold cyan]--conservative[/] or [bold cyan]--relaxed[/]). "
                    "Please modify your command to remove either the preset or the "
                    "filtering option(s).",
                    border_style="red",
                    title="Error",
                    title_align="left",
                ),
                (0, 0, 0, 0),
            )
        )
        sys.exit(1)
    if value is False:
        set_preset_values(ctx, 0, 1, 0, -100, -100, 0, 0, 0, 0, 100)
    elif value is True:
        set_preset_values(ctx, 0.8, 0.05, 1, 1.5, 1.5, 1, 1, 1, 1, 2)


def set_preset_values(
    ctx,
    min_score,
    max_fdr,
    min_number_genes,
    min_plasmid_marker_enrichment,
    min_virus_marker_enrichment,
    min_plasmid_hallmarks,
    min_plasmid_hallmarks_short_seqs,
    min_virus_hallmarks,
    min_virus_hallmarks_short_seqs,
    max_uscg,
):
    ctx.params["min_score"] = min_score
    ctx.params["max_fdr"] = max_fdr
    ctx.params["min_number_genes"] = min_number_genes
    ctx.params["min_plasmid_marker_enrichment"] = min_plasmid_marker_enrichment
    ctx.params["min_virus_marker_enrichment"] = min_virus_marker_enrichment
    ctx.params["min_plasmid_hallmarks"] = min_plasmid_hallmarks
    ctx.params["min_plasmid_hallmarks_short_seqs"] = min_plasmid_hallmarks_short_seqs
    ctx.params["min_virus_hallmarks"] = min_virus_hallmarks
    ctx.params["min_virus_hallmarks_short_seqs"] = min_virus_hallmarks_short_seqs
    ctx.params["max_uscg"] = max_uscg


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(prog_name="MetaHViC")
def cli():
    """
    [cyan]MetaHViC[/cyan]: a Host-informed Virus Classifier for Metagenomic contigs
    Read the documentation at:
    [link=https://github.com/geekep/MetaHViC]https://github.com/geekep/MetaHViC[/link]
    """


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("input", type=click.Path(path_type=Path, exists=True))
@click.argument("output", type=click.Path(path_type=Path))
@click.version_option(prog_name="MetaHViC")
@click.option(
    "--restart",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite existing intermediate files.",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=get_n_available_cpus(),
    show_default=True,
    help="Number of threads to use.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    is_flag=True,
    default=True,
    show_default=True,
    help="Display the execution log.",
)
@click.option(
    "--cleanup",
    is_flag=True,
    default=False,
    show_default=True,
    help="Delete intermediate files after execution.",
)
@click.option(
    "--single-window",
    is_flag=True,
    default=False,
    show_default=True,
    help="""Use only the first window (6,000 bases) of each sequence to perform the classification.
            This will make execution faster and reduce memory usage, but the classification accuracy will decrease.""",
)
@click.option(
    "--batch-size",
    type=int,
    default=256,
    show_default=True,
    help="""Number of data points per batch of prediction. Use a smaller value to reduce memory consumption at the cost of speed.""",
)
def host_classification(input, output, single_window, batch_size, restart, threads, verbose, cleanup):
    """
    Classify the sequences in the [u]INPUT[/u] file (FASTA format) using the IGLOO neural network of MetaHViC and write the results to the [u]OUTPUT[/u] directory.
    """
    MetaHViC.host_classification.main(input, output, single_window, batch_size, restart, threads, verbose, cleanup)


if __name__ == "__main__":
    cli()
