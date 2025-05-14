import click
from ..version import __version__

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="hap.py")
def main():
    """Run hap.py: haplotype comparison tool."""
    click.echo(main.get_help(click.Context(main)))
