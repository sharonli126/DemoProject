from pathlib import Path

import click

# added
from cogent3.core.alignment import SequenceCollection
from scitrack import CachingLogger


__author__ = "Tsz Ching Li"
__copyright__ = ["Tsz Ching Li"]
__license__ = "BSD"
__version__ = "2024.5.21"  # A DATE BASED VERSION
__maintainer__ = "YOUR NAME"
__email__ = "u7630977@anu.edu.au"
__status__ = "alpha"


LOGGER = CachingLogger()


@click.group()
@click.version_option(__version__)  # add version option
def main():
    """This is a bioinformatic CLI tool to find unique kmers among sequences."""
    pass


def unique_kmers(seqs: SequenceCollection, k: int = 2):
    # Returns a dictionary {seqname: set(str, ...)} of unique k-mers for each sequence.
    result = {}
    kmers = dict()
    seqs_dict = seqs.to_dict()
    for name, seq in seqs_dict.items():
        kmer_set = set()
        for i in range(len(seq) - k + 1):
            kmer_set.add(seq[i : i + k])
        kmers[name] = kmer_set

    for name, kmer in kmers.items():
        other_kmers = set()
        for n in kmers:
            if n != name:
                other_kmers.update(kmers[n])
        result[name] = kmer - other_kmers

    return result


# custom parsers / validators
def parse_sequences(ctx, param, value):
    seqs = eval(value)
    return SequenceCollection(seqs, moltype="dna")


# the no_args_is_help=True means help is displayed if a
# user doesn't provide any arguments to a subcommand.
# Should be a click default I think!
@main.command(
    name="unique_kmers_for_cli", no_args_is_help=True
)  # Ensure command is registered correctly
@click.option(
    "-s",
    "--seqs",
    required=True,
    callback=parse_sequences,
    help='Input sequences in the format {"s1": "ATAATCC", "s2": "ATGATCC", "s3": "ATACTCC"}',
)
@click.option(
    "-k",
    "--kmer_size",
    required=True,
    type=int,
    help="Size of the k-mers",
)
def unique_kmers_for_cli(seqs, kmer_size):
    """CLI command to find unique k-mers."""
    unique_kmers_result = unique_kmers(seqs, kmer_size)

    for seq_name, kmers in unique_kmers_result.items():
        click.echo(f"Unique k-mers of {seq_name}: {', '.join(kmers)}")


if __name__ == "__main__":
    main()
