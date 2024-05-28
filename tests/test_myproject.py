from pathlib import Path

import pytest

from click.testing import CliRunner
from cogent3.core.alignment import SequenceCollection

from myproject.cli import main, unique_kmers, unique_kmers_for_cli


DATADIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def runner():
    return CliRunner()


def test_unique_kmers_for_cli(runner):
    seqs_input = '{"s1": "AT", "s2": "AC", "s3": "AA"}'
    kmer_size = 2
    args = ["unique_kmers_for_cli", "-s", seqs_input, "-k", str(kmer_size)]

    result = runner.invoke(main, args, catch_exceptions=False)
    assert result.exit_code == 0, result.output

    expected_output = (
        "Unique k-mers of s1: AT\n"
        "Unique k-mers of s2: AC\n"
        "Unique k-mers of s3: AA\n"
    )

    for line in expected_output.strip().split("\n"):
        assert line in result.output


def test_unique_kmers():
    seqs = SequenceCollection(
        {"s1": "ATAATCC", "s2": "ATGATCC", "s3": "ATACTCC"}, moltype="dna"
    )
    k = 3
    result = unique_kmers(seqs, k)

    for seq_name, kmers in result.items():
        for kmer in kmers:
            assert len(kmer) == k
