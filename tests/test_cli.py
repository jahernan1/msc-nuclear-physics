import numpy as np
import pytest

from nucth.cli import build_parser, run_convert


def test_parser_exposes_all_five_subcommands():
    parser = build_parser()
    for command in ("convert", "structure", "formfactors", "cevns", "all"):
        args = parser.parse_args([command, "--nucleus", "Ca40", "--model", "fsugold",
                                  "--neutrons", "20", "--protons", "20"])
        assert args.command == command


def test_parser_requires_a_nucleus():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["structure"])


def test_run_convert_writes_a_five_column_potential(tmp_path):
    raw = tmp_path / "raw_Test_model.dat"
    raw.write_text(
        "0 424.906 346.759 -3.82645 10.5182\n"
        "0.025 424.592 346.265 -3.82123 10.5176\n"
        "0.05 424.100 345.900 -3.80000 10.5100\n"
    )
    out = tmp_path / "Test_model_potential.dat"

    run_convert(raw, out)

    assert out.is_file()
    assert np.loadtxt(out).shape == (3, 5)
