from pathlib import Path

from nucth import data_root


def test_data_root_contains_required_files():
    root = data_root()
    assert isinstance(root, Path)
    assert (root / "GPN.txt").is_file()
    assert (root / "Masses2016.txt").is_file()
    assert (root / "parameters").is_dir()


def test_parameter_sets_are_all_present():
    expected = {
        "fsu020", "fsugold2", "iufsu", "parameter_garnet", "parameter_gold",
        "parameter_linear", "parameter_nl3", "rmf022", "rmf028", "rmf032",
        "tamufsua", "tamufsub", "tamufsuc",
    }
    found = {p.stem for p in (data_root() / "parameters").glob("*.txt")}
    assert found == expected
