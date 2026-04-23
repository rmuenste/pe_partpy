import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
PE_PARTPY_ROOT = REPO_ROOT / "pe_partpy"
if str(PE_PARTPY_ROOT) not in sys.path:
    sys.path.insert(0, str(PE_PARTPY_ROOT))

from partitioner.part import AxisBasedPartitioning, axis_cuts_partitioning, get_grid
from partitioner.part_main import AxisCutsSpec, checkParameters, parse_strategy


UNIT_CUBE_PROJECT = REPO_ROOT / "workflow" / "generators" / "unit_cube_27_case" / "file.prj"
UNIT_CUBE_GRID = REPO_ROOT / "workflow" / "generators" / "unit_cube_27_case" / "unit_cube_27.tri"


def _build_args(strategy: str, subdivision_spec: str, n_part: int, project_file: Path = UNIT_CUBE_PROJECT):
    return SimpleNamespace(
        n_part=n_part,
        strategy=strategy,
        subdivision_spec=subdivision_spec,
        mesh_name="TESTCASE",
        project_file=str(project_file),
        format="v2",
    )


def _single_axis_grid(x_planes):
    coords = []
    for x in x_planes:
        coords.extend(
            [
                (x, 0.0, 0.0),
                (x, 1.0, 0.0),
                (x, 1.0, 1.0),
                (x, 0.0, 1.0),
            ]
        )

    elements = []
    for plane_idx in range(len(x_planes) - 1):
        start = plane_idx * 4 + 1
        next_start = (plane_idx + 1) * 4 + 1
        elements.append(
            (
                start,
                start + 1,
                start + 2,
                start + 3,
                next_start,
                next_start + 1,
                next_start + 2,
                next_start + 3,
            )
        )

    return (
        len(elements),
        len(coords),
        tuple(coords),
        tuple(elements),
        tuple(0 for _ in coords),
    )


def test_axis_uniform_accepts_count_syntax():
    n_part, strategy, n_sub_part, mesh_name, project_file = checkParameters(
        _build_args("axis_uniform", "x3-y3-z3", 27)
    )
    assert n_part == 27
    assert strategy.name == "axis_uniform"
    assert n_sub_part == (3, 3, 3)


def test_axis_cuts_accepts_normalized_syntax():
    _, strategy, n_sub_part, _, _ = checkParameters(
        _build_args("axis_cuts", "x[0.2,0.5]-y[]-z[]", 3)
    )
    assert strategy.name == "axis_cuts"
    assert isinstance(n_sub_part, AxisCutsSpec)
    assert n_sub_part.cuts == ((0.2, 0.5), (), ())
    assert n_sub_part.absolute_axes == (False, False, False)


def test_axis_cuts_accepts_absolute_syntax():
    _, _, n_sub_part, _, _ = checkParameters(
        _build_args("axis_cuts", "x@[0.2,0.5]-y[]-z[]", 3)
    )
    assert isinstance(n_sub_part, AxisCutsSpec)
    assert n_sub_part.absolute_axes == (True, False, False)


def test_wrong_axis_order_is_rejected():
    with pytest.raises(SystemExit, match="axis order"):
        checkParameters(_build_args("axis_uniform", "y3-x3-z3", 27))


def test_axis_uniform_rejects_bracket_syntax():
    with pytest.raises(SystemExit, match="axis_uniform"):
        checkParameters(_build_args("axis_uniform", "x[0.2]-y[]-z[]", 2))


def test_axis_cuts_rejects_count_syntax():
    with pytest.raises(SystemExit, match="axis_cuts"):
        checkParameters(_build_args("axis_cuts", "x3-y1-z1", 3))


def test_axis_cuts_rejects_non_monotone_values():
    with pytest.raises(SystemExit, match="strictly increasing"):
        checkParameters(_build_args("axis_cuts", "x[0.5,0.2]-y[]-z[]", 3))


def test_axis_cuts_rejects_out_of_range_normalized_values():
    with pytest.raises(SystemExit, match="strictly inside"):
        checkParameters(_build_args("axis_cuts", "x[0.0,0.5]-y[]-z[]", 3))


def test_partition_count_mismatch_is_rejected():
    with pytest.raises(SystemExit, match="does not match"):
        checkParameters(_build_args("axis_uniform", "x3-y3-z3", 8))


def test_strategy_mapping_for_metis_and_planes():
    assert parse_strategy("metis_recursive").primary_metis_method == 1
    assert parse_strategy("metis_recursive_reversed").reverse_order is True
    assert parse_strategy("plane_dual").primary_strategy == "plane_dual"


def test_axis_cuts_matches_uniform_partitioning_when_cuts_align():
    grid = get_grid(UNIT_CUBE_GRID)
    uniform = AxisBasedPartitioning(grid, (3, 3, 3), -4)
    cuts = AxisCutsSpec(
        cuts=((1.0 / 3.0, 2.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0), (1.0 / 3.0, 2.0 / 3.0)),
        absolute_axes=(False, False, False),
    )
    explicit = axis_cuts_partitioning(grid, cuts)
    assert explicit == uniform


def test_axis_cuts_geometry_allocates_all_three_regions():
    grid = _single_axis_grid((0.0, 0.2, 0.5, 1.0))
    cuts = AxisCutsSpec(cuts=((0.2, 0.5), (), ()), absolute_axes=(True, False, False))
    partition = axis_cuts_partitioning(grid, cuts)
    assert partition == ([1, 1, 1], [1, 1, 2], [1, 1, 3])


def test_help_lists_named_strategies_and_legacy_codes():
    result = subprocess.run(
        [sys.executable, str(PE_PARTPY_ROOT / "PyPartitioner.py"), "--help"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "axis_cuts" in result.stdout
    assert "metis_recursive_reversed" in result.stdout
    assert "1 -> metis_recursive" in result.stdout
    assert "-7 -> plane_ring" in result.stdout
