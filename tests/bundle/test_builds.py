"""Guards of the build variants, read from the feature bundle (bars: tools/bundle/bars.yaml).

Noise floor: cold and cold2 (same binary, same input) are byte-identical except logstats.dat.
S6  lean build (EXTRA_TRACERS=0): the flow and oxygen outputs are byte-identical to the default build.
S7  MPI (4 ranks, lean): velocity within 5 % of the serial lean run.
"""
from pathlib import Path

from tests.bundle import _series as s
from tests.bundle.conftest import case


def test_noise_floor_two_identical_runs_are_bitwise_equal(bundle, bars):
    bar = bars["cold2"]
    a, b = case(bundle, "cold")["files"], case(bundle, "cold2")["files"]
    a = {k: v for k, v in a.items() if k not in bar["exclude"] and k != "params.json"}
    b = {k: v for k, v in b.items() if k not in bar["exclude"] and k != "params.json"}
    assert set(a) == set(b), f"different output files: {sorted(set(a) ^ set(b))}"
    diff = sorted(k for k in a if a[k] != b[k])
    assert not diff, f"cold vs cold2 differ in {diff}"


def test_s6_lean_build_leaves_flow_and_oxygen_unchanged(bundle, bars):
    bar = bars["lean"]
    cold, lean = case(bundle, "cold"), case(bundle, "lean")
    diff = [f for f in bar["files"] if cold["files"].get(f) != lean["files"].get(f)]
    assert not diff, f"lean differs from cold in {diff}"
    for f, cols in bar["columns_of"].items():
        ha, da = s.table(cold["dir"], f)
        hb, db = s.table(lean["dir"], f)
        a = da[:, [ha.index(c) for c in cols]]
        b = db[:, [hb.index(c) for c in cols]]
        assert a.shape == b.shape and (a == b).all(), f"{f} columns {cols} differ between cold and lean"
    print(f"\nwall time: cold {cold['wall_s']:.1f} s, lean {lean['wall_s']:.1f} s "
          f"({1 - lean['wall_s'] / cold['wall_s']:.0%} saved)")


def test_s7_mpi_matches_serial(bundle, bars):
    lean, mpi = case(bundle, "lean"), case(bundle, "mpi4")
    r = s.compare(bars["mpi4"], Path(lean["dir"]), Path(mpi["dir"]))
    print("\nmpi4:", {m: f"{d:.2e}" for m, (_, _, d, _) in r.items()},
          f"wall: serial {lean['wall_s']:.1f} s, 4 ranks {mpi['wall_s']:.1f} s")
    assert not s.failures(r), "\n".join(s.failures(r))
