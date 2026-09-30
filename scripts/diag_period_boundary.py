"""Where does build/BioReactor stop for a request of exactly k periods?

Same cases as tests/verification/test_period_boundary.py, run directly
(medium tests are run by CI, not on OSCAR). Prints the dump time in periods.
Usage:  uv run python scripts/diag_period_boundary.py OUTDIR
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.conftest import run_bioreactor                                 # noqa: E402
from tests.verification.test_restart_continue import _params, T_PER_ND   # noqa: E402
from scripts.dump_fields import fields as dump_fields                     # noqa: E402

out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
for cycles in (10, 10.4, 20):
    d = run_bioreactor(_params(f"pb_{cycles:g}", cycles), out, timeout=900)
    t = dump_fields(d / "checkpoint.dump")[0]["t"]
    print(f"requested {cycles:5g} periods -> dumped at {t / T_PER_ND:.6f} periods", flush=True)
