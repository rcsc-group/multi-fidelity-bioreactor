"""Both SLURM templates chain segments: the next segment restarts from this one's dump in ITS OWN folder. The
dump's label (checkpoint.phase) must go with the dump, or a binary that reads the label refuses the restart
(src/forcing_state.h restart_source). Static checks: the templates are shell, run only by SLURM."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_mpi_template_copies_the_label_to_the_next_segment():
    txt = (ROOT / "config" / "slurm_mpi_template.sh").read_text()
    assert re.search(r'cp "\$\(dirname "\$CURR_CKPT"\)/checkpoint\.phase"\s*\\\s*"\$NEXT_SCRATCH/checkpoint\.phase"', txt)


def test_serial_template_copies_the_label_to_the_next_segment():
    txt = (ROOT / "config" / "slurm_template.sh").read_text()
    chain = txt[txt.index('if [ -n "$NEXT_RUN" ]'):]
    copy = re.search(r'cp "\$RUN_DIR/checkpoint\.phase" "\$\(dirname "\$NEXT_PARAMS"\)/checkpoint\.phase"', chain)
    assert copy, "the serial chain passes the dump but not its label"
    assert copy.start() < chain.index("sbatch"), "the label must be in place before the next segment is submitted"
