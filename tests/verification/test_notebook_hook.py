"""scripts/hooks/pre-commit must block code commits when the notebook is stale.

Each case builds a throwaway repo, stages a change, and runs the hook the way
git would. Ages are set with os.utime / GIT_COMMITTER_DATE rather than waited
for, so every branch of the hook is exercised in well under a second.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path

import pytest

HOOK = Path(__file__).resolve().parents[2] / "scripts" / "hooks" / "pre-commit"
HOUR, DAY = 3600, 86400


def git(repo: Path, *args: str, when: float | None = None) -> None:
    env = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}
    if when is not None:
        stamp = f"@{int(when)} +0000"
        env |= {"GIT_AUTHOR_DATE": stamp, "GIT_COMMITTER_DATE": stamp}
    subprocess.run(["git", *args], cwd=repo, env=env, check=True,
                   capture_output=True)


@pytest.fixture
def repo(tmp_path: Path):
    """Repo whose last commit (and BACKLOG.md) is `commit_age` old."""
    def make(commit_age: float = 3 * HOUR, backlog_age: float | None = None):
        r = tmp_path / "r"
        r.mkdir()
        git(r, "init", "-q")
        (r / "scripts").mkdir()
        (r / "scripts" / "a.py").write_text("x = 1\n")
        (r / "BACKLOG.md").write_text("# Backlog\n")
        (r / ".gitignore").write_text("diary.md\n")
        now = time.time()
        git(r, "add", "-A")
        git(r, "commit", "-qm", "init",
            when=now - (backlog_age if backlog_age is not None else commit_age))
        (r / "docs.md").write_text("d\n")
        git(r, "add", "docs.md")
        git(r, "commit", "-qm", "later", when=now - commit_age)
        (r / "diary.md").write_text("notes\n")
        return r
    return make


def run_hook(r: Path, **env: str) -> subprocess.CompletedProcess:
    return subprocess.run([str(HOOK)], cwd=r, capture_output=True, text=True,
                          env={**os.environ, **env})


def age(path: Path, seconds: float) -> None:
    t = time.time() - seconds
    os.utime(path, (t, t))


def stage_code(r: Path) -> None:
    (r / "scripts" / "a.py").write_text("x = 2\n")
    git(r, "add", "scripts/a.py")


def test_hook_is_executable():
    assert HOOK.exists() and os.access(HOOK, os.X_OK)


def test_stale_diary_blocks_a_code_commit(repo):
    r = repo()
    age(r / "diary.md", 5 * HOUR)          # older than HEAD (3 h) and > 2 h
    stage_code(r)
    res = run_hook(r)
    assert res.returncode == 1
    assert "diary.md" in res.stderr


def test_diary_edited_since_last_commit_passes(repo):
    r = repo()
    age(r / "diary.md", 1 * HOUR)          # newer than HEAD
    stage_code(r)
    assert run_hook(r).returncode == 0


def test_diary_within_grace_window_passes_even_before_head(repo):
    r = repo(commit_age=10 * 60)           # HEAD 10 min ago
    age(r / "diary.md", 30 * 60)           # before HEAD, but only 30 min old
    stage_code(r)
    assert run_hook(r).returncode == 0


def test_missing_diary_blocks(repo):
    r = repo()
    (r / "diary.md").unlink()
    stage_code(r)
    assert run_hook(r).returncode == 1


def test_docs_only_commit_is_not_checked(repo):
    r = repo()
    age(r / "diary.md", 50 * HOUR)
    (r / "docs.md").write_text("changed\n")
    git(r, "add", "docs.md")
    assert run_hook(r).returncode == 0


def test_stale_backlog_blocks(repo):
    r = repo(backlog_age=10 * DAY)
    age(r / "diary.md", 60)
    stage_code(r)
    res = run_hook(r)
    assert res.returncode == 1
    assert "BACKLOG.md" in res.stderr


def test_staging_the_backlog_satisfies_it(repo):
    r = repo(backlog_age=10 * DAY)
    age(r / "diary.md", 60)
    stage_code(r)
    (r / "BACKLOG.md").write_text("# Backlog\n- new item\n")
    git(r, "add", "BACKLOG.md")
    assert run_hook(r).returncode == 0


def test_override_passes_and_says_why(repo):
    r = repo(backlog_age=10 * DAY)
    age(r / "diary.md", 50 * HOUR)
    stage_code(r)
    res = run_hook(r, NOTEBOOK_OK="typo fix")
    assert res.returncode == 0
    assert "typo fix" in res.stderr


@pytest.mark.skipif(shutil.which("git") is None, reason="needs git")
def test_installed_copy_matches_tracked_hook():
    root = HOOK.parents[2]
    installed = subprocess.run(["git", "rev-parse", "--git-path", "hooks"],
                               cwd=root, capture_output=True, text=True).stdout.strip()
    target = (root / installed / "pre-commit")
    if not target.exists():
        pytest.skip("hook not installed in this clone (run `make hooks`)")
    assert target.read_bytes() == HOOK.read_bytes(), "re-run `make hooks`"
