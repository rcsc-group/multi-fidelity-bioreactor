#!/bin/bash
# Daily CI watcher, run by a SessionStart hook with async + asyncRewake.
#
# Silent while CI on main is green. Exits 2 (which wakes Claude) the first time
# it sees a failed run it has not reported yet; the message goes to stderr.
# Exits 0 without doing anything if another watcher is already alive on this
# host. It dies with the session or the node, which is intended.
#
#   CI_WATCH_INTERVAL  seconds between checks (default 86400)
#   CI_WATCH_STATE     state directory (default ~/.cache/ci_watch)
set -u
REPO_DIR=/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor
INTERVAL=${CI_WATCH_INTERVAL:-86400}
STATE=${CI_WATCH_STATE:-$HOME/.cache/ci_watch}
mkdir -p "$STATE"
LOCK="$STATE/watcher.$(hostname -s).pid"

if [ -f "$LOCK" ] && kill -0 "$(cat "$LOCK")" 2>/dev/null; then
    exit 0
fi
echo $$ > "$LOCK"
trap 'rm -f "$LOCK"' EXIT

cd "$REPO_DIR" || exit 0
while true; do
    # latest COMPLETED run of the CI workflow on main
    line=$(gh run list -w CI -b main -s completed -L 1 \
             --json databaseId,conclusion,headSha,displayTitle,url \
             -q '.[0] | "\(.databaseId) \(.conclusion) \(.headSha[0:7]) \(.url) \(.displayTitle)"' 2>/dev/null)
    if [ -n "$line" ]; then
        read -r id concl sha url title <<< "$line"
        if [ "$concl" = "failure" ] && [ "$id" != "$(cat "$STATE/last_red" 2>/dev/null)" ]; then
            echo "$id" > "$STATE/last_red"
            failed=$(gh run view "$id" --log-failed 2>/dev/null | grep -E "FAILED " | sed 's/.*FAILED /FAILED /' | sort -u | head -5)
            echo "CI is RED on main: run $id ($sha, \"$title\") $url" >&2
            [ -n "$failed" ] && echo "$failed" >&2
            exit 2
        fi
    fi
    sleep "$INTERVAL"
done
