import json, math
R = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/ks_l8_rpm%s_th7/results.json"
S = {  # script seconds (RPM=32.5 in script): k1 plain, k1 strict, k2 plain, k2 strict
 "22.4": (74.3277, 71.0961, 75.3604, 71.6810),
 "22.6": (64.6328, 61.4011, 64.5946, 61.0515),
 "29.9": (float("nan"), 122.8023, 132.9449, 124.3897),
 "30.1": (93.7175, 87.2542, 94.7300, 86.4275)}
d = {}
for r, v in S.items():
    ours = json.load(open(R % r))["dtmix_0.95"]; sc = 32.5 / float(r)
    d[r] = [ours] + [x * sc for x in v]
    print(r, "ours %.2f | KP1 plain %.2f strict %.2f | KP2 plain %.2f strict %.2f (true s)" % tuple(d[r]))
names = ["ours", "K1 plain", "K1 strict", "K2 plain", "K2 strict"]
for a, b in (("22.4", "22.6"), ("29.9", "30.1")):
    print("pair", a, b)
    base = abs(math.log(d[b][0] / d[a][0]))
    for i, n in enumerate(names):
        x = abs(math.log(d[b][i] / d[a][i])); print("  %-10s |dlog|=%.3f  ratio to ours %.2f" % (n, x, x / base))
