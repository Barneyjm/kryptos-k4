"""Running-key test: key = some text at some offset (plain, reversed). Exact match on all 24 crib letters."""
import numpy as np, itertools, sys, k4
from texts import *
SOURCES = {"K1":K1,"K2":K2,"K2orig":K2_ORIG_END,"K3":K3,"K123":K1+K2+K3,"K123orig":K1+K2_ORIG_END+K3,
           "PANEL_CT":PANEL,"K3K2K1":K3+K2+K1}
for n,s in list(SOURCES.items()): SOURCES[n+"_rev"] = s[::-1]
def test(dest, label, hits):
    dc = dest[k4.CPOS]
    for an, X in k4.ALPHS.items():
        ci = np.array([X.index(c) for c in k4.CT])[dc]
        pi = np.array([X.index(c) for c in k4.CPL])
        for kn, need in (("vig",(ci-pi)%26),("beau",(ci+pi)%26),("vbeau",(pi-ci)%26)):
            for sn, s in SOURCES.items():
                kidx = np.array([X.index(c) for c in s])
                for model in "AB":
                    pos = k4.CPOS if model=="A" else dc
                    # key[t] = s[off + t]; need kidx[off+pos] == need for all cribs
                    maxoff = len(s) - 97
                    for off in range(-96, max(maxoff,0)+1):
                        p = off + pos
                        if p.min() < 0 or p.max() >= len(s): continue
                        m = (kidx[p] == need).sum()
                        if m >= 14: hits.append((label, an, kn, sn, model, off, int(m)))
hits = []
test(np.arange(97), "none", hits)
print("no-transposition hits:", hits)
W = int(sys.argv[1]) if len(sys.argv)>1 else 6
for w in range(2, W+1):
    for p in itertools.permutations(range(w)):
        order = np.array(k4.order_from(w, p)); dest = np.empty(97,int); dest[order]=np.arange(97)
        test(dest, f"w{w} enc {p}", hits); test(order, f"w{w} dec {p}", hits)
    print("w",w,"hits",len(hits), file=sys.stderr, flush=True)
for h in sorted(hits, key=lambda h:-h[-1])[:20]: print(h)
