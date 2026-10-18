"""Sanborn-phrase alphabets (613k) x all linear key models (repeating, progressive, grid, combined, formula; >=8 checks)
and x autokey (both sides, lags 1-96). Fills the last cheap cells of the coverage matrix."""
import numpy as np, time, k4, gap1, phrases
from gromark_dict import straight, transposed
PH = phrases.phrase_list(); al = {}
for p in PH:
    for f, a in (("s", straight(p)), ("t", transposed(p))):
        al.setdefault(a, f"{p}/{f}"); al.setdefault(a[::-1], f"{p}/{f}-rev")
AL = list(al.items()); print(f"phrase alphabets: {len(AL):,}; linear models: {len(gap1.P)}", flush=True)
# planted check: PALIMPSEST/s alphabet, L=9 repeating key, vig
import random; random.seed(3)
pa = straight("PALIMPSEST"); key = [random.randrange(26) for _ in range(9)]
pt = list(k4.CT)
for i, p in zip(k4.CPOS.tolist(), k4.CPL): pt[i] = p
ct = "".join(pa[(pa.index(c) + key[i % 9]) % 26] for i, c in enumerate(pt))  # Q3-style vig, same alphabet both sides
h = gap1.scan(ct, [(pa, "PLANT")], "plant"); print("planted:", [x for x in h if "L=9" in x[-1]][:3], flush=True)
t = time.time(); h = gap1.scan(k4.CT, AL, "K4"); print(f"K4 linear models: {len(h)} hits ({time.time()-t:.0f}s)", flush=True)
for x in h[:20]: print("  ", x, flush=True)
rng = np.random.default_rng(7)
for i in range(2):
    hn = gap1.scan("".join(rng.permutation(list(k4.CT))), AL, "null"); print(f"shuffled {i}: {len(hn)} hits", flush=True)
# autokey: reuse autokey_dict machinery if importable
try:
    import autokey_dict as ak
    t = time.time(); r = ak.run(AL) if hasattr(ak, "run") else None
    print("autokey:", r if r is not None else "no run() entry point; skipped", f"({time.time()-t:.0f}s)", flush=True)
except Exception as e:
    print("autokey skipped:", repr(e), flush=True)
