"""Expected chance hit count for keyperm's models, and full-text quadgram score of each K4 hit."""
import pickle, numpy as np, collections, k4, cokernel, keyperm, quad, rk_mixed, gap2
from round3 import ALPHA
hits = pickle.load(open("keyperm_hits.pkl", "rb"))
M = keyperm.models(); exp = 0.0; D = {}
for name, Mx in M.items():
    N2, N13 = cokernel.left_null(Mx, 2), cokernel.left_null(Mx, 13)
    if min(len(N2), len(N13)) >= 5: exp += 13.0 ** -len(N13) * 2.0 ** -len(N2); D[name] = Mx
print(f"expected chance hits = 1,801,030 x 10 x sum = {1801030 * 10 * exp:,.0f}; observed {len(hits)}")
print("hits by L:", sorted(collections.Counter(int(h[-1].split('L=')[1]) for h in hits).items()))
# decrypt each hit fully: solve the key values from the cribs, fill the rest of the key where determined
AL = dict(rk_mixed.alph_list() + gap2.build()[0]); AL = {v: k for k, v in AL.items()}
POS = k4.CPOS.tolist(); best = []
for a, pn, cn, kn, name in hits:
    alph = AL[a]; PA = {"AZ": k4.AZ, "KA": k4.KA, "A": alph}[pn]; CA = {"AZ": k4.AZ, "KA": k4.KA, "A": alph}[cn]
    Mx = D[name]; w = int(name.split("w=")[1].split()[0]); L = int(name.split("L=")[1])
    # rebuild the full 97-position class map by re-running models() logic for all positions
    import itertools
    rows = -(-97 // w); cells = [(r, c) for r in range(rows) for c in range(w) if r * w + c < 97]
    seq = sorted(cells, key=lambda rc: (rc[1], rc[0])) if "R2L" not in name else sorted(cells, key=lambda rc: (-rc[1], rc[0]))
    cell_of_n = {n: r * w + c for n, (r, c) in enumerate(seq)}; n_of_cell = {v: k for k, v in cell_of_n.items()}
    cls = [(n_of_cell[p] if "inv" in name else cell_of_n[p]) % L for p in range(97)]
    key = {}
    for i, p in zip(POS, k4.CPL):
        c = CA.index(k4.CT[i]); q = PA.index(p); key[cls[i]] = (c - q) % 26 if kn == "vig" else (c + q) % 26
    out = []
    for i in range(97):
        if cls[i] in key:
            c = CA.index(k4.CT[i]); kv = key[cls[i]]; q = (c - kv) % 26 if kn == "vig" else (kv - c) % 26; out.append(PA[q])
        else: out.append(".")
    txt = "".join(out); known = txt.replace(".", ""); best.append((len(known), txt, a, pn, cn, kn, name))
best.sort(reverse=True)
for n, txt, a, pn, cn, kn, name in best[:12]: print(n, txt, a, pn, cn, kn, name)
