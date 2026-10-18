"""Two-stage polyalphabetic: X_i = A1[(PA(P_i) ± k1[i%L1])], C_i = A2[(B1(X_i) ± k2[i%L2])]. Nonlinear when A1 != B1.
Each crib is a bijection between k1[r1] and k2[r2]; solve per connected component (26 candidates each)."""
import itertools, numpy as np, time, multiprocessing as mp, k4
from round3 import ALPHA
CR = list(zip(k4.CPOS.tolist(), k4.CPL))
def count(ct, L1, L2, PA, A1, B1, A2, s1, s2):
    PAi = {c: i for i, c in enumerate(PA)}; A1i = {c: i for i, c in enumerate(A1)}; B1i = {c: i for i, c in enumerate(B1)}; A2i = {c: i for i, c in enumerate(A2)}
    edges = []
    for i, p in CR:
        pidx, cidx = PAi[p], A2i[ct[i]]
        def fwd(a, pidx=pidx, cidx=cidx):
            xi = (pidx + a) % 26 if s1 == 1 else (a - pidx) % 26; bx = B1i[A1[xi]]
            return (cidx - bx) % 26 if s2 == 1 else (cidx + bx) % 26
        def bwd(v, pidx=pidx, cidx=cidx):
            bx = (cidx - v) % 26 if s2 == 1 else (v - cidx) % 26; xi = A1i[B1[bx]]
            return (xi - pidx) % 26 if s1 == 1 else (xi + pidx) % 26
        edges.append((("k1", i % L1), ("k2", i % L2), fwd, bwd))
    adj = {}
    for u, v, f, b in edges: adj.setdefault(u, []).append((v, f)); adj.setdefault(v, []).append((u, b))
    seen = set(); total = 1; nodes = len(adj)
    for root in adj:
        if root in seen: continue
        comp_ok = 0; comp = set()
        for val in range(26):
            assign = {root: val}; stack = [root]; ok = True
            while stack and ok:
                u = stack.pop()
                for w, f in adj[u]:
                    x = f(assign[u])
                    if w in assign:
                        if assign[w] != x: ok = False; break
                    else: assign[w] = x; stack.append(w)
            comp = set(assign) if not comp else comp
            if ok: comp_ok += 1
        seen |= comp if comp else {root}
        if not comp:  # isolated (cannot happen: every node has an edge)
            continue
        total *= comp_ok
        if total == 0: return 0, nodes
    return total, nodes
def job(args):
    ct, combos = args; out = []
    for (pn, PA), (a1n, A1), (b1n, B1), (a2n, A2), s1, s2, L1, L2 in combos:
        sols, nodes = count(ct, L1, L2, PA, A1, B1, A2, s1, s2)
        if sols: 
            net = 24 - nodes                      # independent checks beyond the unknowns
            if net >= 4: out.append((net, sols, pn, a1n, b1n, a2n, s1, s2, L1, L2))
    return out
if __name__ == "__main__":
    std = [("AZ", k4.AZ), ("KA", k4.KA)]
    combos = [(p, a1, b1, a2, s1, s2, L1, L2) for p in std for a1 in std for b1 in std for a2 in std
              for s1 in (1, -1) for s2 in (1, -1) for L1 in range(1, 27) for L2 in range(1, 27)]
    t = time.time(); r = job((k4.CT, combos)); print(f"std/KA all slots, L1,L2<=26: {len(combos):,} configs, significant solutions: {len(r)} ({time.time()-t:.0f}s)", flush=True)
    for x in sorted(r, reverse=True)[:10]: print("  ", x, flush=True)
    rng = np.random.default_rng(0); sh = "".join(rng.permutation(list(k4.CT))); rs = job((sh, combos)); print("shuffled baseline:", len(rs), flush=True)
    # stage 1 = Quagmire on alphabet X (PA=A1=X), stage 2 on a different alphabet Y (B1=A2=Y); clue-word alphabets
    items = list(ALPHA.items()); pairs = [(x, y) for x in items for y in items if x[1] != y[1]]
    combos2 = [(x, x, y, y, s1, s2, L1, L2) for x, y in pairs for s1 in (1, -1) for s2 in (1, -1) for L1 in range(1, 11) for L2 in range(1, 11)]
    print(f"clue alphabets, two different tableaux, L1,L2<=10: {len(combos2):,} configs", flush=True)
    chunks = [combos2[i:i + 20000] for i in range(0, len(combos2), 20000)]
    t = time.time()
    with mp.Pool(6) as pool: res = [x for r in pool.imap_unordered(job, [(k4.CT, c) for c in chunks]) for x in r]
    print(f"significant solutions: {len(res)} ({time.time()-t:.0f}s)", flush=True)
    for x in sorted(res, reverse=True)[:15]: print("  ", x, flush=True)
