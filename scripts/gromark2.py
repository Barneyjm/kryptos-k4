"""Gromark with BOTH alphabets unknown. Crib i: c(C_i) - p(P_i) = k_i (mod 26).
Graph: plain-letter nodes and cipher-letter nodes, edge per crib. Potentials along a spanning tree are linear in k.
Survive iff (a) every non-tree edge closes consistently, (b) same-side nodes in one component get distinct values."""
import numpy as np, itertools, sys, time, multiprocessing as mp, k4
POS = k4.CPOS.tolist(); PL = k4.CPL; CL = [k4.CT[i] for i in POS]
nodes = sorted({"p" + c for c in PL} | {"c" + c for c in CL}); NI = {n: i for i, n in enumerate(nodes)}
edges = [(NI["p" + PL[j]], NI["c" + CL[j]], j) for j in range(24)]          # c - p = k_j
adj = {i: [] for i in range(len(nodes))}
for a, b, j in edges: adj[a].append((b, j, +1)); adj[b].append((a, j, -1))
pot = {}; comp = {}; tree = set()
for s in range(len(nodes)):
    if s in pot: continue
    pot[s] = np.zeros(24, int); comp[s] = s; stack = [s]
    while stack:
        u = stack.pop()
        for v, j, sg in adj[u]:
            if v not in pot:
                pot[v] = pot[u].copy(); pot[v][j] += sg; comp[v] = s; tree.add(j); stack.append(v)
cyc = []          # non-tree edges: pot[c] - pot[p] - k_j == 0
for a, b, j in edges:
    if j not in tree:
        v = pot[b] - pot[a]; v = v.copy(); v[j] -= 1; cyc.append(v)
dist = []         # same-side nodes in same component must differ
for x, y in itertools.combinations(range(len(nodes)), 2):
    if comp[x] == comp[y] and nodes[x][0] == nodes[y][0]: dist.append(pot[x] - pot[y])
CYC = np.array(cyc).T if cyc else np.zeros((24, 0), int); DIST = np.array(dist).T
def expand(prim, b, lags, upto=97):
    m, n = prim.shape; d = np.zeros((m, upto), np.int64); d[:, :n] = prim; p, q = lags
    for i in range(n, upto): d[:, i] = (d[:, i - p] + d[:, i - q]) % b
    return d
def survivors(args):
    b, n, lags = args; total = b ** n; out = []
    for start in range(0, total, 1_000_000):
        idx = np.arange(start, min(total, start + 1_000_000))
        prim = np.stack([(idx // b ** (n - 1 - t)) % b for t in range(n)], axis=1)
        K = expand(prim, b, lags, 74)[:, POS]
        ok = np.all((K @ CYC) % 26 == 0, axis=1) & np.all((K @ DIST) % 26 != 0, axis=1)
        out += ["".join(map(str, prim[i])) if b <= 10 else tuple(prim[i].tolist()) for i in np.nonzero(ok)[0]]
    return (b, n, lags, out)
if __name__ == "__main__":
    print(f"graph: {len(nodes)} letter-nodes, {len(cyc)} loop checks, {len(dist)} distinctness checks", flush=True)
    b10 = survivors((10, 5, (4, 5))); print("base 10, 5-digit primer, standard rule:", len(b10[3]), "survivors (Bean: 39); 26717 in?", "26717" in b10[3], "84393 in?", "84393" in b10[3], flush=True)
    b8 = survivors((8, 5, (4, 5))); print("base 8: includes 00351?", "00351" in b8[3], "00537?", "00537" in b8[3], "count", len(b8[3]), flush=True)
