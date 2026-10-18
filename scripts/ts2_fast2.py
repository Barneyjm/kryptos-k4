"""ts2_fast with flat dense forward AND backward maps (single gathers)."""
import numpy as np, ts2_vec as V
POS = V.POS
EDGES = {}
def edges_of(L1, L2):
    if (L1, L2) not in EDGES: EDGES[(L1, L2)] = [(("1", p % L1), ("2", p % L2)) for p in POS]
    return EDGES[(L1, L2)]

def maps(P, s1, s2):
    pidx, cidx, M, Minv = P
    a = np.arange(26, dtype=np.int64)[None, :]
    F, G = [], []
    for j in range(24):
        F.append(((cidx[:, j:j + 1] - s2 * np.take_along_axis(M, (a + s1 * pidx[:, j:j + 1]) % 26, 1)) % 26).astype(np.int64).ravel())
        G.append(((np.take_along_axis(Minv, (s2 * (cidx[:, j:j + 1] - a)) % 26, 1) - s1 * pidx[:, j:j + 1]) % 26).astype(np.int64).ravel())
    return F, G

def solve_cfg(F, G, n, L1, L2, pl):
    nodes, comps, ncomp = pl
    edges = edges_of(L1, L2)
    sols = np.full(n, 26 ** (ncomp - len(comps)), dtype=np.int64)
    alive_mask = np.ones(n, bool)
    for root, ops, Vn, E in comps:
        if not alive_mask.any(): break
        if len(ops) >= 2 and ops[0][0] == "prop" and ops[1][0] == "chk" and set(edges[ops[0][1]]) == set(edges[ops[1][1]]):
            j, j2 = ops[0][1], ops[1][1]
            m = (F[j] == F[j2]); m &= np.repeat(alive_mask, 26)
            flat = np.nonzero(m)[0]; base = flat - flat % 26
            u, v = edges[j]; vals = {u: flat % 26, v: F[j][flat]}; rest = ops[2:]
        else:
            alive = np.nonzero(alive_mask)[0]
            base = np.repeat(alive * 26, 26); vals = {root: np.tile(np.arange(26), alive.size)}; rest = ops
        for op in rest:
            if base.size == 0: break
            if op[0] == "prop":
                _, j, forward = op; u, v = edges[j]
                if forward: vals[v] = F[j][base + vals[u]]
                else: vals[u] = G[j][base + vals[v]]
            else:
                j = op[1]; u, v = edges[j]
                m = F[j][base + vals[u]] == vals[v]
                base = base[m]; vals = {k: x[m] for k, x in vals.items()}
        sols *= np.bincount(base // 26, minlength=n)
        alive_mask = sols > 0
    return sols

def run_chunk(args):
    """(ct, role, X, strs, names, Lmax, netmin); role 'i': PA=A1=X,B1=A2=D; 'ii': PA=A1=D,B1=A2=X."""
    ct, role, X, strs, names, Lmax, netmin = args
    D = V.to_arr(strs); n = len(strs); Xa = np.broadcast_to(V.to_arr([X]), D.shape)
    P = V.prep(Xa, Xa, D, D, ct) if role == "i" else V.prep(D, D, Xa, Xa, ct)
    P = tuple(x.astype(np.int64) for x in P)
    out = []; ncfg = 0
    for s1 in (1, -1):
        for s2 in (1, -1):
            F, G = maps(P, s1, s2)
            for L1 in range(1, Lmax + 1):
                for L2 in range(1, Lmax + 1):
                    pl = V.get_plan(L1, L2); net = 24 - pl[0]
                    if net < netmin: continue
                    ncfg += n
                    sols = solve_cfg(F, G, n, L1, L2, pl)
                    for a in np.nonzero(sols)[0]:
                        out.append((net, int(sols[a]), names[a], strs[a], L1, L2, s1, s2))
    return out, ncfg
