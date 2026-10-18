"""Vectorized two-stage solver across many alphabets (see twostage.py for the model).
stage1: xi = (PA(P) + s1*k1) if s1=+1 else (k1 - PA(P));  X = A1[xi]
stage2: bx = B1(X); C = A2[(bx + k2)] if s2=+1 else A2[(k2 - bx)]   (=> k2 = cidx - s2*bx ... see code)
Edge per crib j: k1[pos%L1] <-> k2[pos%L2].  forward: v = (cidx - s2*M[(a + s1*pidx)])  with M = B1inv o A1."""
import numpy as np, k4
POS = k4.CPOS.tolist(); PL = [ord(c) - 65 for c in k4.CPL]
AZ, KA = k4.AZ, k4.KA

def plan(L1, L2):
    """Return (nodes_total, list of cyclic components as op lists). node ids: ('1',r) / ('2',r)."""
    edges = [(("1", p % L1), ("2", p % L2)) for p in POS]
    nodes = sorted({u for e in edges for u in e})
    par = {u: u for u in nodes}
    def f(u):
        while par[u] != u: par[u] = par[par[u]]; u = par[u]
        return u
    for u, v in edges: par[f(u)] = f(v)
    comps = {}
    for j, (u, v) in enumerate(edges): comps.setdefault(f(u), []).append(j)
    out = []
    for js in comps.values():
        cn = {u for j in js for u in edges[j]}
        if len(js) <= len(cn) - 1: continue           # tree: always 26 solutions, no constraint
        best = None
        for root in cn:
            ops, cost = sim(root, js, edges)
            if best is None or cost < best[0]: best = (cost, root, ops)
        out.append((best[1], best[2], len(cn), len(js)))
    return len(nodes), out, len(comps)

def sim(root, js, edges):
    assigned = {root}; pend = set(js); ops = []; checks = 0; cost = 0.0
    while pend:
        ch = [j for j in pend if edges[j][0] in assigned and edges[j][1] in assigned]
        if ch:
            for j in sorted(ch): ops.append(("chk", j)); pend.discard(j)
            checks += len(ch); continue
        cands = []
        for j in pend:
            u, v = edges[j]
            if (u in assigned) != (v in assigned):
                new = v if u in assigned else u
                score = sum(1 for k in pend if new in edges[k] and (edges[k][0] in assigned or edges[k][1] in assigned))
                cands.append((-score, j, new))
        cands.sort(); _, j, new = cands[0]
        ops.append(("prop", j, new[0] == "2"))           # True: forward k1->k2
        pend.discard(j)
        assigned.add(new); cost += 26.0 ** (-checks)
    return ops, cost

def inv(A):   # A (n,26) letter idx at position -> position of letter
    n = A.shape[0]; out = np.empty_like(A); out[np.arange(n)[:, None], A] = np.arange(26, dtype=A.dtype)[None, :]; return out

def prep(PA, A1, B1, A2, ct):
    """each arg (n,26) int arrays of letter indices; returns pidx(n,24), cidx(n,24), M(n,26), Minv(n,26)."""
    PAi, B1i, A2i, A1i = inv(PA), inv(B1), inv(A2), inv(A1)
    ctl = np.array([ord(ct[p]) - 65 for p in POS])
    pidx = PAi[:, PL]; cidx = A2i[:, ctl]
    M = np.take_along_axis(B1i, A1.astype(np.int64), 1); Minv = np.take_along_axis(A1i, B1.astype(np.int64), 1)
    return pidx.astype(np.int16), cidx.astype(np.int16), M.astype(np.int16), Minv.astype(np.int16)

def fwd(a, rows, j, P, s1, s2):
    pidx, cidx, M, Minv = P
    xi = (a + s1 * pidx[rows, j]) % 26; bx = M[rows, xi]; return (cidx[rows, j] - s2 * bx) % 26
def bwd(v, rows, j, P, s1, s2):
    pidx, cidx, M, Minv = P
    bx = (s2 * (cidx[rows, j] - v)) % 26; xi = Minv[rows, bx]; return (xi - s1 * pidx[rows, j]) % 26

def solve_cfg(P, n, L1, L2, s1, s2, pl):
    """returns sols array (n,) int64 = number of consistent assignments per alphabet (0 = inconsistent)."""
    nodes, comps, ncomp = pl
    edges = [(("1", p % L1), ("2", p % L2)) for p in POS]
    sols = np.ones(n, dtype=np.int64) * (26 ** (ncomp - len(comps)))
    alive = np.arange(n)
    for root, ops, V, E in comps:
        if alive.size == 0: break
        rows = np.repeat(alive, 26); vals = {root: np.tile(np.arange(26, dtype=np.int16), alive.size)}
        for op in ops:
            if rows.size == 0: break
            if op[0] == "prop":
                _, j, forward = op; u, v = edges[j]
                if forward: vals[v] = fwd(vals[u], rows, j, P, s1, s2)
                else: vals[u] = bwd(vals[v], rows, j, P, s1, s2)
            else:
                j = op[1]; u, v = edges[j]
                m = fwd(vals[u], rows, j, P, s1, s2) == vals[v]
                rows = rows[m]; vals = {k: x[m] for k, x in vals.items()}
        cnt = np.bincount(rows, minlength=n)
        sols *= cnt
        alive = np.nonzero(sols)[0]
    sols[np.setdiff1d(np.arange(n), alive)] = 0
    return sols

def to_arr(strs):
    return np.frombuffer("".join(strs).encode(), dtype=np.uint8).reshape(len(strs), 26).astype(np.int16) - 65

PLANS = {}
def get_plan(L1, L2):
    if (L1, L2) not in PLANS: PLANS[(L1, L2)] = plan(L1, L2)
    return PLANS[(L1, L2)]

def run_chunk(args):
    """args: (ct, role, X, alph_strs, names, Lmax, netmin). role 'i': PA=A1=X, B1=A2=D ; 'ii': PA=A1=D, B1=A2=X."""
    ct, role, X, strs, names, Lmax, netmin = args
    D = to_arr(strs); n = len(strs); Xa = np.broadcast_to(to_arr([X]), D.shape)
    P = prep(Xa, Xa, D, D, ct) if role == "i" else prep(D, D, Xa, Xa, ct)
    out = []; ncfg = 0
    for L1 in range(1, Lmax + 1):
        for L2 in range(1, Lmax + 1):
            pl = get_plan(L1, L2); net = 24 - pl[0]
            if net < netmin: continue
            for s1 in (1, -1):
                for s2 in (1, -1):
                    ncfg += n
                    sols = solve_cfg(P, n, L1, L2, s1, s2, pl)
                    for a in np.nonzero(sols)[0]:
                        out.append((net, int(sols[a]), names[a], strs[a], L1, L2, s1, s2))
    return out, ncfg

# ---------- scalar reconstruction ----------
def assignments(ct, L1, L2, PA, A1, B1, A2, s1, s2):
    """all components (cyclic and tree) with list of consistent assignments dict node->value."""
    PAi = {c: i for i, c in enumerate(PA)}; A1i = {c: i for i, c in enumerate(A1)}; B1i = {c: i for i, c in enumerate(B1)}; A2i = {c: i for i, c in enumerate(A2)}
    adj = {}
    for p, pc in zip(POS, k4.CPL):
        pidx, cidx = PAi[pc], A2i[ct[p]]
        f = lambda a, pidx=pidx, cidx=cidx: (cidx - s2 * B1i[A1[(a + s1 * pidx) % 26]]) % 26
        b = lambda v, pidx=pidx, cidx=cidx: (A1i[B1[(s2 * (cidx - v)) % 26]] - s1 * pidx) % 26
        u, v = ("1", p % L1), ("2", p % L2)
        adj.setdefault(u, []).append((v, f)); adj.setdefault(v, []).append((u, b))
    seen = set(); comps = []
    for root in adj:
        if root in seen: continue
        good = []; cset = None
        for val in range(26):
            asg = {root: val}; st = [root]; ok = True
            while st and ok:
                u = st.pop()
                for w, fn in adj[u]:
                    x = fn(asg[u])
                    if w in asg:
                        if asg[w] != x: ok = False; break
                    else: asg[w] = x; st.append(w)
            if cset is None and ok: cset = set(asg)
            if ok: good.append(asg)
        if cset is None:   # inconsistent; find node set by BFS
            cset = set(); st = [root]
            while st:
                u = st.pop()
                if u in cset: continue
                cset.add(u); st += [w for w, _ in adj[u]]
        seen |= cset; comps.append(good)
    return comps

def decrypt(ct, L1, L2, PA, A1, B1, A2, s1, s2, k1, k2):
    out = []
    for i, c in enumerate(ct):
        a, v = k1.get(i % L1), k2.get(i % L2)
        if a is None or v is None: out.append("?"); continue
        cidx = A2.index(c); bx = (s2 * (cidx - v)) % 26; xi = A1.index(B1[bx]); p = (s1 * (xi - a)) % 26
        out.append(PA[p])
    return "".join(out)

def encrypt(pt, L1, L2, PA, A1, B1, A2, s1, s2, k1, k2):
    out = []
    for i, p in enumerate(pt):
        a, v = k1[i % L1], k2[i % L2]
        xi = (PA.index(p) + a) % 26 if s1 == 1 else (a - PA.index(p)) % 26
        bx = B1.index(A1[xi]); ci = (bx + v) % 26 if s2 == 1 else (v - bx) % 26
        out.append(A2[ci])
    return "".join(out)
