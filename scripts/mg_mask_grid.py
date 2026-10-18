"""Arbitrary mask + structured keys (extends mask.py). Exact linear algebra over Z2 and Z13 (=> Z26 by CRT).
plain side : x[P_i] + z_i = num(C_i)      (x = unknown injective mask on the 13 crib plaintext letters)
cipher side: y[C_i] - z_i = num(P_i)      (y = unknown injective mask on the crib ciphertext letters)
(vig/beau/variant are absorbed by negating the arbitrary mask/key; numbering AZ or KA of the fixed side matters.)
Key families:
 grid    z_i = r[i%W] + c[i//W]                     W=2..30
 gridH   z_i = r[i%W] + c[(i//W)%H]                 W=2..30, H=2..13
 comb    z_i = u[i%L1] + v[i%L2]                    2<=L1<L2<=26
 prog    z_i = u[i%L] + s*(i//L)   (s unknown)      L=2..26
 kgrid   z_i = u[f(i)%L], f = keyword written in w-wide grid row-wise, read col-wise (keyperm.py), 4 variants, w=2..30, L=2..26
Gauge: the shift x+t, z-t is fixed by x[first letter]=0. checks = 24 - rank(24 crib rows) per prime.
Classes: inconsistent | noinj (consistent but no injective mask) | inj (consistent, has checks, injective mask exists)
         | trivial (0 checks at both primes; injective exists -> undecidable) | skip (x-nullity > 6 -> undecidable)."""
import sys, itertools, random, numpy as np, k4
POS = k4.CPOS.tolist(); CR = dict(zip(POS, k4.CPL)); N = 97

def rref(A, p):
    A = A.copy() % p; m, n = A.shape; piv = []; r = 0
    for c in range(n):
        if r == m: break
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        pr = r + nz[0]; A[[r, pr]] = A[[pr, r]]
        A[r] = A[r] * pow(int(A[r, c]), -1, p) % p
        f = A[:, c].copy(); f[r] = 0
        A = (A - np.outer(f, A[r])) % p
        piv.append(c); r += 1
    return A, piv

def solve(A, b, p):
    """A x = b mod p -> (consistent, particular, nullbasis (k,n), rank of A)."""
    m, n = A.shape
    M, piv = rref(np.hstack([A, b[:, None]]) % p, p)
    if n in piv: return False, None, None, len(piv) - 1
    x = np.zeros(n, np.int64)
    for i, c in enumerate(piv): x[c] = M[i, n]
    free = [c for c in range(n) if c not in piv]; B = np.zeros((len(free), n), np.int64)
    for k, f in enumerate(free):
        B[k, f] = 1
        for i, c in enumerate(piv): B[k, c] = (-M[i, f]) % p
    return True, x, B, len(piv)

def proj_basis(B, p):
    if len(B) == 0: return B
    M, piv = rref(B, p); return M[:len(piv)]

def points(x0, B, p):
    """all points of affine space x0 + span(B) mod p, shape (p^k, nx)."""
    k = len(B)
    if k == 0: return x0[None, :] % p
    coef = np.array(list(itertools.product(range(p), repeat=k)), np.int64)
    return (x0[None, :] + coef @ B) % p

def injective_exists(x2, B2, x13, B13):
    nx = len(x2); I, J = np.triu_indices(nx, 1)
    # 1) a pair forced equal at both primes on the whole solution space -> no injective mask
    f2 = ((x2[I] - x2[J]) % 2 == 0) & ((B2[:, I] - B2[:, J]) % 2 == 0).all(0) if len(B2) else (x2[I] - x2[J]) % 2 == 0
    f13 = ((x13[I] - x13[J]) % 13 == 0) & ((B13[:, I] - B13[:, J]) % 13 == 0).all(0) if len(B13) else (x13[I] - x13[J]) % 13 == 0
    if (f2 & f13).any(): return False
    # 2) random sampling (fast witness)
    rs = np.random.default_rng(1); n = 4096
    for _ in range(5):
        P2 = (x2[None] + (rs.integers(0, 2, (n, len(B2))) @ B2 if len(B2) else 0)) % 2
        P13 = (x13[None] + (rs.integers(0, 13, (n, len(B13))) @ B13 if len(B13) else 0)) % 13
        if (~((P2[:, I] == P2[:, J]) & (P13[:, I] == P13[:, J])).any(1)).any(): return True
    # 3) exhaustive
    P2 = points(x2, B2, 2); E2 = np.unique(P2[:, I] == P2[:, J], axis=0)
    for s in range(0, 13 ** len(B13), 1 << 20):            # chunk mod-13 enumeration
        k = len(B13)
        if k == 0: P13 = x13[None, :] % 13
        else:
            idx = np.arange(s, min(s + (1 << 20), 13 ** k)); coef = np.stack([(idx // 13 ** t) % 13 for t in range(k)], 1)
            P13 = (x13[None, :] + coef @ B13) % 13
        E13 = P13[:, I] == P13[:, J]
        for e2 in E2:
            if (~(E13 & e2[None, :]).any(1)).any(): return True
    return False

def grid_index(w):
    rows = -(-N // w); cells = [(r, c) for r in range(rows) for c in range(w) if r * w + c < N]
    out = {}
    for name, seq in (("L2R", sorted(cells, key=lambda rc: (rc[1], rc[0]))), ("R2L", sorted(cells, key=lambda rc: (-rc[1], rc[0])))):
        cell_of_n = [r * w + c for (r, c) in seq]; n_of_cell = [0] * N
        for n_, cc in enumerate(cell_of_n): n_of_cell[cc] = n_
        out[name] = cell_of_n; out[name + " inv"] = n_of_cell
    return out

def models():
    """yield (family, name, keyfun) ; keyfun(i) -> list of (key component label, integer coefficient)."""
    for W in range(2, 31): yield "grid", f"W={W}", (lambda i, W=W: [(("r", i % W), 1), (("c", i // W), 1)])
    for W in range(2, 31):
        for H in range(2, 14): yield "gridH", f"W={W} H={H}", (lambda i, W=W, H=H: [(("r", i % W), 1), (("c", (i // W) % H), 1)])
    for L1 in range(2, 27):
        for L2 in range(L1 + 1, 27): yield "comb", f"L1={L1} L2={L2}", (lambda i, a=L1, b=L2: [(("u", i % a), 1), (("v", i % b), 1)])
    for L in range(2, 27): yield "prog", f"L={L}", (lambda i, L=L: [(("u", i % L), 1)] + ([(("s",), i // L)] if i // L else []))
    for w in range(2, 31):
        G = grid_index(w)
        for vn, f in G.items():
            for L in range(2, 27): yield "kgrid", f"w={w} {vn} L={L}", (lambda i, f=f, L=L: [(("u", f[i] % L), 1)])

def build(ct, side, X, keyfun):
    letters = sorted({(CR[i] if side == "plain" else ct[i]) for i in POS}); LI = {c: j for j, c in enumerate(letters)}
    terms = [keyfun(i) for i in POS]; labs = sorted({t for tl in terms for t, _ in tl}, key=str); KI = {t: j for j, t in enumerate(labs)}
    nx = len(letters); A = np.zeros((25, nx + len(labs)), np.int64); b = np.zeros(25, np.int64)
    sg = 1 if side == "plain" else -1
    for r, i in enumerate(POS):
        A[r, LI[CR[i] if side == "plain" else ct[i]]] = 1
        for t, co in terms[r]: A[r, nx + KI[t]] += sg * co
        b[r] = X.index(ct[i]) if side == "plain" else X.index(CR[i])
    A[24, 0] = 1                                            # gauge x[letters[0]] = 0
    return A, b, nx, letters

def test(ct, side, X, keyfun, maxnull=6):
    A, b, nx, letters = build(ct, side, X, keyfun); res = {"unknowns": A.shape[1] - 1}
    sols = {}
    for p in (2, 13):
        ok, x, B, rk = solve(A, b, p); res[f"checks{p}"] = 25 - rk
        if not ok: res["cls"] = "inconsistent"; continue
        res[f"null{p}"] = len(B); Bx = proj_basis(B[:, :nx], p); res[f"xnull{p}"] = len(Bx); sols[p] = (x[:nx], Bx)
    if res.get("cls") == "inconsistent": return res
    if max(res["xnull2"], res["xnull13"]) > maxnull: res["cls"] = "skip"; return res
    inj = injective_exists(sols[2][0], sols[2][1], sols[13][0], sols[13][1])
    if not inj: res["cls"] = "noinj"
    else: res["cls"] = "trivial" if res["checks2"] == 0 and res["checks13"] == 0 else "inj"
    return res

def check_support(ct, side, X, kf):
    """crib positions used by the parity checks (left null vectors) at each prime, gauge row excluded."""
    A, b, nx, letters = build(ct, side, X, kf); sup = set()
    for p in (2, 13):
        ok, v0, Bn, rk = solve(A.T.copy(), np.zeros(A.shape[1], np.int64), p)
        for v in Bn: sup.add(tuple(POS[r] for r in range(24) if v[r]))
    return sorted(sup, key=len)

def run(ct, label, out, list_consistent=False):
    CLS = ["inconsistent", "noinj", "inj", "trivial", "skip"]; tab = {}; listed = []
    for fam, name, kf in models():
        for side in ("plain", "cipher"):
            for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
                r = test(ct, side, X, kf); t = tab.setdefault(fam, dict.fromkeys(CLS, 0)); t[r["cls"]] += 1
                if r["cls"] in ("inj", "trivial", "skip"):
                    if r["cls"] == "inj": r["support"] = check_support(ct, side, X, kf)
                    listed.append((fam, name, side, an, r))
    print(f"\n== {label} ==", file=out)
    print(f"{'family':7s} {'models':>6s} " + " ".join(f"{c:>12s}" for c in CLS), file=out)
    for fam, t in tab.items(): print(f"{fam:7s} {sum(t.values()):6d} " + " ".join(f"{t[c]:12d}" for c in CLS), file=out)
    if list_consistent:
        print(f"-- every non-ruled-out model for {label} (cls, checks2/13, nullity2/13, x-nullity2/13, unknowns) --", file=out)
        for fam, name, side, an, r in listed:
            print(f"  {r['cls']:8s} {fam:6s} {name:22s} {side:6s} {an}  checks {r['checks2']}/{r['checks13']}  null {r.get('null2')}/{r.get('null13')}"
                  f"  xnull {r.get('xnull2')}/{r.get('xnull13')}  unk {r['unknowns']}" + (f"  checks-on {r['support']}" if 'support' in r else ""), file=out)
        inj = [x for x in listed if x[4]["cls"] == "inj"]
        only = [x for x in inj if all(set(sp) == {27, 65} for sp in x[4]["support"])]
        print(f"-- {label}: {len(inj)} consistent+injective models with checks; {len(only)} of them are checked ONLY by the pair 27/65 "
              f"(PT R -> CT P twice, i.e. key27 = key65 forced); remaining {len(inj) - len(only)} by family: "
              + str({f: sum(1 for x in inj if x[0] == f and x not in only) for f in tab}), file=out)
    out.flush(); return tab, listed

def plant(side, W=9, seed=5):
    rng = random.Random(seed); P = [rng.choice(k4.AZ) for _ in range(N)]
    for s, w in k4.CRIBS.items(): P[s:s + len(w)] = w
    X = k4.KA; mask = rng.sample(range(26), 26); r = [rng.randrange(26) for _ in range(W)]; c = [rng.randrange(26) for _ in range(-(-N // W))]
    if side == "plain": ct = "".join(X[(mask[X.index(p)] + r[i % W] + c[i // W]) % 26] for i, p in enumerate(P))
    else:  # y[C] = num(P) + z  ->  C = y^-1(num(P)+z)
        inv = {v: k for k, v in enumerate(mask)}; ct = "".join(X[inv[(X.index(p) + r[i % W] + c[i // W]) % 26]] for i, p in enumerate(P))
    return ct

if __name__ == "__main__":
    out = open("mg_results.out", "w")
    class Tee:
        def write(self, s): out.write(s); sys.stdout.write(s)
        def flush(self): out.flush(); sys.stdout.flush()
    T = Tee()
    print("VALIDATION: planted mask (KA numbering) + grid key W=9", file=T)
    for side in ("plain", "cipher"):
        ct = plant(side); res = []
        for W in range(2, 31):
            kf = lambda i, W=W: [(("r", i % W), 1), (("c", i // W), 1)]
            res.append((W, test(ct, side, k4.KA, kf)))
        print(f" planted {side}-side, tested {side}/KA:", file=T)
        for W, r in res: print(f"   W={W:2d} {r['cls']:12s} checks {r['checks2']}/{r['checks13']} xnull {r.get('xnull2')}/{r.get('xnull13')}", file=T)
        oth = [r['cls'] for W, r in res if W != 9]; print("   other W summary:", {k: oth.count(k) for k in set(oth)}, file=T)
        # also the planted text through the whole gridH family (H repeats) for W=9
    T.flush()
    run(k4.CT, "K4", T, list_consistent=True)
    rng = np.random.default_rng(0)
    for k in range(2):
        perm = rng.permutation(N); sh = "".join(k4.CT[j] for j in perm)
        run(sh, f"shuffle {k+1}", T, list_consistent=True)
    # extra baseline: 10 more shuffles, inj counts per family
    agg = {}
    for k in range(10):
        perm = rng.permutation(N); sh = "".join(k4.CT[j] for j in perm)
        tab, _ = run(sh, f"extra shuffle {k+3}", open("/dev/null", "w"))
        for fam, t in tab.items(): agg.setdefault(fam, []).append(t["inj"])
    print("\n== extra baseline: 'inj' counts per family over 10 further shuffles ==", file=T)
    for fam, v in agg.items(): print(f"  {fam:6s} {v}  mean {np.mean(v):.1f}", file=T)
    T.flush()
