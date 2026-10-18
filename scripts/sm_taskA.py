"""Task A: Chaocipher with keyword alphabets on K4 (+ 2 shuffled controls + planted validation)."""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, multiprocessing as mp
import k4, fe_common as fc
from sm_chao import SIG_L, SIG_R, encrypt, encrypt_swapped, decrypt, decrypt_swapped

TL = ((np.arange(26)[:, None] + SIG_L[None, :]) % 26).astype(np.int32)
TR = ((np.arange(26)[:, None] + SIG_R[None, :]) % 26).astype(np.int32)
ROT = ((np.arange(26)[:, None] + np.arange(26)[None, :]) % 26)          # ROT[d] : rot by d
POS = k4.CPOS; PLC = np.array([ord(c) - 65 for c in k4.CPL], np.uint8)
NST = 74
def codes(s): return np.array([ord(c) - 65 for c in s], np.uint8)

def run(ct, L, R, swapped):
    n = L.shape[0]; base = (np.arange(n, dtype=np.int32) * 26)[:, None]
    out = np.empty((n, NST), np.uint8); Lf = np.ascontiguousarray(L).ravel(); Rf = np.ascontiguousarray(R).ravel()
    for t in range(NST):
        A = (Rf if swapped else Lf).reshape(n, 26)
        i = np.argmax(A == ct[t], axis=1)
        out[:, t] = (Lf if swapped else Rf)[base[:, 0] + i]
        Lf = Lf[base + TL[i]].ravel(); Rf = Rf[base + TR[i]].ravel()
    m = out[:, POS] == PLC[None, :]
    return m.sum(1), m[:, :13].sum(1)

G = {}
def init():
    al = fc.all_alphabets()
    G["names"] = [n for a, n in al]
    G["D"] = np.array([[ord(c) - 65 for c in a] for a, n in al], np.uint8)
    cl = fc.clue_alphabets()
    G["cnames"] = [n for a, n in cl]; G["C"] = np.array([[ord(c) - 65 for c in a] for a, n in cl], np.uint8)
    G["F"] = np.array([codes(x) for x in (k4.AZ, k4.KA, k4.AZ[::-1], k4.KA[::-1])])
    G["Fn"] = ["AZ", "KA", "AZrev", "KArev"]
    rng = np.random.default_rng(5); nd = len(al)
    G["pairs"] = rng.integers(0, nd, size=(1_000_000, 2))

def rows_for(kind, s, e):
    """returns L, R arrays and a descriptor function(row)->str"""
    d26 = np.arange(26)
    if kind == "clue":                     # clue x clue, pairs indices s..e of the flattened 116x116
        C = G["C"]; nc = len(C); pi = np.arange(s, e); a, b = pi // nc, pi % nc
        L = np.repeat(C[a], 26, 0); R = C[b][:, None, :][:, :, None, :].squeeze(2)
        R = C[b][:, ROT].reshape(-1, 26)    # (k,26 rots,26)
        def desc(r): k, d = divmod(r, 26); return f"L={G['cnames'][a[k]]} R={G['cnames'][b[k]]} rotR={d}"
        return L, R, desc
    if kind == "dictF":                    # dict alphabet idx s..e  x 4 fixed, both sides
        D = G["D"][s:e]; F = G["F"]; m = len(D)
        La = np.repeat(np.repeat(D, 4, 0), 26, 0)                                   # (m*4*26)
        Ra = F[:, ROT].reshape(4 * 26, 26)[None].repeat(m, 0).reshape(-1, 26)
        Lb = np.repeat(np.tile(F, (m, 1)), 26, 0)
        Rb = np.repeat(D[:, ROT], 4, 0).reshape(-1, 26)
        Rb = D[:, ROT][:, None].repeat(4, 1).reshape(-1, 26)
        L = np.concatenate([La, Lb]); R = np.concatenate([Ra, Rb]); half = m * 104
        def desc(r):
            side = "Ldict" if r < half else "Rdict"; r %= half
            k, rest = divmod(r, 104); f, d = divmod(rest, 26)
            nm = G["names"][s + k]; fn = G["Fn"][f]
            return f"L={nm} R={fn} rotR={d}" if side == "Ldict" else f"L={fn} R={nm} rotR={d}"
        return L, R, desc
    if kind == "dd":                       # dict x dict random pairs s..e
        P = G["pairs"][s:e]; D = G["D"]
        L = np.repeat(D[P[:, 0]], 26, 0); R = D[P[:, 1]][:, ROT].reshape(-1, 26)
        def desc(r): k, d = divmod(r, 26); return f"L={G['names'][P[k,0]]} R={G['names'][P[k,1]]} rotR={d}"
        return L, R, desc

THR = 8
def job(args):
    ctname, ct, kind, s, e = args
    L, R, desc = rows_for(kind, s, e)
    res = {"ht": np.zeros(25, np.int64), "hb": np.zeros(14, np.int64), "top": [], "exact": [], "n": 0}
    for sw in (False, True):
        tot, b1 = run(ct, L, R, sw)
        res["ht"] += np.bincount(tot, minlength=25); res["hb"] += np.bincount(b1, minlength=14); res["n"] += len(tot)
        for r in np.nonzero(tot >= THR)[0]:
            res["top"].append((int(tot[r]), int(b1[r]), ("swap " if sw else "std ") + desc(int(r))))
        for r in np.nonzero(tot == 24)[0]:
            l = "".join(chr(65 + v) for v in L[r]); rr = "".join(chr(65 + v) for v in R[r])
            ctxt = "".join(chr(65 + v) for v in ct)
            pt = (decrypt_swapped if sw else decrypt)(ctxt, l, rr)
            res["exact"].append((("swap " if sw else "std ") + desc(int(r)), pt))
    return res

def merge(rs):
    out = {"ht": sum(r["ht"] for r in rs), "hb": sum(r["hb"] for r in rs), "n": sum(r["n"] for r in rs),
           "top": sorted([t for r in rs for t in r["top"]], reverse=True), "exact": [x for r in rs for x in r["exact"]]}
    return out

def report(tag, ctname, m, f):
    ht = m["ht"]; mx = int(np.nonzero(ht)[0].max())
    line = (f"[{tag}] {ctname}: rows {m['n']:,}  exact(24/24) {len(m['exact'])}  best partial {mx}/24 "
            f"(#rows at best {int(ht[mx])})  best EASTNORTHEAST-block {int(np.nonzero(m['hb'])[0].max())}/13  "
            f"hist>=7: {{{', '.join(f'{k}:{int(ht[k])}' for k in range(7, 25) if ht[k])}}}")
    print(line, file=f, flush=True); print(line, flush=True)
    for t in m["top"][:3]: print(f"      top: {t}", file=f, flush=True)
    for x in m["exact"][:20]: print(f"      EXACT: {x}", file=f, flush=True)

if __name__ == "__main__":
    init()
    nc = len(G["C"]); nd = len(G["D"])
    print(f"clue alphabets {nc}, dictionary+clue alphabets {nd:,}", flush=True)
    rng = np.random.default_rng(0)
    sh = ["".join(rng.permutation(list(k4.CT))) for _ in range(2)]
    f = open("sm_taskA.out", "a")
    print(f"=== Task A run {time.ctime()} ===", file=f)
    # ---- planted validation ----
    prng = np.random.default_rng(42)
    planted = []
    pt = fc.planted_plaintext(prng)
    di = G["names"].index("SCULPTURE/t") if "SCULPTURE/t" in G["names"] else 123456
    lD = "".join(chr(65 + v) for v in G["D"][di]); rKA = k4.KA[9:] + k4.KA[:9]
    ct1 = encrypt(pt, lD, rKA)                                   # std, L=dict, R=KA rot 9
    ct2 = encrypt_swapped(pt, fc.keyed("BERLIN"), fc.keyed("PALIMPSEST")[::-1][5:] + fc.keyed("PALIMPSEST")[::-1][:5])
    print("planted pt:", pt, file=f)
    vjobs = [("plant1", codes(ct1), "dictF", di - 3, di + 5), ("plant2", codes(ct2), "clue", 0, nc * nc)]
    with mp.Pool(3, initializer=init) as pool:
        for (nm, ct, kind, s, e) in vjobs:
            if kind == "clue":
                rs = pool.map(job, [(nm, ct, kind, a, min(a + 1000, e)) for a in range(s, e, 1000)])
            else:
                rs = [job((nm, ct, kind, s, e))]
            m = merge(rs); report("planted " + kind, nm, m, f)
            print(f"      planted recovered: {any(p == pt for _, p in m['exact'])}", file=f, flush=True)
        # ---- real runs ----
        cts = [("K4", codes(k4.CT)), ("shuf0", codes(sh[0])), ("shuf1", codes(sh[1]))]
        plan = [("clue", 0, nc * nc, 1000), ("dictF", 0, nd, 400), ("dd", 0, len(G["pairs"]), 4000)]
        for kind, s, e, step in plan:
            for nm, ct in cts:
                t0 = time.time()
                rs = pool.map(job, [(nm, ct, kind, a, min(a + step, e)) for a in range(s, e, step)], chunksize=4)
                report(kind, nm, merge(rs), f); print(f"      ({time.time()-t0:.0f}s)", file=f, flush=True)
