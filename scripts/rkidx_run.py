"""Exact-lookup running-key test for K4.
Index every 13-gram / 11-gram of a text corpus (fwd+rev) as exact base-26 int64 codes.
For every alphabet/arrangement/convention/indicator compute the 24 implied key letters at the crib
positions; look up block1 (pos 21..33) and block2 (63..73); a full running-key hit = both blocks occur in
the same text with block2 starting exactly 42 letters after block1 (verified by exact scan).
Usage: python3 rkidx_run.py validate | python3 rkidx_run.py full
"""
import sys, glob, re, time, numpy as np, multiprocessing as mp
import k4, texts

AZ, KA = k4.AZ, k4.KA
CPOS = np.array(k4.CPOS); CPL = k4.CPL
P13 = 26 ** np.arange(12, -1, -1, dtype=np.int64)
P11 = 26 ** np.arange(10, -1, -1, dtype=np.int64)
A2I = lambda s: np.frombuffer(s.encode(), dtype=np.uint8).astype(np.int64) - 65

# ---------------- corpus ----------------
def strip(s): return re.sub("[^A-Z]", "", s.upper())
def strip_de(s):
    s = s.upper()
    for x, y in (("Ä", "AE"), ("Ö", "OE"), ("Ü", "UE"), ("ß", "SS"), ("ẞ", "SS")): s = s.replace(x, y)
    return re.sub("[^A-Z]", "", s)
def rd(f): return open(f, encoding="utf-8", errors="ignore").read()
def corpus():
    T = []
    for f in sorted(glob.glob("g*.txt")): T.append((f, strip(rd(f))))
    for f in sorted(glob.glob("de*.txt")): T.append((f, strip_de(rd(f))))
    T.append(("sanborn_oral_history.txt", strip(rd("sanborn_oral_history.txt"))))
    T.append(("K1+K2end+K3", texts.K1 + texts.K2_ORIG_END + texts.K3))
    T.append(("panel.txt", strip(rd("panel.txt"))))
    T.append(("tableau.txt", strip(rd("tableau.txt"))))
    try: T.append(("rkidx_carter.txt", strip(rd("rkidx_carter.txt"))))
    except FileNotFoundError: pass
    out = []
    for n, t in T:
        out.append((n + ":fwd", t)); out.append((n + ":rev", t[::-1]))
    return out

def grams(a, L, P):
    n = len(a) - L + 1
    if n <= 0: return np.zeros(0, np.int64)
    h = np.zeros(n, np.int64)
    for j in range(L): h = h * 26 + a[j:j + n]
    return h

TXT = None; H13 = None; H11 = None; U13 = None; U11 = None
def build_index():
    global TXT, H13, H11, U13, U11
    TXT = corpus(); H13 = []; H11 = []
    for n, t in TXT:
        a = A2I(t); H13.append(grams(a, 13, P13)); H11.append(grams(a, 11, P11))
    U13 = np.unique(np.concatenate(H13)); U11 = np.unique(np.concatenate(H11))

def member(U, q):
    i = np.searchsorted(U, q); i[i == len(U)] = 0
    return U[i] == q

def full_locate(q1, q2):
    """exact: offsets o with T[o+21..o+33] code q1 and T[o+63..o+73] code q2."""
    res = []
    for ti, (n, t) in enumerate(TXT):
        h13, h11 = H13[ti], H11[ti]
        cand = np.nonzero(h13 == q1)[0]            # cand = o+21
        for p in cand:
            o = p - 21
            if o >= 0 and o + 63 < len(h11) and h11[o + 63] == q2: res.append((ti, int(o)))
    return res

# ---------------- configs ----------------
AZA = np.arange(26); KAA = A2I(KA); KAINV = np.argsort(KAA)
def dec_letter(code, L):
    out = []
    for _ in range(L): out.append(chr(65 + code % 26)); code //= 26
    return "".join(out[::-1])

def run_chunk(args):
    """alphabets: list of (alph, name); cts: list of (label, ct). mode 'std' or 'mix'."""
    alphs, cts, mode = args
    stats = {lab: dict(cfg=0, b1=0, b2=0, both=0, full=0) for lab, _ in cts}
    hits = {lab: [] for lab, _ in cts}         # block-level hit records (capped)
    fulls = {lab: [] for lab, _ in cts}
    PLi = A2I(CPL)
    if mode == "std":
        A = np.stack([AZA, KAA]); names = ["AZ", "KA"]
        arrs = [("AZ", "AZ"), ("KA", "KA"), ("AZ", "KA"), ("KA", "AZ")]
        inds = ["AZ", "KA"]
    else:
        A = np.array([A2I(a) for a, _ in alphs]); names = [nm for _, nm in alphs]
        arrs = [("AZ", "A"), ("A", "AZ"), ("KA", "A"), ("A", "KA"), ("A", "A")]
        inds = ["AZ", "KA", "A"]
    n = len(A); INV = np.argsort(A, axis=1)
    def inv(x):  # letter -> slot, (n,26)
        return {"AZ": np.broadcast_to(AZA, (n, 26)), "KA": np.broadcast_to(KAINV, (n, 26)), "A": INV}[x]
    for lab, ct in cts:
        Ci = A2I(ct)[CPOS]
        for pa, ca in arrs:
            if mode == "std":   # only one row needed for fixed pair; use row 0
                pass
            pc = inv(pa)[:, PLi]; cc = inv(ca)[:, Ci]
            for conv in ("vig", "beau", "vb"):
                k = {"vig": cc - pc, "beau": cc + pc, "vb": pc - cc}[conv] % 26
                for ind in inds:
                    if ind == "AZ": K = k
                    elif ind == "KA": K = KAA[k]
                    else: K = np.take_along_axis(A, k, 1)
                    if mode == "std": K = K[:1]; rows = 1
                    else: rows = n
                    q1 = K[:, :13] @ P13; q2 = K[:, 13:] @ P11
                    m1 = member(U13, q1); m2 = member(U11, q2); mb = m1 & m2
                    st = stats[lab]; st["cfg"] += rows; st["b1"] += int(m1.sum()); st["b2"] += int(m2.sum()); st["both"] += int(mb.sum())
                    for j in np.nonzero(m1 | m2)[0]:
                        if len(hits[lab]) < 2000:
                            hits[lab].append((names[j] if mode == "mix" else "std", pa, ca, conv, ind, bool(m1[j]), bool(m2[j]),
                                              "".join(chr(65 + x) for x in K[j])))
                    for j in np.nonzero(mb)[0]:
                        for ti, o in full_locate(int(q1[j]), int(q2[j])):
                            fulls[lab].append((names[j] if mode == "mix" else "std", alphs[j][0] if mode == "mix" else None,
                                               pa, ca, conv, ind, ti, o))
                            stats[lab]["full"] += 1
    return stats, hits, fulls

def decrypt(ct, alph, pa, ca, conv, ind, ti, o):
    t = TXT[ti][1]
    al = {"AZ": AZ, "KA": KA, "A": alph}
    PA, CA, I = al[pa], al[ca], al[ind]
    out = []
    for i, c in enumerate(ct):
        if o + i >= len(t): out.append("?"); continue
        kv = I.index(t[o + i]); cv = CA.index(c)
        pv = {"vig": cv - kv, "beau": kv - cv, "vb": cv + kv}[conv] % 26
        out.append(PA[pv])
    return "".join(out)

def run_all(cts, alph_groups, out, workers=4, chunk=20000):
    tot = {lab: dict(cfg=0, b1=0, b2=0, both=0, full=0) for lab, _ in cts}
    H = {lab: [] for lab, _ in cts}; F = {lab: [] for lab, _ in cts}
    jobs = [([], cts, "std")]
    for g in alph_groups:
        jobs += [(g[i:i + chunk], cts, "mix") for i in range(0, len(g), chunk)]
    ctx = mp.get_context("fork")
    t0 = time.time()
    with ctx.Pool(workers) as pool:
        for k, (st, hi, fu) in enumerate(pool.imap_unordered(run_chunk, jobs)):
            for lab in st:
                for x in st[lab]: tot[lab][x] += st[lab][x]
                H[lab] += hi[lab]; F[lab] += fu[lab]
            if k % 10 == 0: print(f"  chunk {k+1}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    return tot, H, F

def log(out, s):
    print(s, flush=True); out.write(s + "\n"); out.flush()

if __name__ == "__main__":
    mode = sys.argv[1]
    t0 = time.time(); build_index()
    out = open("rkidx_results.out" if mode == "full" else f"rkidx_{mode}.out", "w")
    log(out, f"corpus: {len(TXT)} texts (fwd+rev); letters per direction:")
    for n, t in TXT[::2]: log(out, f"   {n[:-4]:28s} {len(t):>10,}")
    tot_letters = sum(len(t) for _, t in TXT)
    log(out, f"total letters (fwd+rev) {tot_letters:,}; full 24-tuple offsets {sum(max(0,len(t)-73) for _,t in TXT):,}; "
             f"distinct 13-grams {len(U13):,}; distinct 11-grams {len(U11):,}  [{time.time()-t0:.0f}s]")
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        import round3
    ALPHA = [(a, "r3:" + nm) for nm, a in round3.ALPHA.items()]
    rng = np.random.default_rng(0)
    if mode == "validate":
        import rk_mixed
        # plant 1: KA tableau Vigenere (K1/K2 style), key text g98 at offset 123456
        g98 = [t for n, t in TXT if n == "g98.txt:fwd"][0]; o1 = 123456
        r = np.random.default_rng(42)
        def plain():
            p = list("".join(chr(65 + x) for x in r.integers(0, 26, 97)))
            for s, w in k4.CRIBS.items(): p[s:s + len(w)] = w
            return "".join(p)
        P1 = plain()
        C1 = "".join(KA[(KA.index(p) + KA.index(g98[o1 + i])) % 26] for i, p in enumerate(P1))
        # plant 2: mixed alphabet from rk_mixed, PA=A, CA=KA, beaufort, indicator A, key text g98 at offset 345678
        AL = rk_mixed.alph_list(); jA = 500000; alph, an = AL[jA]; o2 = 345678
        P2 = plain()
        C2 = "".join(KA[(alph.index(g98[o2 + i]) - alph.index(p)) % 26] for i, p in enumerate(P2))
        log(out, f"plant1 KA/KA vig ind KA, g98 off {o1}: CT {C1}")
        log(out, f"plant2 {an} PA=A CA=KA beau ind A, g98 off {o2}: CT {C2}")
        cts = [("plant1", C1), ("plant2", C2)]
        groups = [ALPHA, AL[jA - 30000:jA + 30000]]
    else:
        sys.argv = sys.argv[:1]; import rk_mixed, gap2
        AL = rk_mixed.alph_list()
        with contextlib.redirect_stdout(io.StringIO()):
            G2, nde = gap2.build()
        log(out, f"alphabets: round3.ALPHA {len(ALPHA)}, rk_mixed {len(AL):,}, gap2 {len(G2):,}")
        sh = ["".join(rng.permutation(list(k4.CT))) for _ in range(2)]
        cts = [("K4", k4.CT), ("shuf1", sh[0]), ("shuf2", sh[1])]
        for lab, c in cts: log(out, f"{lab}: {c}")
        groups = [ALPHA, AL, G2]
    tot, H, F = run_all(cts, groups, out)
    ctd = dict(cts)
    for lab, _ in cts:
        s = tot[lab]
        log(out, f"\n== {lab}: configs {s['cfg']:,}  block1(13) hits {s['b1']}  block2(11) hits {s['b2']}  "
                 f"both-blocks {s['both']}  EXACT 24-crib running-key hits {s['full']}")
        uniq = sorted({h[-1] for h in H[lab]})
        log(out, f"   distinct key tuples with a block hit: {len(uniq)}")
        for h in H[lab][:60]: log(out, f"   blockhit {h}")
        for (nm, alph, pa, ca, conv, ind, ti, o) in F[lab]:
            d = decrypt(ctd[lab], alph, pa, ca, conv, ind, ti, o)
            ok = all(d[s + i] == w[i] for s, w in k4.CRIBS.items() for i in range(len(w)))
            log(out, f"   FULL HIT alph={nm} ({alph}) PA={pa} CA={ca} {conv} ind={ind} text={TXT[ti][0]} offset={o} cribs_ok={ok}")
            log(out, f"      key  {TXT[ti][1][o:o+97]}")
            log(out, f"      dec  {d}")
    log(out, f"\ndone in {time.time()-t0:.0f}s")
