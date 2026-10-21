"""Sanborn-phrase alphabets as Chaocipher disks against A-Z, KRYPTOS and the clue alphabets.

Both disk assignments, 26 relative rotations, both directions (ciphertext on the left
disk or the right disk). Scoring kernel is scripts/sm_chao_kernel.c, the same one
used for the clue x dictionary sweep. Phrase x phrase and phrase x dictionary are
estimated from the measured rate and not run.
"""
import os, sys, time, ctypes, subprocess
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import multiprocessing as mp
import phr_common as pc
import k4, fe_common as fc
import sm_taskA as A
from sm_chao import encrypt, encrypt_swapped, decrypt, decrypt_swapped

RES = pc.RES
SCR = pc.SCR
SO = os.path.join(SCR, "sm_chao_kernel.so")
SRC = os.path.join(SCR, "sm_chao_kernel.c")
LOG = os.path.join(RES, "phr_chao.out")
THR = 8
MAXREC = 200000
NPROC = 4
BUDGET_S = 6 * 3600
_lib = None

def build():
    subprocess.check_call(["gcc", "-O3", "-march=native", "-shared", "-fPIC", "-o", SO, SRC])

def lib():
    global _lib
    if _lib is None:
        if not os.path.exists(SO) or os.path.getmtime(SO) < os.path.getmtime(SRC):
            build()
        _lib = ctypes.CDLL(SO)
        P = ctypes.c_void_p
        _lib.sweep.argtypes = [P, ctypes.c_int, P, ctypes.c_int, P, P, ctypes.c_int, P, P, P, ctypes.c_int]
        _lib.sweep.restype = ctypes.c_int
    return _lib

NEED = np.full(97, 255, np.uint8)
NEED[k4.CPOS] = A.PLC

def csweep(clue, dic, ct):
    clue = np.ascontiguousarray(clue, np.uint8)
    dic = np.ascontiguousarray(dic, np.uint8)
    ct = np.ascontiguousarray(ct, np.uint8)
    ht = np.zeros(25, np.int64)
    hb = np.zeros(14, np.int64)
    rec = np.zeros((MAXREC, 7), np.int32)
    nr = lib().sweep(clue.ctypes.data, len(clue), dic.ctypes.data, len(dic),
                     ct.ctypes.data, NEED.ctypes.data, THR, ht.ctypes.data, hb.ctypes.data,
                     rec.ctypes.data, MAXREC)
    return ht, hb, rec[:nr].copy()

# Forked workers read these.
G = {}

def init_worker():
    global _lib
    _lib = None

def job(args):
    """ct, phrase [s,e), clue-index array into the full clue table. Globals are read-only."""
    ct, s, e, cidx = args
    cidx = np.ascontiguousarray(cidx, np.int32)
    ht, hb, rec = csweep(G["C"][cidx], G["P"][s:e], ct)
    if len(rec):
        rec = rec.copy()
        rec[:, 0] = cidx[rec[:, 0]]
        rec[:, 1] += s
    return e - s, ht, hb, rec

def a2s(v):
    return "".join(chr(65 + int(x)) for x in v)

def describe(r):
    a, d, asg, off, sw, tot, b1 = (int(x) for x in r)
    cn, pn = G["cnames"][a], G["pnames"][d]
    L, R = (cn, pn) if asg == 0 else (pn, cn)
    return f"({tot}, {b1}, '{'swap' if sw else 'std'} L={L} R={R} rotR={off}')"

def row_pt(r, ct):
    a, d, asg, off, sw = (int(x) for x in r[:5])
    C, D = a2s(G["C"][a]), a2s(G["P"][d])
    L, R = (C, D) if asg == 0 else (D, C)
    R = R[off:] + R[:off]
    ctxt = "".join(chr(65 + int(v)) for v in ct)
    return (decrypt_swapped if sw else decrypt)(ctxt, L, R)

def summary(tag, ht, hb, recs, extra=""):
    mx = int(np.nonzero(ht)[0].max()) if ht.sum() else -1
    ex = [r for r in recs if int(r[5]) == 24]
    b1 = int(np.nonzero(hb)[0].max()) if hb.sum() else -1
    line = (f"[{tag}] rows {int(ht.sum()):,}  exact(24/24) {len(ex)}  best partial {mx}/24 "
            f"(#rows at best {int(ht[mx]) if mx >= 0 else 0})  best EASTNORTHEAST-block {b1}/13  "
            f"hist>=8: {{{', '.join(f'{k}:{int(ht[k])}' for k in range(8, 25) if ht[k])}}}{extra}")
    return line, ex

def numpy_hist(clue, dic, ct):
    """Independent histogram via sm_taskA.run, sw innermost, matching the C loop order's multiset."""
    na, nd = len(clue), len(dic)
    Ls, Rs = [], []
    for d in range(nd):
        for a in range(na):
            for asg in range(2):
                for off in range(26):
                    L = dic[d] if asg else clue[a]
                    R = clue[a] if asg else dic[d]
                    Ls.append(L)
                    Rs.append(R[A.ROT[off]])
    Ls = np.array(Ls, np.uint8)
    Rs = np.array(Rs, np.uint8)
    ht = np.zeros(25, np.int64)
    hb = np.zeros(14, np.int64)
    for sw in (False, True):
        tot, b1 = A.run(ct, Ls, Rs, sw)
        ht += np.bincount(tot, minlength=25)
        hb += np.bincount(b1, minlength=14)
    return ht, hb

def sweep_all(ct, spans, cidx, pool, log, tag):
    t1 = time.time()
    HT = np.zeros(25, np.int64)
    HB = np.zeros(14, np.int64)
    RECS = []
    done = 0
    nd = sum(e - s for s, e in spans)
    nextp = 0.05
    tasks = [(ct, s, e, cidx) for s, e in spans]
    for m, ht, hb, rec in pool.imap_unordered(job, tasks, chunksize=1):
        HT += ht
        HB += hb
        if len(rec):
            RECS.extend(rec.tolist())
        done += m
        frac = done / nd
        if frac + 1e-12 >= nextp or done == nd:
            el = time.time() - t1
            line, _ = summary(f"{tag} {100 * frac:5.1f}%", HT, HB, RECS,
                              f"  elapsed {el/60:.1f} min ETA {el * (1 - frac) / max(frac, 1e-9) / 60:.1f} min")
            log(line)
            while nextp <= frac + 1e-12:
                nextp += 0.05
    return HT, HB, RECS, time.time() - t1

def main():
    build()
    f = open(LOG, "w")
    def log(s):
        print(s, flush=True)
        print(s, file=f, flush=True)
    P, pnames = pc.phrase_alphabets()
    cl = fc.clue_alphabets()
    C = np.array([[ord(c) - 65 for c in a] for a, _ in cl], np.uint8)
    cnames = [n for _, n in cl]
    log(f"phrase alphabets {len(pnames):,} sha256/16={pc.digest(P)}")
    log(f"clue alphabets {len(cnames)}: {cnames[:8]} ...")
    fixed = [n for n in ("AZ", "KA", "AZrev", "KArev") if n in cnames]
    log(f"A-Z / KA disks inside the clue set: {fixed}")
    # Rubin test vector, so the Python reference the kernel is checked against is intact.
    L0, R0 = "HXUCZVAMDSLKPEFJRIGTWOBNYQ", "PTLNBQDEOYSFAVZKGJRIHWXUMC"
    tv = encrypt("WELLDONEISBETTERTHANWELLSAID", L0, R0)
    log(f"Rubin test vector: {tv == 'OAHQHCNYNXTSZJRRHJBYHQKSOUJY'}")
    if tv != "OAHQHCNYNXTSZJRRHJBYHQKSOUJY":
        log("ABORT: Rubin test vector failed")
        sys.exit(1)
    G["P"], G["pnames"], G["C"], G["cnames"] = P, pnames, C, cnames
    CT4 = A.codes(k4.CT)
    # ---- C vs numpy on 3 phrase x 4 clue ----
    ht_c, hb_c, _ = csweep(C[:4], P[:3], CT4)
    ht_n, hb_n = numpy_hist(C[:4], P[:3], CT4)
    same = np.array_equal(ht_c, ht_n) and np.array_equal(hb_c, hb_n)
    log(f"cross-check C vs numpy on {3*4*2*26*2} rows: histograms identical = {same}")
    if not same:
        log(f"C { {k: int(ht_c[k]) for k in range(25) if ht_c[k]} }")
        log(f"N { {k: int(ht_n[k]) for k in range(25) if ht_n[k]} }")
        log("ABORT: kernel histogram mismatch")
        sys.exit(1)
    # ---- planted: phrase disk x clue disk, both role directions ----
    prng = np.random.default_rng(99)
    pt = fc.planted_plaintext(prng)
    def find_name(pred):
        for i, n in enumerate(pnames):
            if pred(n):
                return i
        return None
    pi = find_name(lambda n: n == "BERLINCLOCK/s")
    if pi is None:
        pi = find_name(lambda n: n.startswith("BERLIN"))
    ci_az = cnames.index("AZ")
    ci_kw = cnames.index("KA")
    plants = [
        ("plant-std", pi, ci_kw, 1, 9, 0, encrypt),          # asg 1: L=phrase R=KA rot 9, std
        ("plant-swap", pi, ci_az, 0, 4, 1, encrypt_swapped),  # asg 0: L=AZ R=phrase rot 4, swap
    ]
    plant_ok = True
    for nm, pii, cii, asg, off, sw, enc in plants:
        phrase, clue = a2s(P[pii]), a2s(C[cii])
        L = phrase if asg == 1 else clue
        R0s = clue if asg == 1 else phrase
        R = R0s[off:] + R0s[:off]
        pct = enc(pt, L, R)
        # window of phrase alphabets around the planted one, full clue set
        s, e = max(0, pii - 2), min(len(P), pii + 3)
        ht, hb, rec = csweep(C, P[s:e], A.codes(pct))
        if len(rec):
            rec[:, 1] += s
        line, ex = summary(f"planted {nm}", ht, hb, rec.tolist() if len(rec) else [])
        log(line)
        recovered = any(row_pt(r, A.codes(pct)) == pt for r in ex)
        want = [r for r in ex if int(r[0]) == cii and int(r[1]) == pii and int(r[2]) == asg and int(r[3]) == off and int(r[4]) == sw]
        log(f"    planted L={L[:12]}... Rrot={off} asg={asg} sw={sw} alph={pnames[pii]} clue={cnames[cii]} "
            f"exact_rows={len(ex)} planted_row={bool(want)} plaintext_recovered={recovered}")
        plant_ok &= recovered and bool(want)
    if not plant_ok:
        log("ABORT: planted Chaocipher message was not recovered")
        sys.exit(1)
    if os.environ.get("PHR_VALIDATE_ONLY"):
        log("PHR_VALIDATE_ONLY: stopping after planted recovery")
        return
    # ---- rate and scope ----
    t0 = time.time()
    csweep(C, P[:150], CT4)
    dt = time.time() - t0
    rate = 150 * len(C) * 104 / dt  # rows/s, one process
    rows_ct = len(P) * len(C) * 104
    est = rows_ct * 3 / rate / NPROC
    log(f"kernel {rate/1e6:.3f} M rows/s/process. Full phrase x {len(C)} clue x 2 assign x 26 off x 2 dir "
        f"= {rows_ct:,} rows/ciphertext, x3 = {rows_ct*3:,}. Estimated {est/3600:.2f} h on {NPROC} processes.")
    # phrase x phrase (ordered, one side rotated, 2 directions; assignments would double-count)
    rows_pp = len(P) * (len(P) - 1) * 26 * 2
    est_pp = rows_pp * 3 / rate / NPROC / 3600
    # phrase x the Session 2 dictionary+clue set (842,054), both assignments
    ndict = 842054
    rows_pd = len(P) * ndict * 104
    est_pd = rows_pd * 3 / rate / NPROC / 3600
    log(f"NOT RUN phrase x phrase (ordered pairs, 26 rot, 2 dir, 3 ct): {rows_pp*3:,.3e} rows, estimate {est_pp:.0f} h.")
    log(f"NOT RUN phrase x dictionary+clue ({ndict:,}, both disks, 3 ct): {rows_pd*3:,.3e} rows, estimate {est_pd:.0f} h.")
    # Passes are (phrase spans, clue-index array). Full run is one pass.
    passes = []
    if est > BUDGET_S:
        chunk = 1000
        nfixed = len(fixed)
        t_fixed = len(P) * nfixed * 104 * 3 / rate / NPROC
        remain = max(0.0, BUDGET_S - t_fixed)
        m_raw = int(remain * rate * NPROC / (max(len(C) - nfixed, 1) * 104 * 3))
        M = max(chunk, (m_raw // chunk) * chunk)
        M = min(M, len(P))
        fixed_idx = np.array([cnames.index(n) for n in fixed], np.int32)
        rest_idx = np.array([i for i, n in enumerate(cnames) if n not in fixed], np.int32)
        spans_all = [(s, min(s + chunk, len(P))) for s in range(0, len(P), chunk)]
        spans_pre = [(s, min(s + chunk, M)) for s in range(0, M, chunk)]
        passes.append((f"fixed {fixed}", spans_all, fixed_idx))
        if len(rest_idx) and M:
            passes.append((f"prefix {M} x {len(rest_idx)} clue", spans_pre, rest_idx))
        cov = (len(P) * nfixed + M * len(rest_idx)) / (len(P) * len(C))
        log(f"SUBSET: full estimate {est/3600:.2f} h > 6 h. "
            f"Pass A: all {len(P):,} phrase alphabets x {fixed}. "
            f"Pass B: the first {M:,} phrase alphabets (gap_phr order) x the other {len(rest_idx)} clue alphabets. "
            f"Same passes for K4 and both shuffles. Coverage {cov:.4f} of phrase x clue.")
    else:
        cov = 1.0
        M = len(P)
        passes.append(("full", [(s, min(s + 1000, len(P))) for s in range(0, len(P), 1000)],
                        np.arange(len(C), dtype=np.int32)))
        log("Running the full phrase x clue set on K4 and both shuffles.")
    rng = np.random.default_rng(0)
    sh = ["".join(rng.permutation(list(k4.CT))) for _ in range(2)]
    log(f"shuf0 {sh[0]}")
    log(f"shuf1 {sh[1]}")
    cts = [("K4", CT4), ("shuf0", A.codes(sh[0])), ("shuf1", A.codes(sh[1]))]
    wall0 = time.time()
    results = {}
    try:
        mp.set_start_method("fork")
    except RuntimeError:
        pass
    with mp.Pool(NPROC, initializer=init_worker) as pool:
        for nm, ct in cts:
            HT = np.zeros(25, np.int64)
            HB = np.zeros(14, np.int64)
            RECS = []
            for tag, spans, cidx in passes:
                ht, hb, recs, dt = sweep_all(ct, spans, cidx, pool, log, f"{nm} {tag}")
                HT += ht
                HB += hb
                RECS.extend(recs)
                log(f"[{nm} {tag} done] {dt/60:.1f} min")
            line, ex = summary(f"FINAL {nm}", HT, HB, RECS)
            log(line)
            RECS.sort(key=lambda r: (-int(r[5]), -int(r[6])))
            for r in RECS[:5]:
                log(f"    top: {describe(r)}")
            for r in ex[:10]:
                log(f"    EXACT: {describe(r)} pt={row_pt(r, ct)}")
            if len(ex) != int(HT[24]):
                log(f"NOTE: recorded exact rows {len(ex)} != histogram {int(HT[24])} (record cap {MAXREC})")
            results[nm] = (int(np.nonzero(HT)[0].max()), int(HT[24]), int(HT.sum()))
    elapsed = time.time() - wall0
    rows = results["K4"][2] + results["shuf0"][2] + results["shuf1"][2]
    log(f"SUMMARY family=chaocipher rows_scored={rows:,} "
        f"exact_k4={results['K4'][1]} exact_s0={results['shuf0'][1]} exact_s1={results['shuf1'][1]} "
        f"partial_k4={results['K4'][0]}/24 partial_s0={results['shuf0'][0]}/24 partial_s1={results['shuf1'][0]}/24 "
        f"planted=pass runtime_s={elapsed:.0f} coverage={cov:.4f} "
        f"phrase_x_phrase_est_h={est_pp:.0f} phrase_x_dict_est_h={est_pd:.0f}")
    log(f"done {time.ctime()}")

if __name__ == "__main__":
    main()
