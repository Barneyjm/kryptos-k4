"""Sanborn-phrase alphabets x the 1.35M clock-derived additive key streams.

Exact test reuses sm_taskB: Q3 / Q1(AZ,KA) / Q2(AZ,KA) x vig/beau, sigma in {+1,-1},
offset b in Z26, by hashing normalized implied keys. Partial score reuses the same
definition as sm_taskB.best_partial (mode of key - sigma*stream over the 24 cribs),
extended from the 116 clue alphabets to every phrase alphabet. Quag IV clue x clue
is not repeated.
"""
import os, sys, time, ctypes, subprocess
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import phr_common as pc
import k4, fe_common as fc
import sm_taskB as B

RES = pc.RES
SO = os.path.join(RES, "phr_clock_kernel.so")
SRC = os.path.join(RES, "phr_clock_kernel.c")
LOG = os.path.join(RES, "phr_clock.out")
BUDGET_S = 6 * 3600

def build():
    subprocess.check_call(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC", "-o", SO, SRC])

def load_lib():
    if not os.path.exists(SO) or os.path.getmtime(SO) < os.path.getmtime(SRC):
        build()
    lib = ctypes.CDLL(SO)
    lib.clock_partial.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                                  ctypes.c_void_p, ctypes.c_void_p]
    lib.clock_partial.restype = None
    return lib

def partial(lib, keys, pat):
    keys = np.ascontiguousarray(keys, np.uint8)
    pat = np.ascontiguousarray(pat, np.uint8)
    hist = np.zeros(25, np.uint64)
    best = np.full(5, -1, np.int32)
    lib.clock_partial(keys.ctypes.data, len(keys), pat.ctypes.data, len(pat), hist.ctypes.data, best.ctypes.data)
    return hist, best

def py_mode(k, p, sigma):
    d = (k.astype(int) - sigma * p.astype(int)) % 26
    return int(np.bincount(d, minlength=26).max()), int(np.bincount(d, minlength=26).argmax())

def self_check(lib, pat):
    rng = np.random.default_rng(2)
    keys = rng.integers(0, 26, size=(6, 24), dtype=np.uint8)
    sub = pat[:80]
    # known scores
    keys[0] = sub[0]
    keys[1] = sub[0].copy(); keys[1, 3] = (int(keys[1, 3]) + 3) % 26
    hist, best = partial(lib, keys, sub)
    ok = True
    for i in range(6):
        for s in range(80):
            for sg, code in ((1, 0), (-1, 1)):
                m, b = py_mode(keys[i], sub[s], sg)
                # recount via a one-row call would be slow; check the forced rows and a sample
    m0, b0 = py_mode(keys[0], sub[0], 1)
    m1, b1 = py_mode(keys[1], sub[0], 1)
    # full python best over the 6 x 80 x 2
    py_best = 0
    py_hist = np.zeros(25, np.int64)
    for i in range(6):
        for s in range(80):
            for sg in (1, -1):
                m, _ = py_mode(keys[i], sub[s], sg)
                py_hist[m] += 1
                if m > py_best:
                    py_best = m
    ok = (m0 == 24) and (m1 == 23) and int(best[0]) == py_best and np.array_equal(hist, py_hist)
    # also match sm_taskB.best_partial on these keys against `sub` as the pattern matrix
    bp, arg = B.best_partial(keys.astype(np.int16), sub.astype(np.int16))
    ok = ok and bp == py_best
    return ok, f"forced 24/23={m0}/{m1} Cbest={int(best[0])} pybest={py_best} hist_match={np.array_equal(hist, py_hist)} taskB={bp}"

def main():
    build()
    lib = load_lib()
    f = open(LOG, "w")
    def log(s):
        print(s, flush=True)
        print(s, file=f, flush=True)
    A, names = pc.phrase_alphabets()
    log(f"phrase alphabets {len(names):,} sha256/16={pc.digest(A)}")
    t0 = time.time()
    S = B.build_streams(); snames = list(S)
    POS = B.POS
    # int16: sigma = -1 must not wrap the way uint8 arithmetic would.
    SC = np.array([S[n][POS] for n in snames], np.int16)
    SCU, uinv = np.unique(SC.astype(np.uint8), axis=0, return_inverse=True)
    urep = {}
    for i, u in enumerate(uinv):
        urep.setdefault(int(u), i)
    log(f"streams {len(snames):,} distinct crib patterns {len(SCU):,} ({time.time()-t0:.0f}s)")
    table = {}
    for sg in B.SIGMAS:
        signed = (sg * SC.astype(np.int16)) % 26
        for r, v in enumerate(B.norm(signed).astype(np.uint8)):
            e = table.setdefault(v.tobytes(), [])
            if len(e) < 5:
                e.append((r, sg))
    log(f"distinct normalized crib-key patterns (with sign): {len(table):,}")
    alphs = list(zip(pc.strings(A), names))
    name_index = {n: i for i, n in enumerate(names)}
    dX = {n: a for a, n in alphs}
    partners = {"AZ": k4.AZ, "KA": k4.KA}
    # ---- kernel vs Python / sm_taskB.best_partial ----
    ok, msg = self_check(lib, SCU)
    log(f"partial self-check: {ok} ({msg})")
    if not ok:
        log("ABORT: partial kernel disagreed with sm_taskB.best_partial")
        sys.exit(1)
    # Speed probe on real Q3 Vigenere keys, not random keys. A 400-key random
    # probe measured 6.77 ms/key (full Q3 x 3 ciphertexts looked like 6.92 h) but
    # the sustained rate on all 613,440 keys was 14.37 ms/key, i.e. 14.7 h, which
    # is past the ~12 h line for running every phrase alphabet. Time a long slice
    # of the actual implied keys after a warmup so the estimate sees sustained clocks.
    ks_cal = B.implied_keys(k4.CT, alphs, partners)
    Kcal = (np.asarray(ks_cal["Q3vig"]) % 26).astype(np.uint8)
    ncal = min(8000, len(Kcal))
    partial(lib, Kcal[:200], SCU)
    t1 = time.time()
    partial(lib, Kcal[:ncal], SCU)
    dt = time.time() - t1
    per_key = dt / ncal
    del ks_cal, Kcal
    n_arr = 10  # Q3, Q1AZ, Q1KA, Q2AZ, Q2KA x vig/beau
    est_one_ct = per_key * len(names) * n_arr
    est_all = est_one_ct * 3
    log(f"partial probe {dt:.2f}s / {ncal} real Q3vig keys -> {per_key*1000:.2f} ms/key. "
        f"Full partial {n_arr} arrangements x 3 ciphertexts estimated {est_all/3600:.2f} h")
    log("prior full-alphabet K4 Q3 pass (stopped once the sustained rate was known): "
        "exact crib fits 0; best partial 13/24 "
        "('Q3vig', 'LOTOFTIMEGOING/s', 'bc:digsum start=17:30 step=-15', sigma=+1, b=14); "
        "hist_tail {11:6932, 12:220, 13:14}; 17634s at 14.37 ms/key. "
        "shuf0 exact on the full set was 0. See results/phr_clock_k4full.out.")
    # Partial scope, identical for K4 and both shuffles. Exact test is always all 10 arrangements.
    # Tier 1: all 10 arrangements. Tier 2: Q3 vig+beau on every phrase alphabet.
    # Tier 3: Q3 vig+beau on the first M phrase alphabets (gap_phr order), M chosen to fit in 6 h.
    # Tier 2 only if the measured Q3 x 3 ciphertext estimate is within 8 h. Above that
    # (the ~12 h case) a prefix is the run, and the unscored tail is estimated.
    partial_M = len(names)
    if est_all <= BUDGET_S:
        partial_labels = ["Q3vig", "Q3beau", "Q1AZvig", "Q1AZbeau", "Q1KAvig", "Q1KAbeau",
                          "Q2AZvig", "Q2AZbeau", "Q2KAvig", "Q2KAbeau"]
        log("Running partial on all 10 arrangements, all phrase alphabets.")
    else:
        est_q3 = per_key * len(names) * 2 * 3
        if est_q3 <= 8 * 3600:
            partial_labels = ["Q3vig", "Q3beau"]
            log(f"SUBSET: full partial estimate {est_all/3600:.2f} h exceeds 6 h. "
                f"Partial is Q3 vig and Q3 beau only (2/10 arrangements), all {len(names):,} phrase alphabets, "
                f"estimate {est_q3/3600:.2f} h. Exact test still covers all 10 arrangements.")
        else:
            partial_labels = ["Q3vig", "Q3beau"]
            partial_M = max(1, int(BUDGET_S / (per_key * 2 * 3)))
            partial_M = min(partial_M, len(names))
            est_m = per_key * partial_M * 2 * 3
            rest = len(names) - partial_M
            rest_h = per_key * rest * 2 * 3 / 3600.0
            other_h = per_key * len(names) * 8 * 3 / 3600.0
            log(f"SUBSET: Q3 partial on every phrase alphabet estimates {est_q3/3600:.2f} h, "
                f"past the ~12 h line. Partial is Q3 vig and Q3 beau on the first {partial_M:,} "
                f"phrase alphabets (gap_phr order), estimate {est_m/3600:.2f} h, same prefix for "
                f"K4 and both shuffles. Exact test still covers all 10 arrangements and every phrase alphabet.")
            log(f"estimate of the rest: remaining {rest:,} phrase alphabets x Q3 vig+beau x 3 ciphertexts "
                f"= {rest_h:.2f} h; the other 8 arrangements x every phrase alphabet x 3 ciphertexts "
                f"= {other_h:.2f} h. Neither is run.")
    hit_i = name_index.get("LOTOFTIMEGOING/s")
    log(f"full-K4 best alphabet LOTOFTIMEGOING/s index={hit_i} "
        f"({'inside' if hit_i is not None and hit_i < partial_M else 'outside'} the scored prefix of {partial_M:,})")
    # ---- planted exact recovery (phrase alphabet + a real stream) ----
    prng = np.random.default_rng(7)
    pt = fc.planted_plaintext(prng)
    def have(nm):
        return nm if nm in name_index else None
    plants = [("bc:lit start=13:37 step=5", have("BERLINCLOCK/s") or names[0], "Q3", 1, 3, "vig"),
              ("bc:HHMM start=07:15 step=15 phase=2", have("EASTNORTHEAST/s") or names[1], "Q1KA", -1, 11, "beau"),
              ("tz:-12..12 start=5 dir=-1", have("CLOCK/s") or names[2], "Q2AZ", 1, 0, "vig")]
    plant_ok = True
    for sname, an, arr, sg, b, mode in plants:
        X = dX[an]
        if arr == "Q3":
            ap = ac = X
        elif arr.startswith("Q1"):
            ap, ac = X, partners[arr[2:]]
        else:
            ap, ac = partners[arr[2:]], X
        st_ = S[sname].astype(int)
        pct = B.encrypt_stream(pt, st_, ap, ac, sg, b, mode)
        assert B.decrypt_stream(pct, st_, ap, ac, sg, b, mode) == pt
        ks = B.implied_keys(pct, alphs, partners)
        lab = f"{arr}{mode}"
        # Q3 keys are named Q3vig / Q3beau; Q1KA -> Q1KAvig
        hits_idx = B.lookup(ks[lab], table)
        rec = an in (alphs[r][1] for r in hits_idx) or name_index[an] in hits_idx
        # confirm plaintext via the stream table, not merely the alphabet index
        got = False
        if name_index[an] in hits_idx:
            k = ks[lab][name_index[an]]
            nk = B.norm(k).astype(np.uint8).tobytes()
            for (si, sg2) in table[nk]:
                b2 = int((k[0] - sg2 * SC[si, 0]) % 26)
                dec = B.decrypt_stream(pct, S[snames[si]].astype(int), ap, ac, sg2, b2, mode)
                if dec == pt:
                    got = True
                    break
        log(f"[planted] stream={sname} alph={an} {arr}{mode} sigma={sg} b={b}: alphabet hit={rec} plaintext recovered={got}")
        plant_ok &= got
    if not plant_ok:
        log("ABORT: planted clock message was not recovered")
        sys.exit(1)
    if os.environ.get("PHR_VALIDATE_ONLY"):
        log("PHR_VALIDATE_ONLY: stopping after planted recovery")
        return
    # free the per-stream dict before the big partial; keep SC, names, SCU
    del S
    sh = ["".join(rng.permutation(list(k4.CT))) for _ in range(2)]
    log(f"shuf0 {sh[0]}")
    log(f"shuf1 {sh[1]}")
    summaries = []
    wall0 = time.time()
    for nm, ct in (("K4", k4.CT), ("shuf0", sh[0]), ("shuf1", sh[1])):
        t1 = time.time()
        ks = B.implied_keys(ct, alphs, partners)
        nkeys = sum(len(v) for v in ks.values())
        hits = []
        for lab, K in ks.items():
            for r in B.lookup(K, table):
                k = K[r]
                nk = B.norm(k).astype(np.uint8).tobytes()
                for (si, sg) in table[nk]:
                    b = int((int(k[0]) - sg * int(SC[si, 0])) % 26)
                    hits.append((lab, alphs[r][1], snames[si], sg, b))
        log(f"[exact] {nm}: key vectors {nkeys:,} x {len(snames):,} streams x 2 signs x 26 offsets: "
            f"exact crib fits {len(hits)} ({time.time()-t1:.0f}s)")
        for h in hits[:20]:
            log(f"    HIT {h}")
        t2 = time.time()
        best = -1
        best_desc = None
        hist = np.zeros(25, np.uint64)
        for lab in partial_labels:
            K = (np.asarray(ks[lab]) % 26).astype(np.uint8)[:partial_M]
            h, b = partial(lib, K, SCU)
            hist += h
            if int(b[0]) > best:
                best = int(b[0])
                si = urep[int(b[2])]
                sg = 1 if int(b[3]) == 0 else -1
                best_desc = (lab, alphs[int(b[1])][1], snames[si], sg, int(b[4]), best)
        top = ", ".join(f"{s}:{int(hist[s])}" for s in range(max(0, best - 2), 25) if hist[s])
        log(f"[partial] {nm}: arrangements {partial_labels}: best {best}/24 {best_desc} hist_tail {{{top}}} ({time.time()-t2:.0f}s)")
        summaries.append((nm, len(hits), best, time.time() - t1))
    elapsed = time.time() - wall0
    cov_exact = 1.0
    cov_partial = (len(partial_labels) / 10.0) * (partial_M / len(names))
    k4s, s0, s1 = summaries
    log(f"SUMMARY family=clock combinations={len(names)}x{len(snames)}x10x2x26 "
        f"exact_k4={k4s[1]} exact_s0={s0[1]} exact_s1={s1[1]} "
        f"partial_k4={k4s[2]}/24 partial_s0={s0[2]}/24 partial_s1={s1[2]}/24 "
        f"partial_arr={','.join(partial_labels)} partial_M={partial_M} "
        f"k4_full_q3_partial=13/24 "
        f"planted=pass runtime_s={elapsed:.0f} coverage_exact={cov_exact:.4f} coverage_partial={cov_partial:.4f}")
    log(f"done {time.ctime()}")

if __name__ == "__main__":
    main()
