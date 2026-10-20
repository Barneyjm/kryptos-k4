"""Sanborn-phrase alphabets x interrupted / event-driven keyword keys.

Same models and Quagmire I-III arrangements as fe_interrupt.py (MIN 8 checks,
5 arrangements x vig/beau). Quagmire IV over phrase pairs is the same order of
search the dictionary run left undone (phrase x phrase); phrase x clue Quag IV
is estimated, not run. Partial score: majority key inside each equality class,
counted out of 24 crib letters. A score of 24 is an exact fit.
"""
import os, sys, time, ctypes, subprocess
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import phr_common as pc
import k4, fe_common as F
import fe_interrupt as I

RES = pc.RES
SO = os.path.join(RES, "phr_maj_kernel.so")
SRC = os.path.join(RES, "phr_maj_kernel.c")
LOG = os.path.join(RES, "phr_interrupt.out")
BUDGET_S = 6 * 3600

def build():
    subprocess.check_call(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC", "-o", SO, SRC])

def load_lib():
    if not os.path.exists(SO) or os.path.getmtime(SO) < os.path.getmtime(SRC):
        build()
    lib = ctypes.CDLL(SO)
    lib.maj_partial.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,
                                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,
                                ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p]
    lib.maj_partial.restype = None
    return lib

def pack_models(eng):
    """Union-find the engine's spanning edges back into classes of size >= 2."""
    parent = list(range(24))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    members = []
    cls_off = [0]
    model_off = [0]
    nfree = []
    for edges in eng.edges:
        parent[:] = range(24)
        for e in edges:
            a, b = F.PAIRS[e]
            union(a, b)
        groups = {}
        for i in range(24):
            groups.setdefault(find(i), []).append(i)
        nf = 0
        for mem in groups.values():
            if len(mem) == 1:
                nf += 1
            else:
                members.extend(mem)
                cls_off.append(len(members))
        model_off.append(cls_off.__len__() - 1)
        nfree.append(nf)
    return (np.array(members, np.uint8), np.array(cls_off, np.int32),
            np.array(model_off, np.int32), np.array(nfree, np.uint8))

def maj_py(key, members, cls_off, model_off, nfree, m):
    sc = int(nfree[m])
    for c in range(int(model_off[m]), int(model_off[m + 1])):
        cnt = np.zeros(26, np.int16)
        for j in range(int(cls_off[c]), int(cls_off[c + 1])):
            cnt[key[members[j]]] += 1
        sc += int(cnt.max())
    return sc

def key_matrix(ct, alph_u8, which):
    """which is one of I.scan's five arrangements, returning vig and beau (n,24) uint8."""
    n = len(alph_u8)
    A = alph_u8.astype(np.int16)
    INV = np.argsort(A, axis=1).astype(np.int16)
    AZi = np.tile(np.arange(26, dtype=np.int16), (n, 1))
    KAi = np.tile(np.array([k4.KA.index(chr(65 + x)) for x in range(26)], dtype=np.int16), (n, 1))
    combos = {"AA": (INV, INV), "AZA": (AZi, INV), "AAZ": (INV, AZi), "KAA": (KAi, INV), "AKA": (INV, KAi)}
    Pi, Ci = combos[which]
    C = np.array([ord(ch) - 65 for ch in ct], np.int16)[F.POS]
    pc_ = Pi[:, F.PL]
    cc = Ci[:, C]
    vig = ((cc - pc_) % 26).astype(np.uint8)
    beau = ((cc + pc_) % 26).astype(np.uint8)
    return vig, beau

def run_maj(lib, keys, pack, max_exact=200):
    members, cls_off, model_off, nfree = pack
    keys = np.ascontiguousarray(keys, np.uint8)
    hist = np.zeros(25, np.uint64)
    best = np.zeros(3, np.int32)
    ex = np.zeros(max_exact * 2, np.int32)
    nexact = ctypes.c_int(0)
    lib.maj_partial(keys.ctypes.data, len(keys), members.ctypes.data, cls_off.ctypes.data,
                    model_off.ctypes.data, nfree.ctypes.data, len(nfree), hist.ctypes.data, best.ctypes.data,
                    ex.ctypes.data, max_exact, ctypes.byref(nexact))
    pairs = [(int(ex[2 * i]), int(ex[2 * i + 1])) for i in range(nexact.value)]
    return hist, best, pairs

def self_check(lib, pack, names_m):
    members, cls_off, model_off, nfree = pack
    rng = np.random.default_rng(1)
    keys = rng.integers(0, 26, size=(40, 24), dtype=np.uint8)
    # force model 0 to be an exact fit on key 0, and a one-letter miss on key 1
    k0 = keys[0].copy()
    for c in range(int(model_off[0]), int(model_off[1])):
        mem = members[int(cls_off[c]):int(cls_off[c + 1])]
        k0[mem] = k0[mem[0]]
    keys[0] = k0
    keys[1] = k0.copy()
    # flip one constrained position
    c0 = int(model_off[0])
    pos = int(members[int(cls_off[c0])])
    keys[1, pos] = (int(keys[1, pos]) + 1) % 26
    hist, best, _ = run_maj(lib, keys, pack)
    ok_py = maj_py(keys[0], members, cls_off, model_off, nfree, 0) == 24
    ok_miss = maj_py(keys[1], members, cls_off, model_off, nfree, 0) == 23
    # C best on the 40 includes the planted 24
    ok = ok_py and ok_miss and int(best[0]) == 24 and int(hist[24]) >= 1
    # spot-check C vs Python on model 0 by restricting the pack to that model
    pack0 = (members, cls_off, model_off[:2].copy(), nfree[:1].copy())
    h0, b0, _ = run_maj(lib, keys[:8], pack0)
    py0 = [maj_py(keys[i], members, cls_off, model_off, nfree, 0) for i in range(8)]
    # hist from a one-model pack: count of each score across 8 keys
    ok_match = all(int(h0[s]) == sum(v == s for v in py0) for s in range(25))
    return ok and ok_match, f"model0 exact py={ok_py} miss={ok_miss} Cbest={int(best[0])} hist24={int(hist[24])} match={ok_match} example={names_m[0]}"

def main():
    build()
    lib = load_lib()
    f = open(LOG, "w")
    def log(s):
        print(s, flush=True)
        print(s, file=f, flush=True)
    A, names = pc.phrase_alphabets()
    log(f"phrase alphabets {len(names):,} sha256/16={pc.digest(A)} first={names[0]} last={names[-1]}")
    log("scope: fe_interrupt models with >= 8 checks; Quagmire I-III (AA, AZ/A, A/AZ, KA/A, A/KA) x vig/beau. "
        "Quag IV phrase x phrase not run (same omission as dictionary x dictionary). "
        "Quag IV phrase x 116 clue is estimated from the Session 2 matmul rate and not run.")
    t0 = time.time()
    AL = list(zip(pc.strings(A), names))
    log(f"alphabet strings built ({time.time()-t0:.0f}s)")
    log("planted-solution validation (fe_interrupt.planted, 4 mechanisms):")
    # planted() prints to stdout; mirror by wrapping is awkward, so require the boolean.
    ok = I.planted(8, AL)
    log(f"planted recovered all 4 mechanisms: {ok}")
    if not ok:
        log("ABORT: planted validation failed")
        sys.exit(1)
    eng_k4 = F.Engine(I.gen_models(k4.CT), 8)
    pack_k4 = pack_models(eng_k4)
    okc, msg = self_check(lib, pack_k4, eng_k4.names)
    log(f"[K4] majority-kernel self-check: {okc} ({msg}) models={len(eng_k4):,}")
    if not okc:
        log("ABORT: majority kernel disagreed with the Python scorer")
        sys.exit(1)
    # Partial-score probe. Exact Quag I-III always covers every alphabet; partial may be a subset.
    combos = ("AA", "AZA", "AAZ", "KAA", "AKA")
    full_plan = [(w, kn) for w in combos for kn in ("vig", "beau")]
    probe_n = 300
    rng_probe = np.random.default_rng(0)
    probe_keys = rng_probe.integers(0, 26, size=(probe_n, 24), dtype=np.uint8)
    tprobe = time.time()
    run_maj(lib, probe_keys, pack_k4)
    per_key = (time.time() - tprobe) / probe_n
    est_full = per_key * len(names) * len(full_plan) * 3
    global PARTIAL_PLAN, PARTIAL_M
    if est_full <= BUDGET_S:
        PARTIAL_PLAN, PARTIAL_M = full_plan, len(names)
        log(f"partial probe {per_key*1000:.2f} ms/key. All {len(full_plan)} arrangements x all alphabets "
            f"x 3 ciphertexts estimated {est_full/3600:.2f} h. Running that.")
    else:
        # Q3 is the AA arrangement (same mixed alphabet both sides), vig and beau.
        est_q3 = per_key * len(names) * 2 * 3
        if est_q3 <= BUDGET_S:
            PARTIAL_PLAN, PARTIAL_M = [("AA", "vig"), ("AA", "beau")], len(names)
            log(f"SUBSET partial: full estimate {est_full/3600:.2f} h exceeds 6 h. "
                f"Partial is AA vig and AA beau (2/10 arrangements, the same-alphabet Quagmire) "
                f"on all {len(names):,} phrase alphabets, estimate {est_q3/3600:.2f} h. "
                f"Exact test still covers all 10 arrangements.")
        else:
            PARTIAL_M = max(1, int(BUDGET_S / (per_key * 2 * 3)))
            PARTIAL_PLAN = [("AA", "vig"), ("AA", "beau")]
            log(f"SUBSET partial: AA vig/beau on the first {PARTIAL_M:,} phrase alphabets "
                f"(gap_phr order). Full estimate {est_full/3600:.2f} h, all-alphabet AA estimate {est_q3/3600:.2f} h. "
                f"Exact test still covers all 10 arrangements and every phrase alphabet.")
    if os.environ.get("PHR_VALIDATE_ONLY"):
        log("PHR_VALIDATE_ONLY: stopping after planted recovery")
        return
    rng = np.random.default_rng(11)
    sh = ["".join(rng.permutation(list(k4.CT))) for _ in range(2)]
    log(f"shuf0 {sh[0]}")
    log(f"shuf1 {sh[1]}")
    summaries = []
    wall0 = time.time()
    for nm, ct in (("K4", k4.CT), ("shuf0", sh[0]), ("shuf1", sh[1])):
        t1 = time.time()
        if nm == "K4":
            eng, pack = eng_k4, pack_k4
        else:
            eng = F.Engine(I.gen_models(ct), 8)
            pack = pack_models(eng)
        nvec = len(AL) * 10
        exp = I.expected_fp(eng, nvec)
        log(f"[{nm}] models {len(eng):,} (>=8 checks) expected exact fits by chance {exp:.3f}")
        if nm == "K4":
            # Score 24 is an exact fit. Check that against fe_interrupt.scan on a slice.
            sl = 4000
            sub = list(zip(pc.strings(A[:sl]), names[:sl]))
            t2 = time.time()
            hits = I.scan(ct, sub, eng, quag4=False)
            c24 = 0
            for which in ("AA", "AZA", "AAZ", "KAA", "AKA"):
                vig, beau = key_matrix(ct, A[:sl], which)
                for K in (vig, beau):
                    h, _, _ = run_maj(lib, K, pack)
                    c24 += int(h[24])
            log(f"[{nm}] slice cross-check n={sl}: scan hits {len(hits)} majority-24 {c24} ({time.time()-t2:.0f}s)")
            if c24 != len(hits):
                log("ABORT: majority score-24 disagreed with fe_interrupt.scan on the slice")
                sys.exit(1)
        t3 = time.time()
        best = -1
        best_desc = None
        hist = np.zeros(25, np.uint64)
        exact_n = 0
        shown = 0
        for which, kn in PARTIAL_PLAN:
            K = key_matrix(ct, A, which)[0 if kn == "vig" else 1][:PARTIAL_M]
            h, b, pairs = run_maj(lib, K, pack)
            hist += h
            exact_n += int(h[24])
            if int(h[24]) and shown < 20:
                for ki, mi in pairs:
                    if shown >= 20:
                        break
                    log(f"    HIT alph={names[ki]} {which} {kn} model={eng.names[mi]}")
                    shown += 1
            if int(b[0]) > best:
                best = int(b[0])
                best_desc = (which, kn, names[int(b[1])], eng.names[int(b[2])], int(b[0]))
        log(f"[{nm}] exact fits {exact_n}  partial best {best}/24 at {best_desc}  "
            f"plan={len(PARTIAL_PLAN)} arr x {PARTIAL_M:,} alph  "
            f"hist>=20: {{{', '.join(f'{s}:{int(hist[s])}' for s in range(20, 25) if hist[s])}}}  ({time.time()-t3:.0f}s)")
        summaries.append((nm, exact_n, best, int(hist[best]) if best >= 0 else 0, time.time() - t1))
    elapsed = time.time() - wall0
    # Quag IV phrase x clue estimate. Session 2: 842,054 alphabets x 10 streams in 314s on K4.
    rate = 842054 * 10 / 314.0
    q4 = 2 * len(names) * 116 * 2  # phrase/clue both directions x vig/beau
    est_q4 = q4 / rate * 3 / 3600.0
    log(f"estimate Quag IV phrase x 116 clue, both directions, vig/beau, 3 ciphertexts: "
        f"{q4 * 3:,.0f} key vectors, ~{est_q4:.1f} h at the Session 2 rate of {rate:,.0f} vectors/s. Not run.")
    k4s, s0, s1 = summaries
    # Exact fits and the partial score come from the same sweep, so they share a coverage.
    cov_exact = (len(PARTIAL_PLAN) / 10.0) * (PARTIAL_M / len(names))
    cov_partial = cov_exact
    log(f"SUMMARY family=interrupt combinations={len(names)}x{len(eng_k4)}x10 "
        f"exact_k4={k4s[1]} exact_s0={s0[1]} exact_s1={s1[1]} "
        f"partial_k4={k4s[2]}/24 partial_s0={s0[2]}/24 partial_s1={s1[2]}/24 "
        f"at_best_k4={k4s[3]} at_best_s0={s0[3]} at_best_s1={s1[3]} "
        f"planted=pass runtime_s={elapsed:.0f} coverage_exact={cov_exact:.4f} coverage_partial={cov_partial:.4f} "
        f"partial_plan={len(PARTIAL_PLAN)}x{PARTIAL_M} "
        f"coverage_note=QuagI-III_exact_complete QuagIV_phrase_x_clue_not_run")
    log(f"done {time.ctime()}")

if __name__ == "__main__":
    main()
