"""T1: Interrupted / gated keyword keys (Friedman, Military Cryptanalytics II: 'interrupted keys').
A keyword of length L (or a non-repeating key, L=inf) is applied, but the key position is driven by events:
  ct_restart_after(T): key restarts at keyword[0] after each ciphertext letter in T (decryptor sees C)
  ct_restart_at(T):    the trigger letter itself takes keyword[0]
  ct_advance(T):       key advances only after a ciphertext letter in T
  ct_pause(T):         key does not advance after ciphertext letters in T
  pt_restart(t):       key restarts after each plaintext letter t (positions before the first t inside a crib are unconstrained)
  pt_advance/pt_pause(t), phase-slip: key driven by plaintext letters; unknown offset delta between the two cribs (enumerated)
  word_restart:        key restarts at each word (crib segmentations EAST|NORTH|EAST etc.)
  line_restart:        key restarts at each carved line of the sculpture (K4 line starts at 0-idx 4, 35, 66)
T = single letters (26) and letter pairs (325). L = 1..26 and inf.
Alphabets: Quagmire I/II/III (+KA) over 842k dictionary + 116 clue alphabets; Quagmire IV over all clue pairs.
Exact test; shuffled-CT runs are reported only to calibrate the chance hit rate.
"""
import sys, time, itertools, numpy as np, k4
import fe_common as F

LMAX = 26
LS = list(range(1, LMAX + 1)) + [None]   # None = non-repeating key (offset itself)
CRIBS = [(21, "EASTNORTHEAST"), (63, "BERLINCLOCK")]
LINE_STARTS = [0, 4, 35, 66]

def labels_from_idx(idx, L, tag=""):
    """idx: dict pos->int or None. Returns 24 labels."""
    out = []
    for p in F.POS:
        v = idx.get(p)
        out.append(None if v is None else (tag, v if L is None else v % L))
    return out

def ct_index(ct, T, kind):
    idx = {}; last = -1; cnt_in = 0; cnt_out = 0
    for i, c in enumerate(ct):
        if kind == "restart_after": idx[i] = i - last - 1
        elif kind == "advance": idx[i] = cnt_in
        elif kind == "pause": idx[i] = cnt_out
        if c in T:
            last = i; cnt_in += 1
        else:
            cnt_out += 1
        if kind == "restart_at": idx[i] = (i - last) if last >= 0 else i
    return idx

def gen_models(ct):
    letters = [chr(65 + i) for i in range(26)]
    trig = [frozenset([a]) for a in letters] + [frozenset(p) for p in itertools.combinations(letters, 2)]
    for T in trig:
        tn = "".join(sorted(T))
        for kind in ("restart_after", "restart_at", "advance", "pause"):
            idx = ct_index(ct, T, kind)
            for L in LS:
                yield (f"ct_{kind} T={tn} L={L}", labels_from_idx(idx, L))
    # plaintext-driven
    crib_letters = sorted(set(k4.CPL))
    for t in crib_letters:
        idx = {}
        for s, w in CRIBS:
            last = None
            for j, ch in enumerate(w):
                idx[s + j] = None if last is None else (s + j) - last - 1
                if ch == t: last = s + j
        for L in LS:
            yield (f"pt_restart_after t={t} L={L}", labels_from_idx(idx, L))
    for t in crib_letters + ["*"]:
        for kind in ("advance", "pause"):
            if t == "*" and kind == "advance": continue
            rel = {}
            for s, w in CRIBS:
                n = 0
                for j, ch in enumerate(w):
                    rel[s + j] = n
                    hit = (ch == t)
                    n += (1 if hit else 0) if kind == "advance" else (0 if hit else 1)
            for L in range(1, LMAX + 1):
                for d in range(L):
                    lab = [((rel[p] + (d if p >= 63 else 0)) % L) for p in F.POS]
                    yield (f"pt_{kind} t={t} L={L} delta={d}" if t != "*" else f"phase_slip L={L} delta={d}", lab)
    # word restart: segmentations; first word in each crib either starts at the crib start (offset 0) or unknown
    seg1 = {"EAST|NORTH|EAST": [0, 4, 9], "EAST|NORTHEAST": [0, 4], "EASTNORTHEAST": [0]}
    seg2 = {"BERLIN|CLOCK": [0, 6], "BERLINCLOCK": [0]}
    for (n1, b1), (n2, b2) in itertools.product(seg1.items(), seg2.items()):
        for known1, known2 in itertools.product((True, False), repeat=2):
            idx = {}
            for (s, w), bounds, known in ((CRIBS[0], b1, known1), (CRIBS[1], b2, known2)):
                for j in range(len(w)):
                    b = max(x for x in bounds if x <= j)
                    idx[s + j] = None if (b == 0 and not known) else j - b
            for L in LS:
                yield (f"word_restart {n1}{'' if known1 else '(?)'} {n2}{'' if known2 else '(?)'} L={L}", labels_from_idx(idx, L))
    idx = {p: p - max(b for b in LINE_STARTS if b <= p) for p in F.POS}
    for L in LS:
        yield (f"line_restart L={L}", labels_from_idx(idx, L))

def key_vectors(ct, PAi, CAi):
    """PAi, CAi: (n,26) position-of-letter arrays. Returns vig, beau implied keys (n,24)."""
    C = np.array([ord(ch) - 65 for ch in ct])[F.POS]
    pc = PAi[:, F.PL]; cc = CAi[:, C]
    return (cc - pc) % 26, (cc + pc) % 26

def scan(ct, alphs, eng, quag4=False):
    """alphs: list of (alphabet, name). Returns hits list."""
    A = F.alph_array([a for a, _ in alphs]); INV = F.inv_index(A).astype(np.int16)
    n = len(A)
    AZi = np.tile(np.arange(26, dtype=np.int16), (n, 1))
    KAi = np.tile(np.array([k4.KA.index(chr(65 + x)) for x in range(26)], dtype=np.int16), (n, 1))
    hits = []
    combos = [("A", INV, "A", INV), ("AZ", AZi, "A", INV), ("A", INV, "AZ", AZi), ("KA", KAi, "A", INV), ("A", INV, "KA", KAi)]
    for pn, Pi, cn, Ci in combos:
        for kn, K in zip(("vig", "beau"), key_vectors(ct, Pi, Ci)):
            for r, m in eng.fits(K):
                hits.append((alphs[r][1], pn, cn, kn, eng.names[m], int(eng.checks[m])))
    if quag4:   # all ordered pairs of alphabets (small sets only)
        for i in range(n):
            Pi = np.repeat(INV[i:i + 1], n, 0)
            for kn, K in zip(("vig", "beau"), key_vectors(ct, Pi, INV)):
                for r, m in eng.fits(K):
                    hits.append((alphs[i][1] + "|" + alphs[r][1], "P", "C", kn, eng.names[m], int(eng.checks[m])))
    return hits

def expected_fp(eng, nvec):
    return float(sum(nvec * 26.0 ** (-c) for c in eng.checks))

def encrypt(pt, PA, CA, kw, mech, T=None, L=None, mode="vig"):
    out = []; idxf = None
    last = -1; cin = cout = 0
    for i, ch in enumerate(pt):
        if mech == "restart_after": ki = i - last - 1
        elif mech == "advance": ki = cin
        elif mech == "pause": ki = cout
        elif mech == "line_restart": ki = i - max(b for b in LINE_STARTS if b <= i)
        k = kw[ki % len(kw)]
        p = PA.index(ch)
        c = CA[(p + k) % 26] if mode == "vig" else CA[(k - p) % 26]
        out.append(c)
        if T is not None:
            if c in T: last = i; cin += 1
            else: cout += 1
    return "".join(out)

def planted(MIN, AL):
    rng = np.random.default_rng(7); ok = True
    cases = [("restart_after", frozenset("K"), 7, "vig", "A", "A"), ("advance", frozenset("EU"), 5, "beau", "AZ", "A"),
             ("pause", frozenset("Q"), 11, "vig", "A", "KA"), ("line_restart", None, 9, "beau", "A", "A")]
    for mech, T, L, mode, pside, cside in cases:
        j = int(rng.integers(len(AL))); X = AL[j][0]
        PA = X if pside == "A" else (k4.AZ if pside == "AZ" else k4.KA); CA = X if cside == "A" else (k4.AZ if cside == "AZ" else k4.KA)
        kw = [int(v) for v in rng.integers(0, 26, L)]
        pt = F.planted_plaintext(rng); ct = encrypt(pt, PA, CA, kw, mech, T, L, mode)
        eng = F.Engine(gen_models(ct), MIN)
        sub = [AL[x] for x in rng.choice(len(AL), 20000, replace=False)] + [AL[j]]
        h = scan(ct, sub, eng)
        want = (f"ct_{mech} T={''.join(sorted(T))} L={L}" if T else f"line_restart L={L}")
        # the planted model may have been deduplicated under another name with identical edges
        tgt_edges = F.model_from_labels(dict(gen_named(ct)).get(want))
        names_with_edges = {eng.names[m] for m in range(len(eng)) if eng.edges[m] == tgt_edges}
        found = [x for x in h if x[0] == AL[j][1] and x[4] in names_with_edges and x[1] == pside and x[2] == cside and x[3] == mode]
        print(f"  planted {want} {mode} P={pside} C={cside} alph={AL[j][1]}: recovered={bool(found)}  "
              f"(model checks={len(tgt_edges)}, other hits in 20k decoys={len(h) - len(found)})", flush=True)
        ok &= bool(found)
    return ok

def gen_named(ct):
    return list(gen_models(ct))

def run_full(ct, AL, CL, eng, label):
    t = time.time()
    h = scan(ct, AL, eng)
    h4 = scan(ct, CL, eng, quag4=True)
    nvec = len(AL) * 10 + len(CL) ** 2 * 2
    print(f"{label}: {len(h)} Quag I-III hits + {len(h4)} Quag IV(clue) hits; expected by chance ~{expected_fp(eng, nvec):.2f}  ({time.time()-t:.0f}s)", flush=True)
    return h + h4

if __name__ == "__main__":
    MIN = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    AL = F.all_alphabets(); CL = F.clue_alphabets()
    eng = F.Engine(gen_models(k4.CT), MIN)
    allm = gen_named(k4.CT)
    print(f"T1 interrupted keys: {len(allm):,} raw models -> {len(eng):,} distinct with >= {MIN} checks; alphabets {len(AL):,} (+{len(CL)} clue for Quag IV)", flush=True)
    print("planted-solution validation:", flush=True)
    assert planted(MIN, AL), "planted test failed"
    h = run_full(k4.CT, AL, CL, eng, "K4")
    for x in sorted(h, key=lambda x: -x[5])[:40]: print("   ", x, flush=True)
    rng = np.random.default_rng(11)
    for t in range(2):
        sh = "".join(rng.permutation(list(k4.CT)))
        e2 = F.Engine(gen_models(sh), MIN)
        hs = run_full(sh, AL, CL, e2, f"shuffled{t}")
