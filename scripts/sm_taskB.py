"""Task B: clock-derived additive key streams (Berlin clock, HHMM digits, time zones, K2 coords/dates)
against AZ/KA/clue/dictionary alphabets. Exact crib fits via hashing of normalized implied keys.
Model: idxC(C_i) = sigma*stream_i + b  -/+ idxP(P_i)  (mod 26), sigma in {+1,-1}, b in Z26,
  vig: k = idxC(C) - idxP(P);  beau: k = idxC(C) + idxP(P).
Arrangements: Q3 (same mixed alphabet both sides), Q1 (plain mixed / cipher AZ or KA), Q2 (plain AZ or KA / cipher mixed),
  Q4 over clue x clue pairs."""
import time, sys, numpy as np, k4, fe_common as fc
POS = np.array(k4.CPOS); PL = k4.CPL; N = 97
SIGMAS = (1, -1)

# ---------------- stream library ----------------
def berlin(t):
    t %= 1440; h, m = divmod(t, 60)
    return h // 5, h % 5, m // 5, m % 5, h, m
def lamp_bits(t):
    r1, r2, r3, r4, h, m = berlin(t)
    sec = [1]
    return sec + [1] * r1 + [0] * (4 - r1) + [1] * r2 + [0] * (4 - r2) + [1] * r3 + [0] * (11 - r3) + [1] * r4 + [0] * (4 - r4)

def build_streams():
    S = {}              # name -> int array (97)
    def add(name, arr):
        S[name] = (np.array(arr[:N], dtype=np.int64) % 26).astype(np.uint8)
    steps = (1, 5, 15, 30, 60, -1, -5, -15, -30, -60)
    scal = {
        "lit": lambda r: r[0] + r[1] + r[2] + r[3],
        "unlit": lambda r: 23 - (r[0] + r[1] + r[2] + r[3]),
        "r1": lambda r: r[0], "r2": lambda r: r[1], "r3": lambda r: r[2], "r4": lambda r: r[3],
        "r3+r4": lambda r: r[2] + r[3], "r1+r2": lambda r: r[0] + r[1],
        "red": lambda r: r[0] + r[1] + r[2] // 3, "yellow": lambda r: r[2] - r[2] // 3 + r[3],
        "base5": lambda r: r[0] * 5 ** 3 + r[1] * 5 ** 2 + r[2] * 5 + r[3],   # mixed base, value mod 26 used
        "h": lambda r: r[4], "m": lambda r: r[5], "minuteofday": lambda r: r[4] * 60 + r[5], "h+m": lambda r: r[4] + r[5],
        "hhmm": lambda r: r[4] * 100 + r[5], "h12": lambda r: (r[4] - 1) % 12 + 1,
        "digsum": lambda r: sum(int(c) for c in f"{r[4]:02d}{r[5]:02d}"),
    }
    groups = {
        "rows_td": lambda r: [r[0], r[1], r[2], r[3]], "rows_bu": lambda r: [r[3], r[2], r[1], r[0]],
        "h,m": lambda r: [r[4], r[5]], "m,h": lambda r: [r[5], r[4]],
        "HHMM": lambda r: [int(c) for c in f"{r[4]:02d}{r[5]:02d}"],
        "HHMM12": lambda r: [int(c) for c in f"{(r[4]-1)%12+1:02d}{r[5]:02d}"],
        "MMHH": lambda r: [int(c) for c in f"{r[5]:02d}{r[4]:02d}"],
        "HHMM_ns": lambda r: [int(c) for c in f"{r[4]}{r[5]:02d}"],     # no leading zero on hour
    }
    for st in steps:
        for t0 in range(1440):
            times = [t0 + st * i for i in range(N)]
            rs = [berlin(t) for t in times]
            for nm, f in scal.items():
                add(f"bc:{nm} start={t0//60:02d}:{t0%60:02d} step={st}", [f(r) for r in rs])
            for nm, f in groups.items():
                g = len(f(rs[0])); flat = []
                for i in range(N // 1 + 2):
                    flat += f(berlin(t0 + st * i))
                    if len(flat) > N + g: break
                for ph in range(g):
                    add(f"bc:{nm} start={t0//60:02d}:{t0%60:02d} step={st} phase={ph}", flat[ph:ph + N])
            # lamp bit stream (top-down, row by row; and bottom-up)
            for order in ("td", "bu"):
                flat = []
                i = 0
                while len(flat) < N + 24:
                    b = lamp_bits(t0 + st * i); flat += b if order == "td" else b[::-1]; i += 1
                for ph in range(24):
                    add(f"bc:lampbits_{order} start={t0//60:02d}:{t0%60:02d} step={st} phase={ph}", flat[ph:ph + N])
    # time zones
    zl = {"-12..12": list(range(-12, 13)), "-12..11": list(range(-12, 12)), "-11..12": list(range(-11, 13))}
    for nm, z in zl.items():
        for d in (1, -1, 2, -2, 5, 7, 11):
            for s in range(len(z)):
                add(f"tz:{nm} start={z[s]} dir={d}", [z[(s + d * i) % len(z)] for i in range(N)])
                add(f"tz:abs{nm} start={z[s]} dir={d}", [abs(z[(s + d * i) % len(z)]) for i in range(N)])
    # K2 coordinates and dates, cyclic digit / number streams
    digs = {"coords": "38576577844", "lat": "385765", "lon": "77844",
            "coords_nodec": "3857677844", "1986": "1986", "1989": "1989", "1990": "1990",
            "19861989": "19861989", "198619891990": "198619891990", "11091989": "11091989", "09111989": "09111989",
            "1191989": "1191989", "9111989": "9111989", "110989": "110989", "091189": "091189", "19891109": "19891109",
            "11031990": "11031990", "03111990": "03111990", "1131990": "1131990", "3111990": "3111990", "19901103": "19901103",
            "1986-1989-1990": "198619891990", "coords+1990": "38576577844" + "1990"}
    for nm, d in digs.items():
        v = [int(c) for c in d]
        for s in range(len(v)):
            add(f"dig:{nm} phase={s}", [v[(s + i) % len(v)] for i in range(N)])
            add(f"digrev:{nm} phase={s}", [v[::-1][(s + i) % len(v)] for i in range(N)])
    nums = {"coordnums": [38, 57, 6, 5, 77, 8, 44], "coordnums2": [38, 57, 65, 77, 8, 44], "dms": [38, 57, 6, 77, 8, 44],
            "years": [1986, 1989, 1990], "wall": [11, 9, 1989], "wall2": [9, 11, 1989], "ded": [11, 3, 1990], "ded2": [3, 11, 1990]}
    for nm, v in nums.items():
        for s in range(len(v)):
            add(f"num:{nm} phase={s}", [v[(s + i) % len(v)] for i in range(N)])
    return S

def norm(v):           # v: (...,24) ints -> normalized mod-26 vector relative to first entry
    v = np.asarray(v) % 26
    return (v - v[..., :1]) % 26

# ---------------- alphabets / implied keys ----------------
def idx_of(alph, s):
    return np.array([alph.index(c) for c in s])

def implied_keys(ct, alphs, partners):
    """yields (label, (24,) implied key) for Q3 / Q1 / Q2 with each partner, vig and beau."""
    cc = [ct[p] for p in POS]
    A = np.array([[ord(c) - 65 for c in a] for a, _ in alphs], np.int16)
    inv = np.argsort(A, axis=1)                                     # inv[a, letter] = index of letter in alphabet a
    ci = inv[:, [ord(c) - 65 for c in cc]]; pi = inv[:, [ord(c) - 65 for c in PL]]
    out = {"Q3vig": ci - pi, "Q3beau": ci + pi}
    for pn, P in partners.items():
        cpi = idx_of(P, cc)[None, :]; ppi = idx_of(P, PL)[None, :]
        out[f"Q1{pn}vig"] = cpi - pi; out[f"Q1{pn}beau"] = cpi + pi        # plain mixed, cipher partner
        out[f"Q2{pn}vig"] = ci - ppi; out[f"Q2{pn}beau"] = ci + ppi        # plain partner, cipher mixed
    return out

def q4_keys(ct, clue):
    cc = [ct[p] for p in POS]
    A = np.array([[ord(c) - 65 for c in a] for a, _ in clue], np.int16); inv = np.argsort(A, axis=1)
    ci = inv[:, [ord(c) - 65 for c in cc]]; pi = inv[:, [ord(c) - 65 for c in PL]]
    n = len(clue)
    return {"Q4vig": (ci[None, :, :] - pi[:, None, :]).reshape(n * n, 24),   # row = pa*n + ca
            "Q4beau": (ci[None, :, :] + pi[:, None, :]).reshape(n * n, 24)}

def lookup(keys, table):
    """keys: (n,24) -> list of row indices whose normalized key is in table"""
    nk = norm(keys).astype(np.uint8)
    return [r for r in range(len(nk)) if nk[r].tobytes() in table]

def best_partial(keys, SC):
    """max over streams, sigma, b of #crib positions with key == sigma*s+b. keys (m,24), SC (S,24). returns (best, argstream)"""
    best = 0; arg = None
    for sg in SIGMAS:
        for r, k in enumerate(keys):
            D = np.sort((k[None, :] - sg * SC) % 26, axis=1).astype(np.int8)
            # longest run of equal values per row
            eq = np.diff(D, axis=1) == 0
            # run length via cumulative trick
            run = np.zeros(D.shape[0], np.int64); cur = np.ones(D.shape[0], np.int64)
            for j in range(23):
                cur = np.where(eq[:, j], cur + 1, 1); run = np.maximum(run, cur)
            run = np.maximum(run, 1); m = int(run.max())
            if m > best: best = m; arg = (r, sg, int(run.argmax()))
    return best, arg

def encrypt_stream(pt, stream, alph_p, alph_c, sigma, b, mode):
    out = []
    for i, p in enumerate(pt):
        k = sigma * stream[i] + b
        j = (k + alph_p.index(p)) if mode == "vig" else (k - alph_p.index(p))
        out.append(alph_c[j % 26])
    return "".join(out)

def decrypt_stream(ct, stream, alph_p, alph_c, sigma, b, mode):
    out = []
    for i, c in enumerate(ct):
        k = sigma * stream[i] + b
        j = (alph_c.index(c) - k) if mode == "vig" else (k - alph_c.index(c))
        out.append(alph_p[j % 26])
    return "".join(out)

if __name__ == "__main__":
    t0 = time.time()
    S = build_streams(); names = list(S)
    SC = np.array([S[n][POS] for n in names]).astype(np.int16)
    SCU, uinv = np.unique(SC, axis=0, return_inverse=True); uinv = uinv.ravel()
    urep = {}
    for i, u in enumerate(uinv): urep.setdefault(int(u), i)
    print(f"streams: {len(names):,}, distinct crib patterns {len(SCU):,}  ({time.time()-t0:.0f}s)", flush=True)
    table = {}
    for sg in SIGMAS:
        for r, v in enumerate(norm(sg * SC).astype(np.uint8)):
            e = table.setdefault(v.tobytes(), [])
            if len(e) < 5: e.append((r, sg))
    print(f"distinct normalized crib-key patterns (with sign): {len(table):,}", flush=True)
    clue = fc.clue_alphabets(); alls = fc.all_alphabets()
    partners = {"AZ": k4.AZ, "KA": k4.KA}
    rng = np.random.default_rng(0); sh = ["".join(rng.permutation(list(k4.CT))) for _ in range(2)]
    f = open("sm_taskB.out", "a"); print(f"=== Task B run {time.ctime()}  streams {len(names):,}, patterns {len(table):,} ===", file=f)

    def full_test(tag, ct, verbose=True):
        hits = []
        ks = implied_keys(ct, alls, partners)
        for lab, K in ks.items():
            for r in lookup(K, table):
                hits.append((lab, alls[r][1], r, K[r]))
        q4 = q4_keys(ct, clue)
        for lab, K in q4.items():
            for r in lookup(K, table):
                pa, ca = divmod(r, len(clue)); hits.append((lab, f"P={clue[pa][1]} C={clue[ca][1]}", r, K[r]))
        out = []
        for lab, an, r, k in hits:
            nk = norm(k).astype(np.uint8).tobytes()
            for (si, sg) in table[nk]:
                b = int((k[0] - sg * SC[si, 0]) % 26)
                out.append((lab, an, names[si], sg, b, len(table[nk])))
        return out, sum(len(v) for v in ks.values()) + sum(len(v) for v in q4.values())

    def decode(ct, hit):
        lab, an, sname, sg, b, _ = hit
        mode = "vig" if lab.endswith("vig") else "beau"
        if lab.startswith("Q4"):
            pn, cn = an.split(" "); d = dict((n, a) for a, n in clue)
            ap, ac = d[pn[2:]], d[cn[2:]]
        else:
            d = dict((n, a) for a, n in alls); X = d[an]
            if lab.startswith("Q3"): ap = ac = X
            elif lab.startswith("Q1"): ap, ac = X, partners[lab[2:4]]
            else: ap, ac = partners[lab[2:4]], X
        return decrypt_stream(ct, S[sname].astype(int), ap, ac, sg, b, mode)

    # ---- planted validation ----
    prng = np.random.default_rng(7)
    pt = fc.planted_plaintext(prng)
    dX = dict((n, a) for a, n in alls)
    plants = [("bc:lit start=13:37 step=5", "CHINAMAN/s", "Q3", 1, 3, "vig"),
              ("bc:HHMM start=07:15 step=15 phase=2", "KRYPTOS/t-rev" if "KRYPTOS/t-rev" in dX else "CHINAMAN/s", "Q1KA", -1, 11, "beau"),
              ("tz:-12..12 start=5 dir=-1", "kw:PALIMPSEST", "Q2AZ", 1, 0, "vig")]
    for sname, an, arr, sg, b, mode in plants:
        X = dX[an] if an in dX else dict((n, a) for a, n in clue)[an]
        ap, ac = (X, X) if arr == "Q3" else ((X, partners[arr[2:]]) if arr.startswith("Q1") else (partners[arr[2:]], X))
        st_ = S[sname].astype(int); pct = encrypt_stream(pt, st_, ap, ac, sg, b, mode)
        assert decrypt_stream(pct, st_, ap, ac, sg, b, mode) == pt
        hits, nkeys = full_test("plant", pct)
        rec = [h for h in hits if decode(pct, h) == pt]
        print(f"[planted] stream={sname} alph={an} {arr}{mode} sigma={sg} b={b}: {len(hits)} hits, "
              f"planted plaintext recovered: {bool(rec)}  e.g. {rec[:1]}", file=f, flush=True)
    # ---- real ----
    for nm, ct in (("K4", k4.CT), ("shuf0", sh[0]), ("shuf1", sh[1])):
        t1 = time.time()
        hits, nkeys = full_test(nm, ct)
        print(f"[exact] {nm}: key vectors tested {nkeys:,} x {len(names):,} streams x 2 signs x 26 offsets: "
              f"exact crib fits {len(hits)}  ({time.time()-t1:.0f}s)", file=f, flush=True)
        for h in hits[:20]:
            print(f"      HIT {h}  pt={decode(ct, h)}", file=f, flush=True)
        # partial calibration: clue alphabets, Q3 vig/beau, all streams, best sigma/b
        ks = implied_keys(ct, clue, {})
        for lab in ("Q3vig", "Q3beau"):
            bp, arg = best_partial(ks[lab] % 26, SCU)
            r, sg, si = arg; si = urep[si]
            print(f"[partial] {nm} {lab} over {len(clue)} clue alphabets x {len(names):,} streams: best {bp}/24 "
                  f"({clue[r][1]}, sigma={sg}, {names[si]})", file=f, flush=True)
    print(f"done {time.time()-t0:.0f}s", file=f, flush=True)
