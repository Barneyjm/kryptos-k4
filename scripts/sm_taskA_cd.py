"""Task A2: Chaocipher, 116 clue alphabets x all 842,054 dictionary+clue alphabets, both disk assignments,
26 relative offsets, both directions. Scoring kernel in C (sm_chao_kernel.c), cross-checked against sm_taskA.run."""
import os, sys, time, ctypes
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np, multiprocessing as mp
import k4, fe_common as fc
import sm_taskA as A
from sm_chao import encrypt, encrypt_swapped, decrypt, decrypt_swapped

HERE = os.path.dirname(os.path.abspath(__file__))
THR = 9; MAXREC = 200000; CHUNK = 1000; NPROC = 4
_lib = None
def lib():
    global _lib
    if _lib is None:
        _lib = ctypes.CDLL(os.path.join(HERE, "sm_chao_kernel.so"))
        P = ctypes.c_void_p
        _lib.sweep.argtypes = [P, ctypes.c_int, P, ctypes.c_int, P, P, ctypes.c_int, P, P, P, ctypes.c_int]
        _lib.sweep.restype = ctypes.c_int
    return _lib

NEED = np.full(97, 255, np.uint8); NEED[k4.CPOS] = A.PLC
def csweep(clue, dic, ct):
    clue = np.ascontiguousarray(clue, np.uint8); dic = np.ascontiguousarray(dic, np.uint8); ct = np.ascontiguousarray(ct, np.uint8)
    ht = np.zeros(25, np.int64); hb = np.zeros(14, np.int64); rec = np.zeros((MAXREC, 7), np.int32)
    nr = lib().sweep(clue.ctypes.data, len(clue), dic.ctypes.data, len(dic), ct.ctypes.data, NEED.ctypes.data, THR,
                     ht.ctypes.data, hb.ctypes.data, rec.ctypes.data, MAXREC)
    return ht, hb, rec[:nr].copy()

def job(args):
    ctname, ct, s, e = args
    ht, hb, rec = csweep(A.G["C"], A.G["D"][s:e], ct)
    rec[:, 1] += s
    return ctname, e - s, ht, hb, rec

def a2s(v): return "".join(chr(65 + x) for x in v)
def describe(r):
    a, d, asg, off, sw, tot, b1 = (int(x) for x in r)
    cn, dn = A.G["cnames"][a], A.G["names"][d]
    L, R = (cn, dn) if asg == 0 else (dn, cn)
    return f"({tot}, {b1}, '{'swap' if sw else 'std'} L={L} R={R} rotR={off}')"
def row_pt(r, ct):
    a, d, asg, off, sw = (int(x) for x in r[:5])
    C, D = a2s(A.G["C"][a]), a2s(A.G["D"][d])
    L, R = (C, D) if asg == 0 else (D, C); R = R[off:] + R[:off]
    return (decrypt_swapped if sw else decrypt)(a2s(ct), L, R)

def summary(tag, n, ht, hb, recs, ctx, f, extra=""):
    mx = int(np.nonzero(ht)[0].max()); ex = [r for r in recs if r[5] == 24]
    line = (f"[{tag}] rows {n:,}  exact(24/24) {len(ex)}  best partial {mx}/24 (#rows at best {int(ht[mx])})  "
            f"best EASTNORTHEAST-block {int(np.nonzero(hb)[0].max())}/13  "
            f"hist>=7: {{{', '.join(f'{k}:{int(ht[k])}' for k in range(7, 25) if ht[k])}}}{extra}")
    print(line, file=f, flush=True)
    return line, ex

if __name__ == "__main__":
    A.init()
    nc, nd = len(A.G["C"]), len(A.G["D"])
    f = open("sm_taskA_cd.out", "a")
    print(f"=== Task A2 (clue x dict Chaocipher) {time.ctime()}: {nc} clue x {nd:,} alphabets x 2 assign x 26 off x 2 dir "
          f"= {nc*nd*104:,} rows per ciphertext ===", file=f, flush=True)
    rng = np.random.default_rng(0); sh = ["".join(rng.permutation(list(k4.CT))) for _ in range(2)]
    CT4 = A.codes(k4.CT)
    # ---- 1. cross-check C kernel against the numpy implementation (sm_taskA.run) ----
    D5 = A.G["D"][1000:1004]; C = A.G["C"]
    Ls, Rs = [], []
    for d in range(len(D5)):
        for a in range(nc):
            for asg in range(2):
                for off in range(26):
                    L, R = (C[a], D5[d]) if asg == 0 else (D5[d], C[a])
                    Ls.append(L); Rs.append(R[A.ROT[off]])
    Ls, Rs = np.array(Ls), np.array(Rs)
    tq = []; bq = []
    for sw in (False, True):
        t_, b_ = A.run(CT4, Ls, Rs, sw); tq.append(t_); bq.append(b_)
    tq = np.stack(tq, 1).ravel(); bq = np.stack(bq, 1).ravel()           # C order: sw innermost
    ht, hb, rec = csweep(C, D5, CT4)
    ok = (np.bincount(tq, minlength=25) == ht).all() and (np.bincount(bq, minlength=14) == hb).all()
    # also row-level check of recorded rows (thr 9 may record none; check via thr-free compare on hist only plus clue x clue)
    htc, hbc, _ = csweep(C, C, CT4)
    print(f"cross-check C vs numpy on {len(tq):,} rows: histograms identical = {ok}; clue x clue via C: hist>=7 "
          f"{ {k: int(htc[k]) for k in range(7, 25) if htc[k]} } (expect exactly 2x the earlier [clue] K4 run {{7:16, 8:2}})",
          file=f, flush=True)
    t0 = time.time(); csweep(C, A.G["D"][:200], CT4); dt = time.time() - t0
    rate = 200 * nc * 104 / dt
    est = nc * nd * 104 * 3 / rate / NPROC / 3600
    print(f"kernel speed {rate/1e6:.2f} M rows/s/process -> estimated {est:.1f} h for K4 + 2 shuffles with {NPROC} processes",
          file=f, flush=True)
    # ---- 2. planted: clue disk x dictionary disk ----
    prng = np.random.default_rng(99); pt = fc.planted_plaintext(prng)
    di = 456789; ci = A.G["cnames"].index("kw:WELTZEITUHR")
    Dd, Cc = a2s(A.G["D"][di]), a2s(C[ci])
    plants = [("plantA2-std", encrypt(pt, Cc, Dd[17:] + Dd[:17])),            # L=clue, R=dict rot 17, std
              ("plantA2-swap", encrypt_swapped(pt, Dd, Cc[4:] + Cc[:4]))]      # L=dict, R=clue rot 4, swap
    for nm, pct in plants:
        ht, hb, rec = csweep(C, A.G["D"][di - 2:di + 3], A.codes(pct)); rec[:, 1] += di - 2
        line, ex = summary(f"planted {nm}", 5 * nc * 104, ht, hb, rec, None, f)
        recovered = any(row_pt(r, A.codes(pct)) == pt for r in ex)
        print(f"      planted keys: clue={A.G['cnames'][ci]} dict={A.G['names'][di]}; exact rows: "
              f"{[describe(r) for r in ex]}; planted plaintext recovered: {recovered}", file=f, flush=True)
    if est > 9:
        print(f"NOTE: estimate {est:.1f} h exceeds 9 h; running in order K4, shuf0, shuf1 and stopping at 9 h wall time.",
              file=f, flush=True)
    # ---- 3. sweep ----
    T_START = time.time(); LIMIT = 9 * 3600
    jobs_all = [(s, min(s + CHUNK, nd)) for s in range(0, nd, CHUNK)]
    stopped = False
    with mp.Pool(NPROC, initializer=A.init) as pool:
        for nm, ct in (("K4", CT4), ("shuf0", A.codes(sh[0])), ("shuf1", A.codes(sh[1]))):
            if stopped: break
            t1 = time.time(); HT = np.zeros(25, np.int64); HB = np.zeros(14, np.int64); RECS = []; done_alph = 0; ndone = 0
            nextp = 0.05
            for _, m, ht, hb, rec in pool.imap_unordered(job, [(nm, ct, s, e) for s, e in jobs_all]):
                HT += ht; HB += hb; RECS.extend(rec.tolist()); done_alph += m; ndone += 1
                frac = done_alph / nd
                if frac >= nextp - 1e-9 or ndone == len(jobs_all):
                    el = time.time() - t1
                    summary(f"A2 {nm} {100*frac:5.1f}%", int(HT.sum()), HT, HB, RECS, ct, f,
                            f"  elapsed {el/60:.0f} min, ETA {el*(1-frac)/max(frac,1e-9)/60:.0f} min")
                    while nextp <= frac: nextp += 0.05
                if time.time() - T_START > LIMIT:
                    stopped = True
                    print(f"STOP: 9 h wall limit reached during {nm} at {100*frac:.1f}% of alphabets", file=f, flush=True)
                    pool.terminate(); break
            line, ex = summary(f"A2 FINAL {nm}{' (PARTIAL)' if stopped else ''}", int(HT.sum()), HT, HB, RECS, ct, f,
                               f"  ({(time.time()-t1)/60:.0f} min)")
            print(line, flush=True)
            RECS.sort(key=lambda r: (-r[5], -r[6]))
            for r in RECS[:3]: print(f"      top: {describe(r)}", file=f, flush=True)
            for r in ex[:20]: print(f"      EXACT: {describe(r)} pt={row_pt(r, ct)}", file=f, flush=True)
    print(f"=== Task A2 done {time.ctime()} ===", file=f, flush=True)
