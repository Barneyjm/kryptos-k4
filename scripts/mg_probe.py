import numpy as np, k4, mg_mask_grid as m, collections, sys
POS = m.POS
def left_null(A, p):
    # rows r with r A = 0 mod p  -> null space of A^T
    ok, x, B, rk = m.solve(A.T.copy(), np.zeros(A.shape[1], np.int64), p); return B
def canon(fam, kf):
    t = [kf(i) for i in POS]; comps = collections.defaultdict(dict); out = []
    for tl in t:
        row = []
        for lab, co in tl:
            d = comps[lab[0]]; d.setdefault(lab, len(d)); row.append((lab[0], d[lab], co))
        out.append(tuple(row))
    return tuple(out)
if __name__ == "__main__":
    want = sys.argv[1:]  # list of model names
    for fam, name, kf in m.models():
        for side in ("plain", "cipher"):
            for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
                key = f"{fam}|{name}|{side}|{an}"
                if key not in want: continue
                A, b, nx, letters = m.build(k4.CT, side, X, kf)
                for p in (2, 13):
                    N_ = left_null(A, p)
                    for v in N_:
                        nzr = [(POS[r] if r < 24 else 'g', int(v[r])) for r in range(25) if v[r]]
                        print(key, p, "check rows:", nzr, " value", int(v @ b % p))
