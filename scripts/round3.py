import numpy as np, itertools, k4
from texts import KEYWORDS
CT, CPOS, CPL = k4.CT, list(k4.CPOS), k4.CPL
P = dict(zip(CPOS, CPL)); AZ = k4.AZ
def keyed(word):
    s = "".join(dict.fromkeys(word + AZ)); return s
ALPHA = {"AZ": AZ, "KA": k4.KA}
for w in KEYWORDS: ALPHA["kw:" + w] = keyed(w)
for w in KEYWORDS: ALPHA["kwrev:" + w] = keyed(w)[::-1]
hits = []
# 1. Quagmire-style: plaintext alphabet PA, cipher alphabet CA, periodic shift
for (pn, PA), (cn, CA) in itertools.product(ALPHA.items(), repeat=2):
    for sign in (1, -1):
        k = [(CA.index(CT[i]) - sign * PA.index(P[i])) % 26 for i in CPOS]
        for L in range(1, 27):
            d = {}; ok = True
            for i, kv in zip(CPOS, k):
                if d.setdefault(i % L, kv) != kv: ok = False; break
            if ok and len(CPOS) - len(d) >= 8: hits.append(("quag", pn, cn, sign, L))
print("quagmire hits", len(hits), hits[:10])
# 2. Ciphertext autokey (key[i] = C[i-L]) and plaintext autokey (key[i]=P[i-L])
ak = []
for an, X in ALPHA.items():
    ix = X.index
    for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26), ("vb", lambda c, p: (p - c) % 26)):
        for L in range(1, 97):
            ck = [f(ix(CT[i]), ix(P[i])) == ix(CT[i - L]) for i in CPOS if i - L >= 0]
            if ck and all(ck): ak.append(("ctak", an, kn, L))
            pk = [f(ix(CT[i]), ix(P[i])) == ix(P[i - L]) for i in CPOS if i - L in P]
            if len(pk) >= 4 and all(pk): ak.append(("ptak", an, kn, L, len(pk)))
print("autokey hits", ak[:20])
# 3. Progressive key: k_i = base[i % L] + d * (i // L)
pg = []
for an, X in ALPHA.items():
    for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26)):
        k = {i: f(X.index(CT[i]), X.index(P[i])) for i in CPOS}
        for L in range(1, 27):
            for d in range(1, 26):
                base = {}; ok = True
                for i in CPOS:
                    b = (k[i] - d * (i // L)) % 26
                    if base.setdefault(i % L, b) != b: ok = False; break
                if ok and len(CPOS) - len(base) >= 8: pg.append((an, kn, L, d))
print("progressive hits", pg[:20])
# 4. Hill cipher 2x2 and 3x3 (mod 26), both block alignments, AZ and KA numbering
def inv_mod(M):
    from sympy import Matrix
    try: return np.array(Matrix(M).inv_mod(26)).astype(int)
    except Exception: return None
hill = []
for an, X in (("AZ", AZ), ("KA", k4.KA)):
    for n in (2, 3):
        for align in range(n):
            blocks = []
            for s in range(align, 97, n):
                if all(j in P for j in range(s, s + n)):
                    blocks.append((np.array([X.index(P[j]) for j in range(s, s + n)]), np.array([X.index(CT[j]) for j in range(s, s + n)])))
            if len(blocks) <= n: continue
            for combo in itertools.combinations(range(len(blocks)), n):
                Pm = np.array([blocks[c][0] for c in combo]).T; Cm = np.array([blocks[c][1] for c in combo]).T
                Pi = inv_mod(Pm.tolist())
                if Pi is None: continue
                K = (Cm @ Pi) % 26
                if all(((K @ p) % 26 == c).all() for p, c in blocks):
                    hill.append((an, n, align, K.tolist())); break
print("hill hits", hill)
