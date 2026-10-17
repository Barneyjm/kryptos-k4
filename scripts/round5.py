import itertools, k4, mpmath
from texts import KEYWORDS
CT, P = k4.CT, dict(zip(k4.CPOS.tolist(), k4.CPL))
ALPHS = {"AZ": k4.AZ, "KA": k4.KA}
F = {"vig": lambda c, p: (c - p) % 26, "beau": lambda c, p: (c + p) % 26, "vb": lambda c, p: (p - c) % 26}
need = {(an, fn): {i: f(X.index(CT[i]), X.index(p)) for i, p in P.items()} for an, X in ALPHS.items() for fn, f in F.items()}
hits = []
# polynomial keys k_i = a i^2 + b i + c (mod 26), also mod 24 / mod 12 "clock" keys reduced into 26
for (an, fn), k in need.items():
    for a, b, c in itertools.product(range(26), repeat=3):
        if all((a*i*i + b*i + c) % 26 == v for i, v in k.items()): hits.append(("poly26", an, fn, a, b, c))
    for m in (12, 24, 60):
        for a, b, c in itertools.product(range(m), range(m), range(m)):
            if all((a*i*i + b*i + c) % m == v for i, v in k.items()): hits.append((f"poly{m}", an, fn, a, b, c))
print("polynomial hits", hits)
# digit streams (Gronsfeld, non-periodic): pi, e, sqrt2, phi, K2 coords; also pairs of digits as 00-25
mpmath.mp.dps = 2100
streams = {n: mpmath.nstr(v, 2000, strip_zeros=False).replace(".", "") for n, v in
           (("pi", mpmath.pi), ("e", mpmath.e), ("sqrt2", mpmath.sqrt(2)), ("phi", mpmath.phi), ("sqrt3", mpmath.sqrt(3)))}
streams["coords"] = "3857657784438576577844" * 10
streams["1986_1989"] = "19861989" * 30
dh = []
for (an, fn), k in need.items():
    for sn, s in streams.items():
        d = [int(x) for x in s]
        for off in range(len(d) - 97):
            if all(d[off + i] == v for i, v in k.items()): dh.append((sn, an, fn, off))
        d2 = [int(s[j:j+2]) % 26 for j in range(0, len(s) - 1, 2)]
        for off in range(len(d2) - 97):
            if all(d2[off + i] == v for i, v in k.items()): dh.append((sn + "-pairs", an, fn, off))
print("digit-stream hits", dh)
