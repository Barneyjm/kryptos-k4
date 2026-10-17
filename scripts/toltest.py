import random, numpy as np, k4, tol
random.seed(4)
P = list("".join(random.choice(k4.AZ) for _ in range(97)))
for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
X = k4.KA; key = "ABSCISSA"
sub = lambda t: "".join(X[(X.index(c) + X.index(key[i % len(key)])) % 26] for i, c in enumerate(t))
w, col = 6, (2, 5, 0, 3, 1, 4); order = k4.order_from(w, col)
C = list("".join(sub("".join(P))[o] for o in order))
dest = np.empty(97, int); dest[order] = np.arange(97)
for e in (1, 2):
    Cx = C[:]
    for pos in random.sample(list(k4.CPOS), e):        # corrupt the carved letter of e crib positions
        j = dest[pos]; Cx[j] = random.choice([c for c in k4.AZ if c != Cx[j]])
    k4.CT = "".join(Cx); hits = []
    tol.tol_check(dest, k4.CT, k4.CPOS, k4.CPL, lambda i: "planted", hits, tol.STD)
    print(f"{e} slip(s):", [h for h in hits if h[4] == 8])
