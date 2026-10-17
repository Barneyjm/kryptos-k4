import random, numpy as np, k4
random.seed(1)
P = list("".join(random.choice(k4.AZ) for _ in range(97)))
for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
X = k4.KA; key = "PALIMPSEST"
def sub(t):  # vigenere over KA, periodic key
    return "".join(X[(X.index(c) + X.index(key[i % len(key)])) % 26] for i, c in enumerate(t))
w, col = 7, (3, 0, 6, 2, 5, 1, 4)
order = k4.order_from(w, col)
for model in "AB":
    if model == "A": C = "".join(sub("".join(P))[o] for o in order)
    else: C = sub("".join(P[o] for o in order))
    k4.CT = C
    hits = k4.run_columnar(7)
    print(model, [h for h in hits if str(col) in h[0] and h[4] == 10][:3], "total", len(hits))
