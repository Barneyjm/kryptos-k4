# plant a solution and make sure the pruned search finds it
import random, k4, prune
random.seed(9)
P = list("".join(random.choice(k4.AZ) for _ in range(97)))
for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
X = k4.KA; key = "PALIMPSEST"
sub = lambda t: "".join(X[(X.index(c) + X.index(key[i % len(key)])) % 26] for i, c in enumerate(t))
for w, col in ((7, (3, 0, 6, 2, 5, 1, 4)), (12, (5, 11, 0, 7, 2, 9, 4, 1, 10, 3, 8, 6))):
    order = k4.order_from(w, col)
    for model in "AB":
        C = "".join(sub("".join(P))[o] for o in order) if model == "A" else sub("".join(P[o] for o in order))
        prune.CT = C
        s = prune.Search(w, "enc"); s.run(col[:2])
        print(w, model, "nodes", s.nodes, "found planted:", [h for h in s.hits if h[2] == col and h[5] == 10])
