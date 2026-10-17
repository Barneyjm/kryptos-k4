import random, k4, rowcol
random.seed(5); P = [random.choice(k4.AZ) for _ in range(97)]
for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
X = k4.KA; col = [3, 17, 8, 0, 21, 5, 12]; row = [random.randrange(26) for _ in range(20)]; off = 2
ct = "".join(X[(X.index(p) + col[(i + off) % 7] + row[(i + off) // 7]) % 26] for i, p in enumerate(P))
k4.CT = ct; rowcol.k4.CT = ct
print("planted W=7 off=2 KA:", [h for h in rowcol.test(("KA", "KA", [7])) if h[4] == 7 and h[5] == 2])
