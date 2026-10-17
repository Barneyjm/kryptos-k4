import numpy as np, k4
rng = np.random.default_rng(42)
A = np.array([ord(c) - 65 for c in k4.CT]); N = 97
def vrep(a, w):          # repeated vertical bigrams at width w: (#bigrams) - (#distinct)
    b = a[:-w] * 26 + a[w:]; return len(b) - len(np.unique(b))
W = range(2, 49)
obs = {w: vrep(A, w) for w in W}
print("width 21 repeated vertical bigrams:", obs[21])
sims = np.array([[vrep(p, w) for w in W] for p in (rng.permutation(A) for _ in range(40000))])
pw = {w: (sims[:, j] >= obs[w]).mean() for j, w in enumerate(W)}
print("p at width 21 (Bean ~1/6750): %.6f  = 1 in %.0f" % (pw[21], 1 / max(pw[21], 1e-9)))
print("most extreme widths:", sorted(pw.items(), key=lambda x: x[1])[:5])
minp_null = np.array([min((sims[:, j] >= sims[t, j]).mean() for j in range(len(W))) for t in range(3000)])
print("corrected for choosing the best of widths 2-48: P(min p <= %.6f) = %.4f" % (pw[21], (minp_null <= pw[21]).mean()))
# minor differences
def cdist(a, b): d = abs(a - b) % 26; return min(d, 26 - d)
PT = {i: ord(c) - 65 for i, c in zip(k4.CPOS.tolist(), k4.CPL)}
sub = [i for i in PT if k4.CPL[list(PT).index(i)] in "KRYPTOS"]
s_sub = sum(cdist(PT[i], A[i]) for i in sub); s_all = sum(cdist(PT[i], A[i]) for i in PT)
null_sub, null_all = [], []
for _ in range(100000):
    p = rng.permutation(A); null_sub.append(sum(cdist(PT[i], p[i]) for i in sub)); null_all.append(sum(cdist(PT[i], p[i]) for i in PT))
print(f"\nminor differences, letters in KRYPTOS ({len(sub)} letters): sum {s_sub}; P = {np.mean(np.array(null_sub) <= s_sub):.5f} (Bean ~1/5520)")
print(f"minor differences, ALL 24 known letters: sum {s_all} (random mean {np.mean(null_all):.1f}); P = {np.mean(np.array(null_all) <= s_all):.4f}")
