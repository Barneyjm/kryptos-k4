"""Given a Gromark key K (97 digits), find plain/cipher alphabets maximizing English quadgram score of the decryption,
with the 24 crib letters forced. Model: cpos(C_i) = ppos(P_i) + K_i (mod 26). Unknowns: ppos, cpos permutations.
Decrypt: P_i = pinv[(cpos[C_i] - K_i) % 26]. Simulated annealing over swaps in either alphabet; cribs as heavy penalty."""
import numpy as np, random, math, re, k4
from quad import Q
CRIB = dict(zip(k4.CPOS.tolist(), k4.CPL))
CRIB_IDX = np.array(list(CRIB)); CRIB_VAL = np.array([ord(c) - 65 for c in CRIB.values()])
def decrypt(ct, K, ppos, cpos):
    pinv = np.argsort(ppos); return pinv[(cpos[ct] - K) % 26]
def fitness(pt):
    q = ((pt[:-3] * 26 + pt[1:-2]) * 26 + pt[2:-1]) * 26 + pt[3:]
    return Q[q].sum() - 40.0 * (pt[CRIB_IDX] != CRIB_VAL).sum()
def solve(ct_str, K, iters=60000, restarts=6, seed=0):
    rnd = random.Random(seed); ct = np.array([ord(c) - 65 for c in ct_str]); K = np.asarray(K)
    best = (-1e9, None)
    for r in range(restarts):
        ppos = np.random.default_rng(seed + r).permutation(26); cpos = np.random.default_rng(seed + 100 + r).permutation(26)
        cur = fitness(decrypt(ct, K, ppos, cpos)); T = 8.0
        for it in range(iters):
            which = cpos if rnd.random() < 0.5 else ppos
            a, b = rnd.randrange(26), rnd.randrange(26)
            which[a], which[b] = which[b], which[a]
            f = fitness(decrypt(ct, K, ppos, cpos))
            if f >= cur or rnd.random() < math.exp((f - cur) / T): cur = f
            else: which[a], which[b] = which[b], which[a]
            T = max(0.2, T * 0.9998)
        if cur > best[0]: best = (cur, "".join(chr(65 + x) for x in decrypt(ct, K, ppos, cpos)), ppos.copy(), cpos.copy())
    return best
if __name__ == "__main__":
    import gromark2 as G, time
    # plant: English plaintext with cribs, primer 31415 base 10 std rule, random alphabets
    eng = re.sub("[^A-Z]", "", open("g98.txt").read().upper())[200000:200097]
    P = list(eng)
    for s, w in k4.CRIBS.items(): P[s:s+len(w)] = w
    K = G.expand(np.array([[3, 1, 4, 1, 5]]), 10, (4, 5), 97)[0]
    rng = np.random.default_rng(7); pp, cp = rng.permutation(26), rng.permutation(26)
    ct = "".join(chr(65 + int(np.argsort(cp)[(pp[ord(c) - 65] + K[i]) % 26])) for i, c in enumerate(P))
    t = time.time(); sc, txt, _, _ = solve(ct, K, iters=40000, restarts=4)
    print(f"planted truth: {''.join(P)}\nrecovered:     {txt}\nscore {sc:.1f}  ({time.time()-t:.0f}s)")
    agree = sum(a == b for a, b in zip(txt, P)); print("letters correct:", agree, "/ 97")
