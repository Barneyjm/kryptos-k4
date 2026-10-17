import numpy as np, itertools, k4, gromark_free as G
AZ, KA = k4.AZ, k4.KA
def as_digits(word, base, alph):
    return [alph.index(c) % base for c in word]
words = ["DYAHR", "YAR", "RAY", "DYAR", "YAHR", "HRDYA", "RHAYD", "DIGETAL", "KRYPTOS", "BERLIN", "CLOCK", "EAST", "NORTHEAST"]
nums = ["1986", "1989", "1990", "19861989", "11091989", "198619891990", "6469", "7074", "2634", "2225", "385765778", "3857657784",
        "38576577844", "22252634", "2225263464697074", "2134", "6373", "2133", "6374", "1969", "0930", "24", "148", "146"]
cands = []
for w in words:
    for an, A in (("AZ", AZ), ("KA", KA)):
        for base in (5, 8, 10, 12, 24, 26):
            cands.append((f"{w}/{an}/b{base}", base, as_digits(w, base, A)))
for s in nums:
    for base in (10, 12, 24, 26):
        if all(int(ch) < base for ch in s): cands.append((f"{s}/b{base}", base, [int(ch) for ch in s]))
print("candidate primers:", len(cands))
hits = []
for name, b, prim in cands:
    n = len(prim)
    lagsets = {(n - 1, n), (1, n), (1, 2)} if n > 2 else {(1, 2)}
    for lags in lagsets:
        for pa, PA in (("AZ", AZ), ("KA", KA)):
            for sign in (1, -1):
                d = G.expand(np.array([prim]), b, n, lags)[:, G.POS]
                pidx = np.array([PA.index(c) for c in G.PL])
                slot = (pidx[None, :] + sign * d) % 26
                ok = np.all(slot[:, G.SA] == slot[:, G.SB]) and np.all(slot[:, G.DA] != slot[:, G.DB])
                if ok: hits.append((name, lags, pa, sign))
print("primers consistent with a permutation cipher alphabet:", hits or "none")
# also: base-26 'Vimark' with the cipher alphabet fixed to AZ or KA (no freedom): direct check
