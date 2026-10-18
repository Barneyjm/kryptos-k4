"""Side check (idea 4): ANY periodic polyalphabetic cipher with arbitrary independent alphabets per residue class
(strip cipher M-138-A, DRYAD-style sheets, unknown tableaux): within each class i = r (mod L), plain->cipher must be a
partial bijection at the crib positions. Also the same with an unknown phase slip between the cribs.
Exact, alphabet-free. Validated on a planted period-L cipher with random alphabets; shuffled CTs for calibration."""
import numpy as np, k4, fe_common as F
def consistent(ct, L, slip=0):
    f, g = {}, {}
    for p, ch in zip(F.POS, k4.CPL):
        r = ((p + (slip if p >= 63 else 0)) % L); c = ct[p]
        if f.setdefault((r, ch), c) != c or g.setdefault((r, c), ch) != ch: return False
    return True
def checks(L):
    return 24 - len({p % L for p in F.POS})
rng = np.random.default_rng(2)
for L in (7, 11, 13, 19):
    al = [rng.permutation(26) for _ in range(L)]; pt = F.planted_plaintext(rng)
    ct = "".join(chr(65 + al[i % L][ord(c) - 65]) for i, c in enumerate(pt))
    assert consistent(ct, L), L
print("planted independent-alphabet ciphers L=7,11,13,19 recovered: True")
Ls = range(1, 49)
k = [L for L in Ls if consistent(k4.CT, L)]
print("K4 consistent periods (no slip), L=1..48:", k)
print("   class-collision count (24 - #classes) for those:", [checks(L) for L in k])
ks = sorted({L for L in range(1, 27) for s in range(L) if consistent(k4.CT, L, s)})
print("K4 consistent with some phase slip between cribs, L=1..26:", ks)
sh = []
for t in range(200):
    c = "".join(rng.permutation(list(k4.CT))); sh.append([L for L in range(1, 27) if consistent(c, L)])
import collections
cnt = collections.Counter(L for s in sh for L in s)
print("shuffled CT (200): fraction consistent per L<=26:", {L: round(cnt[L] / 200, 2) for L in range(1, 27)})
print("K4 smallest consistent L:", min(k), " shuffled median smallest L:", int(np.median([min(s) if s else 99 for s in sh])))
