"""Shared helpers for the fe_* testers: alphabets, planted plaintexts, and a fast
'equal-key partition' engine.

Partition engine: a key model assigns each of the 24 crib positions a class label
(positions in the same class must carry the same key value). A model is reduced to
a spanning set of edges (pairs of crib indices, 0..23). For a batch of implied-key
vectors k (n x 24), we compute which of the 276 pairs are UNEQUAL, then the number
of violated edges per model with one float32 matmul. A model fits iff 0 violations.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "3"); os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "3")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "3")
import re, itertools, numpy as np, k4
from texts import KEYWORDS

AZ, KA, CT = k4.AZ, k4.KA, k4.CT
POS = [int(p) for p in k4.CPOS]               # 0-indexed crib positions
PL = np.array([ord(c) - 65 for c in k4.CPL])
PAIRS = list(itertools.combinations(range(24), 2))
PIDX = {p: j for j, p in enumerate(PAIRS)}
PA = np.array([a for a, b in PAIRS]); PB = np.array([b for a, b in PAIRS])

def keyed(w):
    return "".join(dict.fromkeys(w + AZ))

def clue_alphabets():
    al = {AZ: "AZ", KA: "KA"}
    al.setdefault(AZ[::-1], "AZrev"); al.setdefault(KA[::-1], "KArev")
    for w in KEYWORDS:
        al.setdefault(keyed(w), "kw:" + w); al.setdefault(keyed(w)[::-1], "kwrev:" + w)
    return list(al.items())

def all_alphabets():
    import rk_mixed
    d = {a: n for a, n in rk_mixed.alph_list()}
    for a, n in clue_alphabets():
        d.setdefault(a, n)
    return list(d.items())

_ENG = None
def english(n, rng):
    global _ENG
    if _ENG is None:
        t = open("g98.txt", encoding="utf-8", errors="ignore").read().upper()
        _ENG = re.sub("[^A-Z]", "", t)[20000:]
    s = int(rng.integers(0, len(_ENG) - n)); return _ENG[s:s + n]

def planted_plaintext(rng):
    p = list(english(97, rng))
    for i, ch in zip(POS, k4.CPL): p[i] = ch
    return "".join(p)

def model_from_labels(labels):
    """labels: list of 24 hashables (None = unconstrained singleton). Returns (edge tuple, checks)."""
    cls = {}
    for j, l in enumerate(labels):
        if l is None: continue
        cls.setdefault(l, []).append(j)
    edges = []
    for mem in cls.values():
        for a, b in zip(mem, mem[1:]): edges.append(PIDX[(a, b)])
    return tuple(sorted(edges))

class Engine:
    def __init__(self, named_models, min_checks):
        """named_models: iterable of (name, labels). Dedupe identical edge sets; keep >= min_checks."""
        seen = {}
        for name, labels in named_models:
            e = model_from_labels(labels)
            if len(e) < min_checks: continue
            if e in seen: seen[e][1] += 1
            else: seen[e] = [name, 1]
        self.edges = list(seen.keys()); self.names = [v[0] for v in seen.values()]
        self.alias = [v[1] for v in seen.values()]
        self.checks = np.array([len(e) for e in self.edges])
        M = np.zeros((276, len(self.edges)), np.float32)
        for m, e in enumerate(self.edges): M[list(e), m] = 1
        self.M = M
    def __len__(self): return len(self.edges)
    def fits(self, K, chunk=40000):
        """K: (n,24) int implied keys. Returns list of (row, model) with zero violations."""
        out = []
        for s in range(0, len(K), chunk):
            k = K[s:s + chunk]
            neq = (k[:, PA] != k[:, PB]).astype(np.float32)
            v = neq @ self.M
            r, m = np.nonzero(v == 0)
            out.extend(zip((r + s).tolist(), m.tolist()))
        return out

def inv_index(alphs):
    """alphs: (n,26) letter indices -> (n,26) position of each letter."""
    return np.argsort(alphs, axis=1)

def alph_array(strings):
    return np.array([[ord(c) - 65 for c in a] for a in strings], dtype=np.int16)
