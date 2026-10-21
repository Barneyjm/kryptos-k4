"""Shared loader for the 613,440 Sanborn-phrase alphabets (same construction as gap_phr.py)."""
import hashlib, os, pickle, sys

RES = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(RES)
SCR = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCR)
os.chdir(SCR)

import numpy as np


def phrase_alphabets():
    """Return (A uint8 (n,26), names). Deduped straight/transposed/reversed, gap_phr order."""
    cache = os.path.join(SCR, "phrase_alph.pkl")
    if os.path.exists(cache):
        d = pickle.load(open(cache, "rb"))
        return d["A"], d["names"]
    import phrases
    from gromark_dict import straight, transposed
    al = {}
    for p in phrases.phrase_list():
        for f, a in (("s", straight(p)), ("t", transposed(p))):
            al.setdefault(a, f"{p}/{f}")
            al.setdefault(a[::-1], f"{p}/{f}-rev")
    items = list(al.items())
    A = np.array([[ord(c) - 65 for c in a] for a, _ in items], np.uint8)
    names = [n for _, n in items]
    pickle.dump({"A": A, "names": names}, open(cache, "wb"), protocol=4)
    return A, names


def digest(A):
    return hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest()[:16]


def strings(A):
    return ["".join(chr(65 + int(c)) for c in row) for row in A]
