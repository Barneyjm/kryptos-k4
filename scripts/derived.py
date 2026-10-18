"""Gap 6: keys derived from the solved panels by rule (initials, finals, after-X, word lengths, every n-th letter, misspelling
positions), tested as running key (exact/near, all offsets), as keyword alphabet + linear models, and as the periodic key."""
import re, numpy as np, itertools, k4, cokernel, gap1
K1 = "BETWEEN SUBTLE SHADING AND THE ABSENCE OF LIGHT LIES THE NUANCE OF IQLUSION"
K2 = "IT WAS TOTALLY INVISIBLE HOWS THAT POSSIBLE THEY USED THE EARTHS MAGNETIC FIELD X THE INFORMATION WAS GATHERED AND TRANSMITTED UNDERGRUUND TO AN UNKNOWN LOCATION X DOES LANGLEY KNOW ABOUT THIS THEY SHOULD ITS BURIED OUT THERE SOMEWHERE X WHO KNOWS THE EXACT LOCATION ONLY WW THIS WAS HIS LAST MESSAGE X THIRTY EIGHT DEGREES FIFTY SEVEN MINUTES SIX POINT FIVE SECONDS NORTH SEVENTY SEVEN DEGREES EIGHT MINUTES FORTY FOUR SECONDS WEST X LAYER TWO"
K3 = "SLOWLY DESPARATLY SLOWLY THE REMAINS OF PASSAGE DEBRIS THAT ENCUMBERED THE LOWER PART OF THE DOORWAY WAS REMOVED WITH TREMBLING HANDS I MADE A TINY BREACH IN THE UPPER LEFT HAND CORNER AND THEN WIDENING THE HOLE A LITTLE I INSERTED THE CANDLE AND PEERED IN THE HOT AIR ESCAPING FROM THE CHAMBER CAUSED THE FLAME TO FLICKER BUT PRESENTLY DETAILS OF THE ROOM WITHIN EMERGED FROM THE MIST X CAN YOU SEE ANYTHING Q"
ALL = K1 + " " + K2 + " " + K3
D = {}
for name, txt in (("K1", K1), ("K2", K2), ("K3", K3), ("ALL", ALL)):
    w = txt.split(); D[f"{name} initials"] = "".join(x[0] for x in w); D[f"{name} finals"] = "".join(x[-1] for x in w)
    D[f"{name} seconds"] = "".join(x[1] for x in w if len(x) > 1); D[f"{name} wordlen"] = "".join(chr(65 + (len(x) % 26)) for x in w)
    D[f"{name} after-X"] = "".join(w[i + 1][0] for i in range(len(w) - 1) if w[i] == "X")
    flat = txt.replace(" ", "")
    for n in range(2, 8): D[f"{name} every{n}"] = flat[::n]
D["misspell"] = "QUAE"; D["X-words"] = "".join(x for x in K2.split() if x == "X")
D["K2 numbers"] = "".join(chr(65 + int(d)) for d in "385765778844")
for k in list(D): D[k + " rev"] = D[k][::-1]
D = {k: v for k, v in D.items() if len(v) >= 4}
CR = list(zip(k4.CPOS.tolist(), k4.CPL)); POS = k4.CPOS
print("derived strings:", len(D))
# (a) running key / periodic key with std alphabets
hits = []
for name, s in D.items():
    L = len(s)
    for an, X in (("AZ", k4.AZ), ("KA", k4.KA)):
        kv = [X.index(c) for c in s]
        for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26), ("vb", lambda c, p: (p - c) % 26)):
            need = [f(X.index(k4.CT[i]), X.index(p)) for i, p in CR]
            for off in range(L):
                m = sum(kv[(i + off) % L] == v for (i, _), v in zip(CR, need))
                if m >= 14: hits.append((m, name, an, kn, off))
print("(a) derived string as (repeating/running) key, >= 14/24:", sorted(hits, reverse=True)[:5] or "none")
# (b) derived string as keyword alphabet, all linear models
keyed = lambda w: "".join(dict.fromkeys(w + k4.AZ))
AL = []
for name, s in D.items():
    a = keyed(s); AL += [(a, name + "/alph"), (a[::-1], name + "/alph-rev")]
P = cokernel.prepared(5); gap1.P = P; gap1.names = [p[0] for p in P]
gap1.B13 = np.vstack([p[2] for p in P]).astype(np.float32); gap1.B2 = np.vstack([p[1] for p in P]).astype(np.float32)
gap1.seg13 = np.cumsum([0] + [len(p[2]) for p in P])[:-1]; gap1.seg2 = np.cumsum([0] + [len(p[1]) for p in P])[:-1]
h = gap1.scan(k4.CT, AL, "K4"); print(f"(b) derived-keyword alphabets x {len(P)} linear models: {len(h)} hits", h[:5])
