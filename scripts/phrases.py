"""Phrase-derived keys from Sanborn's own words (2009 oral history), K1-K3 plaintext, Morse, clue quotes."""
import re, numpy as np, itertools, multiprocessing as mp, k4, texts
from gromark_dict import straight, transposed
import rk_mixed
CLUES = """BERLIN CLOCK EAST NORTHEAST DELVE INTO THAT PARTICULAR CLOCK THERE ARE SEVERAL REALLY INTERESTING CLOCKS IN BERLIN
WORLD CLOCK DELIVERING A MESSAGE EGYPT NINETEEN EIGHTY SIX FALL OF THE BERLIN WALL NINETEEN EIGHTY NINE
MASKING TECHNIQUE CHANGE THE LANGUAGE BASE MORE THAN ONE STAGE WHAT WOULD BE PROVIDED AGENTS OR PILOTS IN CASE OF CAPTURE
IT IS SIMPLE CAN BE REMEMBERED AND EXECUTED YEARS LATER WHEN USED WITH THE CORRECT KEY WORDS
VIRTUALLY INVISIBLE DIGETAL INTERPRETATIU SHADOW FORCES LUCID MEMORY T IS YOUR POSITION SOS RQ
WONDERFUL THINGS YES WONDERFUL THINGS CAN YOU SEE ANYTHING"""
def phrase_list():
    srcs = [open("sanborn_oral_history.txt").read(), CLUES,
            "BETWEEN SUBTLE SHADING AND THE ABSENCE OF LIGHT LIES THE NUANCE OF IQLUSION",
            "IT WAS TOTALLY INVISIBLE HOWS THAT POSSIBLE THEY USED THE EARTHS MAGNETIC FIELD X THE INFORMATION WAS GATHERED AND TRANSMITTED UNDERGRUUND TO AN UNKNOWN LOCATION X DOES LANGLEY KNOW ABOUT THIS THEY SHOULD ITS BURIED OUT THERE SOMEWHERE X WHO KNOWS THE EXACT LOCATION ONLY WW THIS WAS HIS LAST MESSAGE X THIRTY EIGHT DEGREES FIFTY SEVEN MINUTES SIX POINT FIVE SECONDS NORTH SEVENTY SEVEN DEGREES EIGHT MINUTES FORTY FOUR SECONDS WEST X LAYER TWO",
            "SLOWLY DESPARATLY SLOWLY THE REMAINS OF PASSAGE DEBRIS THAT ENCUMBERED THE LOWER PART OF THE DOORWAY WAS REMOVED WITH TREMBLING HANDS I MADE A TINY BREACH IN THE UPPER LEFT HAND CORNER AND THEN WIDENING THE HOLE A LITTLE I INSERTED THE CANDLE AND PEERED IN THE HOT AIR ESCAPING FROM THE CHAMBER CAUSED THE FLAME TO FLICKER BUT PRESENTLY DETAILS OF THE ROOM WITHIN EMERGED FROM THE MIST X CAN YOU SEE ANYTHING Q"]
    out = set()
    for s in srcs:
        w = re.findall(r"[A-Z]+", s.upper())
        for n in range(1, 6):
            for i in range(len(w) - n + 1):
                p = "".join(w[i:i + n])
                if 4 <= len(p) <= 40: out.add(p)
    return sorted(out)
CR = list(zip(k4.CPOS.tolist(), k4.CPL)); POS = k4.CPOS
PAIRS = {L: [(a, b) for a, b in itertools.combinations(range(24), 2) if POS[a] % L == POS[b] % L] for L in range(1, 27)}
def periodic_and_digits(args):
    names, alphs = args
    A = np.array([[ord(c) - 65 for c in a] for a in alphs]); INV = np.argsort(A, axis=1); n = len(A)
    C = np.array([ord(c) - 65 for c in k4.CT])[POS]; P = np.array([ord(c) - 65 for c in k4.CPL])
    AZ = np.tile(np.arange(26), (n, 1)); hits = []; small = []
    for pn, Pi in (("AZ", AZ), ("A", INV)):
        for cn, Ci in (("AZ", AZ), ("A", INV)):
            if pn == cn == "AZ": continue
            pc = np.take_along_axis(Pi, np.tile(P, (n, 1)), 1); cc = np.take_along_axis(Ci, np.tile(C, (n, 1)), 1)
            for kn, k in (("vig", (cc - pc) % 26), ("beau", (cc + pc) % 26)):
                for L, pr in PAIRS.items():
                    if 24 - len({POS[i] % L for i in range(24)}) < 8: continue
                    ok = np.all(k[:, [a for a, _ in pr]] == k[:, [b for _, b in pr]], axis=1)
                    for j in np.nonzero(ok)[0]: hits.append((names[j], pn, cn, kn, L))
                mx = k.max(axis=1)
                for j in np.nonzero(mx <= 12)[0]: small.append((int(mx[j]), names[j], pn, cn, kn))
    return hits, small
def phrase_as_key(args):
    phrases = args; out = []
    for ph in phrases:
        L = len(ph)
        for an, X in (("AZ", k4.AZ), ("KA", k4.KA), ("own", straight(ph))):
            kv = [X.index(c) for c in ph]
            for kn, f in (("vig", lambda c, p: (c - p) % 26), ("beau", lambda c, p: (c + p) % 26), ("vb", lambda c, p: (p - c) % 26)):
                need = [f(X.index(k4.CT[i]), X.index(p)) for i, p in CR]
                for off in range(L):
                    m = sum(kv[(i + off) % L] == v for (i, _), v in zip(CR, need))
                    if m >= 15: out.append((m, ph, an, kn, off))
    return out
if __name__ == "__main__":
    PH = phrase_list(); print(f"phrases: {len(PH):,}", flush=True)
    al = {}
    for p in PH:
        for f, a in (("s", straight(p)), ("t", transposed(p))):
            al.setdefault(a, f"{p}/{f}"); al.setdefault(a[::-1], f"{p}/{f}-rev")
    items = list(al.items()); print(f"distinct phrase alphabets: {len(items):,}", flush=True)
    chunks = [items[i:i + 5000] for i in range(0, len(items), 5000)]
    with mp.Pool(11) as pool:
        res = pool.map(periodic_and_digits, [([n for _, n in ch], [a for a, _ in ch]) for ch in chunks])
        hits = [h for r in res for h in r[0]]; small = sorted([s for r in res for s in r[1]])
        print("(a) periodic key (L<=26) with phrase alphabets, exact:", hits[:10] or "none", flush=True)
        print("(c) smallest max key value with phrase alphabets:", small[:5] or "none <= 12", flush=True)
        rk = sorted([x for r in pool.map(rk_mixed.job, [("en", k4.CT, ch) for ch in chunks]) for x in r], reverse=True)
        rng = np.random.default_rng(0); nulls = []
        for t in range(3):
            ct = "".join(rng.permutation(list(k4.CT)))
            nulls.append(max(x[0] for r in pool.map(rk_mixed.job, [("en", ct, ch) for ch in chunks]) for x in r))
        print(f"(b) running-key English score with phrase alphabets: K4 best {rk[0][0]:.2f} {rk[0][1:]} | shuffled nulls {[round(x, 2) for x in nulls]}", flush=True)
        long_ph = [p for p in PH if len(p) >= 5]
        pk = sorted([x for r in pool.map(phrase_as_key, [long_ph[i:i + 3000] for i in range(0, len(long_ph), 3000)]) for x in r], reverse=True)
        print("(d) phrase itself as the key (full length, every phase), >= 15/24:", pk[:8] or "none", flush=True)
