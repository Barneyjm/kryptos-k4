import re, glob, numpy as np, time, pickle, k4, gap1
from gromark_dict import straight, transposed
AZ = k4.AZ
def rev_remainder(w): k = "".join(dict.fromkeys(w)); return k + "".join(c for c in AZ[::-1] if c not in k)
def continue_from_last(w):
    k = "".join(dict.fromkeys(w)); start = AZ.index(k[-1]); rest = [AZ[(start + 1 + i) % 26] for i in range(26)]
    return k + "".join(c for c in rest if c not in k)
def build():
    al = {}
    for base, bn in ((AZ, "AZ"), (k4.KA, "KA")):
        for a in (1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25):
            for b in range(26): al.setdefault("".join(base[(a * i + b) % 26] for i in range(26)), f"decimate {bn} x{a}+{b}")
    en = sorted({re.sub("[^A-Z]", "", w.upper()) for w in open("/usr/share/dict/words")} - {""})
    t = " ".join(open(f, encoding="utf-8", errors="ignore").read() for f in glob.glob("de*.txt")).upper()
    for x, y in (("Ä", "AE"), ("Ö", "OE"), ("Ü", "UE"), ("ß", "SS")): t = t.replace(x, y)
    de = sorted({w for w in re.findall(r"[A-Z]+", t) if len(w) >= 3})
    for w in en:
        for f, a in (("revrem", rev_remainder(w)), ("cont", continue_from_last(w))):
            al.setdefault(a, f"{w}/{f}"); al.setdefault(a[::-1], f"{w}/{f}-rev")
    for w in de:
        for f, a in (("s", straight(w)), ("t", transposed(w)), ("revrem", rev_remainder(w)), ("cont", continue_from_last(w))):
            al.setdefault(a, f"DE:{w}/{f}"); al.setdefault(a[::-1], f"DE:{w}/{f}-rev")
    return list(al.items()), len(de)
if __name__ == "__main__":
    AL, nde = build(); print(f"new alphabets: {len(AL):,} (German words {nde:,})", flush=True)
    t = time.time(); h = gap1.scan(k4.CT, AL, "K4"); print(f"K4: {len(h)} hits in {time.time()-t:.0f}s", flush=True)
    for x in h[:20]: print("  ", x, flush=True)
    pickle.dump(h, open("gap2_hits.pkl", "wb"))
