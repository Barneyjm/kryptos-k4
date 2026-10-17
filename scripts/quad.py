import numpy as np, re, glob, os
def build():
    if os.path.exists("quad.npy"): return np.load("quad.npy")
    t = re.sub("[^A-Z]", "", " ".join(open(f, encoding="utf-8", errors="ignore").read().upper() for f in glob.glob("g*.txt")))
    a = np.frombuffer(t.encode(), dtype=np.uint8).astype(np.int64) - 65
    idx = ((a[:-3] * 26 + a[1:-2]) * 26 + a[2:-1]) * 26 + a[3:]
    c = np.bincount(idx, minlength=26**4).astype(float) + 0.01
    q = np.log10(c / c.sum()); np.save("quad.npy", q); return q
Q = build()
def score(s):
    a = np.array([ord(ch) - 65 for ch in s]); i = ((a[:-3]*26 + a[1:-2])*26 + a[2:-1])*26 + a[3:]
    return Q[i].mean()
if __name__ == "__main__":
    import random
    print("english sample", score("THEREMAINSOFPASSAGE"), score("SHADOWFORCESLUCID"))
    rs = [score("".join(random.choice("ABCDEFGHIJKLMNOPQRSTUVWXYZ") for _ in range(13))) for _ in range(2000)]
    print("random mean %.2f sd %.2f max %.2f" % (np.mean(rs), np.std(rs), max(rs)))
    eng = open("g98.txt").read().upper(); eng = re.sub("[^A-Z]", "", eng)
    es = [score(eng[i:i+13]) for i in range(10000, 400000, 197)]
    print("english 13-letter mean %.2f sd %.2f min %.2f" % (np.mean(es), np.std(es), min(es)))
