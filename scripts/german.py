import numpy as np, re, glob, os, random, itertools, multiprocessing as mp, k4, rkscan, dumb
def de_text():
    t = " ".join(open(f, encoding="utf-8", errors="ignore").read() for f in glob.glob("de*.txt")).upper()
    for a, b in (("Ä", "AE"), ("Ö", "OE"), ("Ü", "UE"), ("ß", "SS")): t = t.replace(a, b)
    return re.sub("[^A-Z]", "", t)
DE = de_text()
a = np.frombuffer(DE.encode(), dtype=np.uint8).astype(np.int64) - 65
idx = ((a[:-3]*26 + a[1:-2])*26 + a[2:-1])*26 + a[3:]
c = np.bincount(idx, minlength=26**4).astype(float) + 0.01; QDE = np.log10(c / c.sum())
rkscan.Q = QDE                                   # rescore with German stats
def sc(s):
    x = np.array([ord(ch)-65 for ch in s]); i = ((x[:-3]*26 + x[1:-2])*26 + x[2:-1])*26 + x[3:]; return QDE[i].mean()
if __name__ == "__main__":
    random.seed(2)
    g = [(sc(DE[i:i+13]) + sc(DE[i+500:i+511])) / 2 for i in range(5000, 300000, 1500)]
    r = [(sc("".join(random.choice(k4.AZ) for _ in range(13))) + sc("".join(random.choice(k4.AZ) for _ in range(11)))) / 2 for _ in range(2000)]
    print("calibration: German fragments mean %.2f sd %.2f | random mean %.2f sd %.2f max %.2f" % (np.mean(g), np.std(g), np.mean(r), np.std(r), max(r)))
    keep = []
    rkscan.scan(np.arange(97)[None, :], lambda i: "no transposition", keep)
    lib = dumb.library(); names = list(lib); orders = np.array([lib[k] for k in names])
    dests = np.empty_like(orders); dests[np.arange(len(names))[:, None], orders] = np.arange(97)
    rkscan.scan(dests, lambda i: names[i] + " enc", keep); rkscan.scan(orders, lambda i: names[i] + " dec", keep)
    with mp.Pool(4) as pool: res = pool.map(rkscan.worker, [(w, f) for w in range(2, 10) for f in range(w)])
    keep += [x for r in res for x in r]
    keep.sort(reverse=True)
    print("best German-key candidates (German ~ above, noise ~ random):")
    for x in keep[:12]: print("  %.2f" % x[0], x[1:])
