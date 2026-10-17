import itertools, numpy as np, k4, round2, tol
from texts import KEYWORDS
DEKW = ["WELTZEITUHR","URANIAWELTZEITUHR","URANIA","ALEXANDERPLATZ","ALEX","MAUER","BERLINERMAUER","MAUERFALL","DIEMAUER",
 "UHR","ZEIT","WELTZEIT","WELTUHR","OSTBERLIN","WESTBERLIN","DDR","CHECKPOINTCHARLIE","WIEDERVEREINIGUNG","BRANDENBURGERTOR",
 "FERNSEHTURM","NEUNTERNOVEMBER","MITTERNACHT","SCHATTEN","NORDOST","NORDOSTEN","OSTEN","OST","NORDEN","GEHEIM","SCHLUESSEL",
 "GEHEIMNIS","VERSCHLUESSELUNG","UHRZEIT","ZEITZONE","ZEITZONEN","WINDROSE","KOMPASS","ERICHJOHN","JOHN","SONNENSYSTEM",
 "PLANETEN","STASI","SPREE","ICHBINEINBERLINER","BERLINERLUFT","LICHT","DUNKEL","UNSICHTBAR","SCHICHT","SCHICHTZWEI",
 "KRYPTOS","BERLINUHR","MENGENLEHREUHR","WELTZEITUHRBERLIN","NEUNZEHNHUNDERTNEUNUNDACHTZIG","ACHTZEHN","EINS","ZWEI"]
def keyed(w): return "".join(dict.fromkeys(w + k4.AZ))
alphs = {"AZ": k4.AZ, "KA": k4.KA}
for w in DEKW: alphs["de:" + w] = keyed(w); alphs["derev:" + w] = keyed(w)[::-1]
pairs = {f"{p}|{c}": (alphs[p], alphs[c]) for p in alphs for c in alphs if p.startswith("de") or c.startswith("de")}
h = []; tol.tol_check(np.arange(97), k4.CT, k4.CPOS, k4.CPL, lambda i: "none", h, pairs, maxerr=2, mineff=8)
ex = []
for name, (PA, CA) in pairs.items():   # exact (0-error) check
    for sign in (1, -1):
        k = [(CA.index(k4.CT[i]) - sign * PA.index(c)) % 26 for i, c in zip(k4.CPOS, k4.CPL)]
        for L in range(1, 27):
            d = {}
            if all(d.setdefault(i % L, v) == v for i, v in zip(k4.CPOS, k)) and 24 - len(d) >= 8: ex.append((name, sign, L))
print("German-keyword alphabets (no transposition): exact", ex[:10], "| 1-2 errors", h[:10])
# keyword columnar widths 2-30 and double columnar (German x all keywords), exact + tolerant
allkw = DEKW + KEYWORDS
ords = {}
for kw in allkw:
    for w in sorted({len(kw)} | set(range(2, 31))) if kw in DEKW else [len(kw)]:
        ords[f"{kw}@{w}"] = np.array(k4.order_from(w, round2.kw_order(kw, w)))
names = list(ords); B, Lb, hx, ht = [], [], [], []
def flush():
    D = np.array(B); round2.check_generic(D, k4.CT, k4.CPOS, k4.CPL, lambda i, L=list(Lb): L[i], hx)
    tol.tol_check(D, k4.CT, k4.CPOS, k4.CPL, lambda i, L=list(Lb): L[i], ht, tol.STD); B.clear(); Lb.clear()
for a in names:
    o = ords[a]
    for d, arr in (("enc", round2.dest_from_order(o)), ("dec", o)): B.append(arr); Lb.append(f"{a} {d}")
for a in [n for n in names if n.split("@")[0] in DEKW and int(n.split("@")[1]) == len(n.split("@")[0])]:
    for b in names:
        o = ords[a][ords[b]]
        for d, arr in (("enc", round2.dest_from_order(o)), ("dec", o)): B.append(arr); Lb.append(f"dbl {a}|{b} {d}")
        if len(B) >= 4000: flush()
flush()
print("German keyword single+double columnar:", len(B), "exact hits", hx[:10], "| 1-2 error hits", ht[:10])
