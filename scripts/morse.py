import numpy as np, itertools, k4, round2, tol, dumb
MORSE = "VIRTUALLYINVISIBLEDIGETALINTERPRETATIUSHADOWFORCESLUCIDMEMORYTISYOURPOSITIONSOSRQ"
orders_ = list(itertools.permutations(["VIRTUALLYINVISIBLE","DIGETALINTERPRETATIU","SHADOWFORCES","LUCIDMEMORY","TISYOURPOSITION","SOS","RQ"]))
SRC = {"morse": MORSE, "morse_rev": MORSE[::-1], "morse_alt": "SHADOWFORCESLUCIDMEMORYTISYOURPOSITIONSOSRQVIRTUALLYINVISIBLEDIGETALINTERPRETATIU"}
SRC["morse_YR"] = MORSE.replace("RQ", "YR")
# running-key exact check, no transposition + simple library (model A & B)
lib = dumb.library(); hits = []
names = list(lib); O = np.array([lib[n] for n in names]); Dd = np.empty_like(O); Dd[np.arange(len(names))[:, None], O] = np.arange(97)
DEST = np.vstack([Dd, O]); LAB = [n + " enc" for n in names] + [n + " dec" for n in names]
for an, X in k4.ALPHS.items():
    cta = np.array([X.index(c) for c in k4.CT]); pi = np.array([X.index(c) for c in k4.CPL]); ci = cta[DEST[:, k4.CPOS]]
    for kn, need in (("vig", (ci - pi) % 26), ("beau", (ci + pi) % 26), ("vb", (pi - ci) % 26)):
        for sn, s_ in SRC.items():
            kk = np.array([X.index(c) for c in s_ * 4]); offs = np.arange(len(s_))
            for model in "AB":
                pos = np.broadcast_to(k4.CPOS, need.shape) if model == "A" else DEST[:, k4.CPOS]
                m = (kk[offs[:, None, None] + pos[None]] == need[None]).sum(axis=2)
                for o_, t in zip(*np.nonzero(m >= 20)): hits.append((LAB[t], an, kn, sn, model, int(o_), int(m[o_, t])))
print("Morse running key hits:", hits[:10])
# Morse phrases as keywords for alphabets and columnar orders
PH = ["VIRTUALLYINVISIBLE","DIGETALINTERPRETATIU","SHADOWFORCES","LUCIDMEMORY","TISYOURPOSITION","SOSRQ","YOURPOSITION","SHADOWFORCESLUCIDMEMORY"]
keyed = lambda w: "".join(dict.fromkeys(w + k4.AZ))
alph = {"AZ": (k4.AZ, k4.AZ), "KA": (k4.KA, k4.KA)}
for p in PH: alph[p] = (keyed(p), keyed(p)); alph["AZ|" + p] = (k4.AZ, keyed(p)); alph["KA|" + p] = (k4.KA, keyed(p))
h2 = []; tol.tol_check(np.arange(97), k4.CT, k4.CPOS, k4.CPL, lambda i: "none", h2, alph, maxerr=2, mineff=8)
h3 = []
for p in PH:
    for w in range(2, 31):
        o = np.array(k4.order_from(w, round2.kw_order(p, w)))
        for d, arr in (("enc", round2.dest_from_order(o)), ("dec", o)):
            round2.check_generic(arr, k4.CT, k4.CPOS, k4.CPL, lambda i, p=p, w=w, d=d: f"{p}@{w} {d}", h3)
print("Morse-phrase alphabets (0-2 errors):", h2[:5], "| Morse-phrase columnar:", h3[:5])
