"""Key = (n unknown letters) + WORD + (unknown) repeating with period L. Test against cribs, with optional transposition."""
import numpy as np, k4, dumb
WORDS = ["QUA","QUAE","QUAEU","EQUA","AQU","UQA","AUQ","EAUQ","QUAU","LOE","LOEI","LOEION","IQLUSION","UNDERGRUUND","DESPARATLY",
         "DIGETAL","INTERPRETATIU","XLAYERTWO","IDBYROWS","PALIMPSEST","ABSCISSA","KRYPTOS","BERLINCLOCK","WELTZEITUHR"]
F = {"vig": lambda c, p: (c - p) % 26, "beau": lambda c, p: (c + p) % 26, "vb": lambda c, p: (p - c) % 26}
def run(dest, label, out):
    for an, X in k4.ALPHS.items():
        for fn, f in F.items():
            need = [(int(dest[i]), f(X.index(k4.CT[dest[i]]), X.index(p))) for i, p in zip(k4.CPOS, k4.CPL)]
            for model in "AB":                    # A: key indexed by plaintext pos; B: by ciphertext pos
                req = [((i if model == "A" else d), v) for (d, v), i in zip(need, k4.CPOS.tolist())]
                for W in WORDS:
                    wv = [X.index(c) for c in W]
                    for L in range(len(W), 97):
                        base = {}; ok = True; cc = 0
                        for pos, v in req:
                            r = pos % L
                            if r in base:
                                if base[r] != v: ok = False; break
                                cc += 1
                            else: base[r] = v
                        if not ok: continue
                        for n in range(L - len(W) + 1):
                            m = 0; bad = False
                            for j, v in enumerate(wv):
                                r = (n + j) % L
                                if r in base:
                                    if base[r] != v: bad = True; break
                                    m += 1
                            if not bad and m >= min(len(W), 4) and cc + m >= 7:
                                out.append((cc + m, label, an, fn, model, W, L, n, cc, m))
if __name__ == "__main__":
  out = []
  run(np.arange(97), "none", out)
  print("no transposition:", sorted(out, reverse=True)[:10])
  lib = dumb.library(); hits = []
  for k, o in lib.items():
      o = np.array(o); d = np.empty(97, int); d[o] = np.arange(97)
      run(d, k + " enc", hits); run(o, k + " dec", hits)
  hits.sort(reverse=True); print("with simple transpositions:", len(hits)); [print(h) for h in hits[:15]]
