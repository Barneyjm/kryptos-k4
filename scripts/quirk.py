import numpy as np, itertools, k4, dumb
def variants(base, tag):
    out = {}
    for i, j in itertools.permutations(range(26), 2):          # move letter at i to position j
        l = list(base); ch = l.pop(i); l.insert(j, ch); out[f"{tag} move {base[i]}->{j}"] = "".join(l)
    for i, j in itertools.combinations(range(26), 2):          # swap two letters
        l = list(base); l[i], l[j] = l[j], l[i]; out[f"{tag} swap {base[i]}{base[j]}"] = "".join(l)
    return out
ALPH = {**variants(k4.KA, "KA"), **variants(k4.AZ, "AZ")}
uniq = {}
for k, v in ALPH.items(): uniq.setdefault(v, k)
ALPH = {k: v for v, k in uniq.items() if v not in (k4.KA, k4.AZ)}
print("quirky alphabets:", len(ALPH), flush=True)
lib = dumb.library(); names = ["identity"] + [n for n in lib if n != "identity"]
orders = np.array([lib[n] for n in names]); dests = np.empty_like(orders); dests[np.arange(len(names))[:, None], orders] = np.arange(97)
D = np.vstack([dests, orders]); labels = [n + " enc" for n in names] + [n + " dec" for n in names]
dc = D[:, k4.CPOS]; hits = []
for an, X in ALPH.items():
    idx = np.array([X.index(c) for c in k4.AZ])
    ci = idx[np.array([ord(c) - 65 for c in k4.CT])][dc]; pi = idx[np.array([ord(c) - 65 for c in k4.CPL])][None, :]
    for kn, k in (("vig", (ci - pi) % 26), ("beau", (ci + pi) % 26)):
        for model in "AB":
            for L in range(1, 27):
                r = np.broadcast_to(k4.CPOS[None, :] % L, k.shape) if model == "A" else dc % L
                ok = ~((r[:, :, None] == r[:, None, :]) & (k[:, :, None] != k[:, None, :])).any(axis=(1, 2))
                for i in np.nonzero(ok)[0]:
                    cons = 24 - len(set(r[i].tolist()))
                    if cons >= 8: hits.append((labels[i], an, kn, model, L, cons))
print("hits:", len(hits)); [print(h) for h in hits[:20]]
