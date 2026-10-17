"""Key constant per carved line (K4 occupies lines 25-28: 4 + 31 + 31 + 31 letters). Optionally after/before simple
or columnar transpositions. Alphabets AZ/KA, vig/beau. Also: line-based key indexed by carved line of the PLAINTEXT slot."""
import numpy as np, itertools, k4, dumb
line = np.array([0]*4 + [1]*31 + [2]*31 + [3]*31)
def check(D, labels, out):
    dc = D[:, k4.CPOS]
    for an, X in k4.ALPHS.items():
        cta = np.array([X.index(c) for c in k4.CT]); ci = cta[dc]; pi = np.array([X.index(c) for c in k4.CPL])[None, :]
        for kn, k in (("vig", (ci - pi) % 26), ("beau", (ci + pi) % 26)):
            for model, r in (("key by carved line of ciphertext", line[dc]), ("key by line of plaintext slot", np.broadcast_to(line[k4.CPOS][None, :], k.shape))):
                ok = ~((r[:, :, None] == r[:, None, :]) & (k[:, :, None] != k[:, None, :])).any(axis=(1, 2))
                for i in np.nonzero(ok)[0]: out.append((labels(i), an, kn, model))
hits = []
check(np.arange(97)[None, :], lambda i: "no transposition", hits); print("no transposition:", hits)
lib = dumb.library(); names = list(lib); orders = np.array([lib[n] for n in names])
dests = np.empty_like(orders); dests[np.arange(len(names))[:, None], orders] = np.arange(97)
check(dests, lambda i: names[i] + " enc", hits); check(orders, lambda i: names[i] + " dec", hits)
print("simple transpositions:", len(hits), hits[:5])
for w in range(2, 10):
    perms = list(itertools.permutations(range(w)))
    for s in range(0, len(perms), 20000):
        ch = perms[s:s+20000]; o = np.array([k4.order_from(w, p) for p in ch]); d = np.empty_like(o); d[np.arange(len(ch))[:, None], o] = np.arange(97)
        check(d, lambda i, ch=ch, w=w: f"w{w} enc {ch[i]}", hits); check(o, lambda i, ch=ch, w=w: f"w{w} dec {ch[i]}", hits)
print("columnar up to 9:", len(hits), hits[:5])
