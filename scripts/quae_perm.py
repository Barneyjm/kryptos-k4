import itertools, quae
quae.WORDS[:] = sorted({"".join(p) for p in itertools.permutations("QUAE")} | {"".join(p) for p in itertools.permutations("QUA")})
exec(open("quae.py").read().split("out = []", 1)[1].join(["out = []", ""]) if False else "")
import numpy as np, dumb
out = []; quae.run(np.arange(97), "none", out); print("words:", quae.WORDS); print("no transposition:", sorted(out, reverse=True)[:10])
hits = []
for k, o in dumb.library().items():
    o = np.array(o); d = np.empty(97, int); d[o] = np.arange(97); quae.run(d, k + " enc", hits); quae.run(o, k + " dec", hits)
print("with simple transpositions:", len(hits), sorted(hits, reverse=True)[:10])
