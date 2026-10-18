"""Autokey (both sides, lags 1-96, 4 value alphabets) x 613k Sanborn-phrase alphabets."""
import time, phrases, autokey_dict as ak
from gromark_dict import straight, transposed
al = {}
for p in phrases.phrase_list():
    for f, a in (("s", straight(p)), ("t", transposed(p))):
        al.setdefault(a, f"{p}/{f}"); al.setdefault(a[::-1], f"{p}/{f}-rev")
AL = list(al.items()); t = time.time(); h = ak.scan(AL)
print(f"autokey x {len(AL):,} phrase alphabets: {len(h)} exact fits ({time.time()-t:.0f}s)", flush=True); [print("  ", x) for x in h[:20]]
