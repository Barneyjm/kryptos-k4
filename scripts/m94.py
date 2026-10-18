"""M-94 (Jefferson wheel) test: 25 published disks, secret disk order, 25-letter lines, per-line generatrix offset g.
C_i = disk[(pos(P_i) + g_line) mod 26] with the disk at column (i+s) % 25. Consistent iff a distinct-disk assignment exists."""
import itertools, random, k4
DISKS = """ABCEIGDJFVUYMHTQKZOLRXSPWN ACDEHFIJKTLMOUVYGZNPQXRWSB ADKOMJUBGEPHSCZINXFYQRTVWL AEDCBIFGJHLKMRUOQVPTNWYXZS
AFNQUKDOPITJBRHCYSLWEMZVXG AGPOCIXLURNDYZHWBJSQFKVMET AHXJEZBNIKPVROGSYDULCFMQTW AIHPJOBWKCVFZLQERYNSUMGTDX
AJDSKQOIVTZEFHGYUNLPMBXWCR AKELBDFJGHONMTPRQSVZUXYWIC ALTMSXVQPNOHUWDIZYCGKRFBEJ AMNFLHQGCUJTBYPZKXISRDVEWO
ANCJILDHBMKGXUZTSWQYVORPFE AODWPKJVIUQHZCTXBLEGNYRSMF APBVHIYKSGUENTCXOWFQDRLJZM AQJNUBTGIMWZRVLXCSHDEOKFPY
ARMYOFTHEUSZJXDPCWGQIBKLNV ASDMCNEQBOZPLGVJRKYTFUIWXH ATOJYLFXNGWHVCMIRBSEKUPDZQ AUTRZXQLYIOVBPESNHJWMDGFCK
AVNKHRGOXEYBFSJMUDQCLZWTIP AWVSFDLIEBHKNRJQZGMXPUCOTY AXKWREVDTUFOYHMLSIQNJCPGBZ AYJPXMVKBQWUGLOSTECHNZFRID AZDNBUHYFWJLVGRCQMPSOEXTKI""".split()
assert all(sorted(d) == list(k4.AZ) for d in DISKS) and len(DISKS) == 25
OFF = [{(p, c): (d.index(c) - d.index(p)) % 26 for p in k4.AZ for c in k4.AZ} for d in DISKS]
def matching(cands):
    match = {}
    def aug(col, seen):
        for d in cands[col]:
            if d in seen: continue
            seen.add(d)
            if d not in match or aug(match[d], seen): match[d] = col; return True
        return False
    return all(aug(col, set()) for col in cands)
def test(ct, cribs, s, gs):
    """gs: dict line -> offset. Returns True if a disk assignment exists."""
    cols = {}
    for i, p in cribs.items():
        line, col = (i + s) // 25, (i + s) % 25
        if line not in gs: return None
        cols.setdefault(col, []).append((p, ct[i], gs[line]))
    cands = {col: [d for d in range(25) if all(OFF[d][(p, c)] == g for p, c, g in req)] for col, req in cols.items()}
    if any(not v for v in cands.values()): return False
    return matching(cands)
def search(ct, cribs):
    hits = []
    for s in range(25):
        lines = sorted({(i + s) // 25 for i in cribs})
        for g in itertools.product(range(1, 26), repeat=len(lines)):
            if test(ct, cribs, s, dict(zip(lines, g))): hits.append((s, dict(zip(lines, g))))
    return hits
if __name__ == "__main__":
    CR = dict(zip(k4.CPOS.tolist(), k4.CPL))
    # planted test: random disk order, random offsets per line, start shift 1
    random.seed(4); order = random.sample(range(25), 25); g = [5, 17, 9, 22, 3]; s = 1
    P = [random.choice(k4.AZ) for _ in range(97)]
    for st, w in k4.CRIBS.items(): P[st:st+len(w)] = w
    ct = "".join(DISKS[order[(i+s) % 25]][(DISKS[order[(i+s) % 25]].index(p) + g[(i+s)//25]) % 26] for i, p in enumerate(P))
    h = search(ct, CR); print("planted found:", any(x[0] == 1 and list(x[1].values()) == [g[(st+s)//25] for st in sorted({i for i in CR})][:0] or x[0] == 1 for x in h), "| hits on planted:", len(h))
    h = search(k4.CT, CR); print("K4: consistent M-94 settings:", len(h), h[:10])
    # also: disks turned the other way (read rows upward) = negative offsets, covered since g ranges 1..25 mod 26 both ways
