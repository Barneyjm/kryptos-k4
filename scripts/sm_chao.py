"""Chaocipher (Byrne 1918; algorithm as published by Rubin 2010).
Left disk = ciphertext alphabet, right disk = plaintext alphabet (standard roles).
Per letter: p found in right at index i, c = left[i]. Then
 left : rotate so c at zenith (index 0); take letter at 1, shift 2..13 up, insert at 13.
 right: rotate so p at zenith, then one more step (p goes to index 25); take letter at 2,
        shift 3..13 up, insert at 13.
Both permutations are written as gathers: new = old[(i + SIG) % 26]."""
import numpy as np
SIG_L = np.array([0] + list(range(2, 14)) + [1] + list(range(14, 26)))
SIG_R = (np.array([0, 1] + list(range(3, 14)) + [2] + list(range(14, 26))) + 1)   # +1 = the extra right shift

def _perm(alph, i, sig):
    return "".join(alph[(i + s) % 26] for s in sig)

def encrypt(pt, left, right):
    out = []
    for p in pt:
        i = right.index(p); c = left[i]; out.append(c)
        left = _perm(left, i, SIG_L); right = _perm(right, i, SIG_R)
    return "".join(out)

def decrypt(ct, left, right):
    out = []
    for c in ct:
        i = left.index(c); p = right[i]; out.append(p)
        left = _perm(left, i, SIG_L); right = _perm(right, i, SIG_R)
    return "".join(out)

# swapped roles: ciphertext read from the RIGHT disk (permuted with the right rule), plaintext on the left.
def decrypt_swapped(ct, left, right):
    out = []
    for c in ct:
        i = right.index(c); p = left[i]; out.append(p)
        left = _perm(left, i, SIG_L); right = _perm(right, i, SIG_R)
    return "".join(out)

def encrypt_swapped(pt, left, right):
    out = []
    for p in pt:
        i = left.index(p); c = right[i]; out.append(c)
        left = _perm(left, i, SIG_L); right = _perm(right, i, SIG_R)
    return "".join(out)

def batch_decrypt(ct_idx, L, R, nsteps, swapped=False):
    """Vectorized. ct_idx: int array of ciphertext letter codes (0..25) in the alphabet-symbol space.
    L, R: (n,26) uint8 arrays of letter codes. Returns (n,nsteps) plaintext codes."""
    n = L.shape[0]; rows = np.arange(n)[:, None]
    out = np.empty((n, nsteps), np.uint8)
    for t in range(nsteps):
        c = ct_idx[t]
        if not swapped:
            i = np.argmax(L == c, axis=1); out[:, t] = R[np.arange(n), i]
        else:
            i = np.argmax(R == c, axis=1); out[:, t] = L[np.arange(n), i]
        L = L[rows, (i[:, None] + SIG_L[None, :]) % 26]
        R = R[rows, (i[:, None] + SIG_R[None, :]) % 26]
    return out

if __name__ == "__main__":
    L0 = "HXUCZVAMDSLKPEFJRIGTWOBNYQ"; R0 = "PTLNBQDEOYSFAVZKGJRIHWXUMC"
    pt = "WELLDONEISBETTERTHANWELLSAID"; exp = "OAHQHCNYNXTSZJRRHJBYHQKSOUJY"
    got = encrypt(pt, L0, R0)
    print("test vector enc:", got, "OK" if got == exp else "MISMATCH")
    print("test vector dec:", decrypt(exp, L0, R0), "OK" if decrypt(exp, L0, R0) == pt else "MISMATCH")
    # vectorized vs reference, both role assignments, random alphabets
    rng = np.random.default_rng(1)
    AZ = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    ok = True
    for trial in range(200):
        l = "".join(rng.permutation(list(AZ))); r = "".join(rng.permutation(list(AZ)))
        msg = "".join(rng.choice(list(AZ), 60))
        for sw, E, D in ((False, encrypt, decrypt), (True, encrypt_swapped, decrypt_swapped)):
            ct = E(msg, l, r); assert D(ct, l, r) == msg
            Lb = np.array([[ord(x) - 65 for x in l]], np.uint8); Rb = np.array([[ord(x) - 65 for x in r]], np.uint8)
            o = batch_decrypt(np.array([ord(x) - 65 for x in ct]), Lb, Rb, 60, sw)
            ok &= "".join(chr(65 + v) for v in o[0]) == msg
    print("vectorized == reference (400 random round trips):", ok)
