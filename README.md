# Kryptos K4: hypothesis testing notes

A record of a structured attempt at K4, the unsolved fourth passage of Jim Sanborn's *Kryptos* sculpture at CIA headquarters. The work was done over a few days in October 2026. **K4 is not solved here.** This repo documents what was tested, how, and what was ruled out, so others don't have to repeat it.

Everything is checked against the 24 plaintext letters Sanborn has released (0-indexed positions 21–33 `EASTNORTHEAST` and 63–73 `BERLINCLOCK`). Every tester was first validated by planting a known solution and confirming it was recovered.

## Data
- `data/panel.txt`: full left-panel ciphertext (Elonka Dunin's transcript). It was verified by decrypting K1 (PALIMPSEST) and K2 (ABSCISSA) and by recovering K3's double-rotation method.
- `data/split.txt`: K1, K2, K3 and K4 ciphertexts, one per line.
- `data/tableau.txt`: the right-panel tableau (867 letters, including the extra `L`).

## Results summary

### Ruled out (no transposition)
- Periodic Vigenère and Beaufort with key lengths 1–26, standard and KRYPTOS alphabets, plus about 13k keyword-alphabet (Quagmire) pairs. Exact matches and matches with 1–2 errors tolerated.
- Hand-modified alphabets: 1,850 variants with one letter moved or two letters swapped.
- Autokey (plaintext and ciphertext), progressive keys, polynomial keys mod 26/24/12/60, and digit streams (π, e, √2, φ, K2 coordinates, 1986/1989).
- Hill ciphers 2×2, 3×3 and 4×4, at every alignment, with and without the `?`. The 4×4 keys forced by the cribs are all non-invertible.
- Any linear rule mod 26 combining up to 5 of: P[i], P[i±1..2], C[i±1..2], i, i², constant (1,544 rules).
- Compound keys made by adding two unknown periodic keys of lengths 2–26 each (Scheidt's "key split" idea). Also clue-word pairs and triples at every phase, and QUAE plus two unknown splits.
- Grid keys (column key plus row key) at widths 2–30 and every offset.
- Running keys from K1–K3 plaintext and ciphertext, the tableau (read 7 different ways), and the Morse text.
- Gromark family with both alphabets free: bases 3–26, primers of 2–8 digits, every two-lag recurrence (2.78 billion primers). This reproduces Bean (2021) exactly: 39 survivors for base 10 with a 5-digit primer, and 00351/00537 for base 8. None of the 804,684 survivors explains the width-21 fingerprint (the best explains 3 of 11 repeated pairs, which is chance level).
- Raised letters (DYAHR, YAR) and clue numbers as primers.
- Running key (English or German text) with keyword-mixed alphabets: 842,006 dictionary-keyword alphabets (straight, transposed, reversed) in every Quagmire arrangement. K4's best score (-5.17 en / -5.10 de) falls inside the range reached by shuffled ciphertext (-4.86 to -5.21).
- Digit keys (Gronsfeld/Gromark/clock digits) with dictionary-keyword alphabets: no alphabet brings all 24 implied key values below 16.
- Chain addition (Vimark/VIC-style) primed by any of 221,577 dictionary words or clue words, every two-lag rule, both alphabets free: 6,976 survive the cribs by chance, and none explains more than 2 of the 11 width-21 repeats.
- Phrase keys from Sanborn's own words: 192,569 runs of 1-5 words from his 2009 Smithsonian oral history (archived copy in data/), K1-K3, Morse and clue quotes, giving 613,440 phrase alphabets. Tested as Quagmire alphabets with periodic keys up to 26 (exact), as alphabets for an English running key (K4 -5.29 vs shuffled -4.99 to -5.22), as digit-key alphabets (none <= 12), and as the full-length key itself (none >= 15/24).
- Arbitrary monoalphabetic mask (any permutation, solved algebraically) on the plaintext or ciphertext side, plus a periodic key: inconsistent for L = 1-12, 14, 15, 17, 21, 22 (cipher side also 18, 25). L = 13, 16, 19, 20, 23, 24, 26 are underdetermined with 24 cribs.
- Crib extension: 3,497 candidate words (top English, Sanborn's vocabulary, direction/measure words) at the 4 crib edges, re-testing the underdetermined mask periods. Real words fit less often than random strings (24 vs 54 of 503,568), and no fit decrypts to readable text.
- "Divided Berlin": separate periodic keys for the two halves, with the same or independent alphabets (13,456 pairs). Neither half alone fits a key of length 1-9 (east) or 1-7 (west).

### Ruled out (with transposition)
- Columnar transposition (both directions) followed or preceded by periodic substitution, exhaustively for widths 2–11. Widths 12–13 were covered by a pruned search with key period ≤ 14.
- 1,282 "simple" hand transpositions: reversal, boustrophedon and column/row reads at widths 2–48, the carved 31-wide layout, decimations, rail fence. Each was tested with periodic substitution, Hill, Quagmire and running keys.
- K3-style double grid rotations of K4, `?K4` and padded K4.
- Double keyword columnar using about 56 English and 56 German clue words.
- English and German running keys behind columnar transposition up to width 11. The best scores equal the null scores obtained from shuffled ciphertext.

### Statistical fingerprints (shuffle tests)
| Property | Raw | Corrected |
|---|---|---|
| Repeated vertical bigrams at width 21 (11) | 1 in ~4,000 | ~1 in 210 (best of widths 2–48) |
| Same, excluding the doublet-stacked pairs (9) | 1 in 175 | ~1 in 10 |
| 5 of 6 doublets at the same position mod 7 | 1 in 2,500 | ~1 in 4 (any modulus) |
| "Minor differences", KRYPTOS-letter subset | 1 in 4,500 | n/a, because the subset was post hoc |
| "Minor differences", all 24 known letters | 0.14 | not significant |
| IC and other ciphertext-only statistics | consistent with random | |

The one property that survives correction is the 7/21 rhythm, and it is carried mostly by the doublets. Within the known plaintext, every visible step from column 5 to column 6 (at width 7) is a ciphertext doublet: NO→QQ, ST→SS, IN→TT. No grid key, compound key or Gromark key reproduces this.

## Note on transposition
Sanborn has described BERLINCLOCK as a direct, position-by-position match with NYPVTTMZFPK ("plain B taken, has the encipherment done to it, and out pops a cipher N"). Taken literally, there is no transposition. The transposition tests above are kept for completeness.

## Open leads
- A masking step applied to the plaintext before substitution (Scheidt: "change the language base", "more than one stage").
- The 7-rhythm combined with Sanborn's statement that the cipher is "one-to-one" (no transposition).

## Running
The scripts live in `scripts/` and expect `panel.txt`, `split.txt`, `tableau.txt` and the corpora in the working directory. Run `scripts/fetch_corpora.sh` once to download the Project Gutenberg texts used for the n-gram scores. Most scripts use numpy and multiprocessing; the Hill tests need sympy, and the digit streams need mpmath.

## References
- R. Bean, "Cryptodiagnosis of Kryptos K4", HistoCrypt 2021. https://github.com/RichardBean/k4testing
- E. Dunin, Kryptos transcript. https://elonka.com/kryptos/transcript.html
- Paradigm, K4 guess verification (2026). https://cn.paradigm.xyz/2026/06/kryptos
