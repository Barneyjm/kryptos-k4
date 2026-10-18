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

## Session 2 (night of 2026-10-17): gap-filling

All exact unless noted; every tester validated on a planted solution first. Scripts in `scripts/`, logs in `results/`.

- **Correction.** Earlier periodic-key tests only reported fits with >= 8 independent checks, which silently excluded key lengths 17, 18, 23-26. Re-tested (`longper.py`): lengths 17, 18, 23, 24 give 0 fits across 1.8M alphabets; 25 and 26 are inconsistent with the cribs for every standard alphabet pairing.
- Linear key models (repeating 1-26, progressive, grid keys W=2-30 all offsets, two combined keys, linear/quadratic formulas; 288 models with >= 8 checks) x 842,006 dictionary-keyword alphabets x 10 arrangements x vig/beau (`cokernel.py`, `gap1.py`): 0 fits; two shuffled baselines also 0.
- The same models x 959,024 further alphabets (decimated A-Z/KA, keyword + reversed remainder, keyword + continue-from-last, 19,535 German words in four constructions) (`gap2.py`): 0.
- Autokey (plaintext and ciphertext, lags 1-96, key value read in plain/cipher/A-Z/KRYPTOS alphabets) x dictionary alphabets (`autokey_dict.py`): 0 exact. Near-miss version (`autokey_near.py`, 707M configurations): K4's best is 10/24; shuffled ciphertext reaches 10-11/24. No signal.
- M-94 wheel cipher (25 published disks, all disk orders via bipartite matching, all line offsets, all alignments) (`m94.py`): 0.
- Arbitrary monoalphabetic mask + autokey, either side, all lags (`mask_autokey.py`): 0.
- Three-term chain addition, bases 3-26, primers <= 8 digits, both alphabets free (`chain3.py`, 2.27B primers): 790,052 chance survivors; the best explains 3 of the 11 width-21 repeats (chance level).
- Physical layout as key source: tableau letter at the same row/column (plain and mirrored, row shifts +-3, column shifts +-6) (`overlay.py`); letter 1-27 rows above/below in the same column on either half (`verticals.py`, 4,368 configurations): 0 exact, best 5/24 (chance). Skip-codes (every k-th letter) through the combined panels: at noise level.
- Keys derived from K1-K3 by rule (initials, finals, letters after X, word lengths, every n-th letter; 90 strings) as running/periodic keys and as keyword alphabets (`derived.py`): 0.
- Backward autokey (key from 1-96 letters later in the ciphertext or plaintext), std + clue + dictionary alphabets (`backak.py`): 0 exact; near-miss tail at noise level.
- Running key taken straight from a text, by exact lookup (`rkidx_run.py`): corpus of 8.9M letters (5 English and 5 German books, Sanborn's 2009 interview, K1-K3, panel, tableau, and Carter & Mace *The Tomb of Tut-ankh-Amen* vol. 1, `data/carter_tomb_vol1.txt`), forward and reversed, x 1,801,146 alphabets x 5 arrangements x 3 conventions x 3 indicator alphabets = 81M configurations: 0 full hits, 0 single-block (13- or 11-letter) hits; shuffled baselines also 0.
- One-slip tolerance (`slip1.py`): repeating (1-26) and progressive keys with any single crib letter treated as an encryption/carving error (866 model x dropped-crib variants) x 1,801,032 alphabets: 0 fits. Validated by recovering a planted key with one corrupted crib.
- Chaocipher / self-modifying alphabets (`sm_chao.py`, `sm_taskA.py`; agent-run, validated on two planted messages, reproduces the Rubin test vector): 116x116 clue alphabets, every dictionary+clue alphabet (842,054) on either disk against AZ/KA and their reverses, and 1M random dictionary pairs, x 26 relative rotations x 2 directions (402M rows): 0 exact fits; K4 best partial 8-9/24 vs shuffled 8-10/24. Not covered: clue x dictionary pairs (~1e10 rows) and non-keyword alphabets.
- Clock-derived additive key streams (`sm_taskB.py`; agent-run, three planted messages recovered): 1.35M streams (Berlin-clock lamp counts/rows/bits/HHMM at every start minute and step, time-zone offset cycles, K2 coordinate and date digit cycles) x 8.45M implied-key vectors (dictionary/clue alphabets as Quagmire I-IV, Vigenere and Beaufort) x sign x offset: 0 exact fits on K4 and on both shuffles; best partial 10/24 for K4 and for the shuffles.
- Arbitrary mask + structured keys (`mg_mask_grid.py`; agent-run, validated on a planted W=9 grid-key message): grid keys (row + column key, W=2-30, column key free or repeating H=2-13), two combined keys (L1<L2<=26), progressive keys, and keyword-through-grid keys, each with an arbitrary injective mask on either side and A-Z/KA numbering (14,408 models). Exactly inconsistent: 9,382. Consistent with an injective mask: 178, of which 126 pass only because positions 27 and 65 are both R->P (one shared key value), and the other 52 carry a single parity check each (chance ~1/26 or 1/2). Shuffled ciphertexts give 4-91 such survivors, so nothing stands out. The rest are undecidable with 24 cribs.
- Sanborn-phrase alphabets (613,440, from the 2009 oral history, clue quotes and K1-K3) x autokey (both sides, lags 1-96) (`gap_phr_ak.py`): 0 exact fits.
- Two-stage polyalphabetic with different tableaux (nonlinear composition; solved by connected components) (`twostage.py`): std/KA alphabets in all roles, L1,L2 <= 26: 0 significant; clue-word alphabets as two different tableaux, L1,L2 <= 10 (5.3M configurations): 0. Shuffled baseline 0.
- Two-stage, dictionary tableaux (`ts2_*.py`, vectorized component solver): stage 1 on AZ/KA + stage 2 on any of 842,006 dictionary alphabets and the reverse (1.8e9 alphabet-configs per role, L1,L2 <= 12), and 56 clue alphabets x 421,003 dictionary alphabets (1.26e10): K4 significant-hit counts equal the shuffled baseline (102 vs 120; 673 vs 665), none at net >= 6. A planted PALIMPSEST/ABSCISSA two-stage message is recovered at net 6 and fully decrypted.
- Fresh-eyes pass (`fe_*.py`, ideas in `results/fe_ideas.md`): interrupted/event-driven keyword keys (restart, advance or pause on letter triggers, words, carved lines; 11,077 distinct models x 842k alphabets x 5 arrangements, plus Quagmire IV over 116^2 clue pairs): 0. Wheatstone cryptograph (27-cell outer ring vs 26-letter inner ring, 21.9M outer rings, inner solved exactly): 0. Base-27 additive ciphers (Z27, Z3^3, Z9xZ3; 22.7M alphabet orders; periodic 1-23 and autokey keys): 0. Independent substitution per position class: impossible for lengths 1-7, 9, 10, 14, 15, 17, 21, 25, because positions 24/66 are both V (from T, L) and 31/73 are both K (from A, K): no scheme may reuse a substitution 42 positions apart, which excludes periods 7, 14, 21, 42. Affine / mixed Vigenere-Beaufort keys: 0 for lengths 1-10; undecidable beyond.
- Information theory (`bm.py`, `bm_dict.py`): K4's ciphertext entropy (4.55 bits/letter in 97-letter windows) is at the random ceiling, above K1/K2 (4.35-4.39, consistent with short repeating keys) and K3 (3.95, English). The implied key runs have random-like linear complexity (Berlekamp-Massey mod 2 and mod 13) under all 8.42M alphabet/type variants; no single linear recurrence generates both runs. If a short memorable key produced K4, the expansion step is nonlinear.
- Paradigm CTF puzzles (`ctf_puzzles.json`, `solve_*.py`): PK1 = KRYPTOS-alphabet Vigenère, key PROVENANCE; PK2 = columnar transposition width 7; PK3 = KRYPTOS-alphabet Vigenère with a 40-letter random key. All ten were already solved by others (prizes claimed June-Oct 2026).

## Running
The scripts live in `scripts/` and expect `panel.txt`, `split.txt`, `tableau.txt` and the corpora in the working directory. Run `scripts/fetch_corpora.sh` once to download the Project Gutenberg texts used for the n-gram scores. Most scripts use numpy and multiprocessing; the Hill tests need sympy, and the digit streams need mpmath.

## References
- R. Bean, "Cryptodiagnosis of Kryptos K4", HistoCrypt 2021. https://github.com/RichardBean/k4testing
- E. Dunin, Kryptos transcript. https://elonka.com/kryptos/transcript.html
- Paradigm, K4 guess verification (2026). https://cn.paradigm.xyz/2026/06/kryptos
