# K4: untested hand-cipher mechanisms (fe_ review, 2026-10-17)

Constraints used as filters: one-to-one, in place (no transposition); C[73] = P[73] = K must be possible;
"more than one stage", "change the language base", memorable key word(s), agent/pilot-style hand
procedure; implied key at the cribs looks random (high linear complexity) under every alphabet tried;
only fingerprint is the 7/21 doublet rhythm. Excluded from this list: anything in the README, the
"closed tonight" list, or the parallel workstreams (two-stage polyalphabetic with different tableaux,
keyword-through-grid, block-constant/drift keys, backward autokey, running-key hash lookup).

Quick fact used below: K4's ciphertext contains all 26 letters, so any 25-cell (5x5 Polybius,
I/J merged) system that outputs letters directly is impossible.

## 1. Interrupted / event-driven keyword keys (Friedman, *Military Cryptanalytics* II) — TESTED (T1)
A keyword is used normally, but the key position is driven by events instead of by the position:
restart the keyword after every ciphertext letter X (the receiver sees C, so it is decipherable),
let the key advance only after certain ciphertext letters, pause on them, restart after a plaintext
letter (word separator) or at every word, or at each carved line of the sculpture. This is the
classic "simple, memorable" way to break up the period, and it is exactly what makes a short keyword
produce an implied key with *random-looking linear complexity*: the period is not fixed, so
periodic/progressive/grid/split models all fail, as observed. Event triggers keyed to the ciphertext
can also correlate key positions with ciphertext letters, which could create doublets. Testability:
exact and cheap. Every variant reduces to "these crib positions share a key letter", a linear model,
so it runs over all 842k dictionary alphabets in Quagmire I-III arrangements.

## 2. Wheatstone cryptograph / 27-vs-26 differential ring — TESTED (T2)
Wheatstone's 1867 device, or its paper form: an outer ring of 26 letters plus a blank (27 cells) and
an inner mixed alphabet of 26. The pointer moves clockwise to each plaintext letter and the inner
ring moves by the same number of cells, so it gains one letter per revolution. That is a literal
"change of language base" (27 against 26), the procedure needs only one keyword (the inner
alphabet), and it is one-to-one and in place. The implied Vigenère key is a cumulative sum of ring
distances. That makes it non-periodic and nonlinear mod 26, so none of the linear or autokey
sweeps covered it. A plaintext letter can encrypt to itself (K->K). Testability: exact. For a
given outer ring, the inner offsets inside each crib are known. The inner alphabet can therefore be
left completely arbitrary and solved as a partial bijection, with one unknown shift between the
cribs. Cheap.

## 3. Additive cipher over a 27-symbol base (Z27, trifid-cube Z3^3, Z9xZ3) — TESTED (T3)
Put the 26 letters plus one symbol (space, '?', '+') into 27 cells, either as a line (mod-27
Vigenère) or as Delastelle's 3x3x3 trifid cube. Then add key coordinates digit by digit mod 3
without fractionating. This is a "Vigenère in base 3": each letter stays in place, the procedure is
hand-friendly with a keyword, and it is a literal change of base. Every Z26 alphabet sweep misses
it, because Z3^3 is not cyclic and no relabeling of Z26 can imitate it. The extra symbol never
appears in K4. With a random-like key that happens with p = (26/27)^97 = 2.6%. With a short
periodic key over English it is common, because each class has only one plaintext letter that
maps to the symbol. The "?" before K4 might also be that symbol. Testability: exact for periodic keys and autokey over every
keyword-built cell order. Cheap.

## 4. Periodic polyalphabetic with *independent* alphabets per position class (strip cipher M-138-A, DRYAD-like sheets)
Agents and pilots carried sheets of random alphabets (US Army M-138-A strips, DRYAD
authentication/numeral sheets). With period L, each residue class uses its own arbitrary
alphabet. Testing the published M-94 disks did not test unpublished strips. The necessary
condition for *any* such cipher is cheap and exact: within each residue class, plain->cipher must be
a partial bijection. Covering every period gives an alphabet-free list of the periods that are
still possible. That list is useful for every hidden-alphabet hypothesis. Testability: trivial
(ran as a side check in T1's output file, see fe_results.out).

## 5. Periodic affine / sign-alternating keys (Vigenère and Beaufort mixed by position)
C = a_i * P + b_i, with (a_i, b_i) periodic, or the special case a_i in {+1, -1}: alternate
Vigenère and Beaufort, for example by a 7-letter pattern ("KRYPTOS" vowels against consonants).
A mask plus a periodic key covers one fixed substitution, not a substitution that changes per
class. The rhythm in the doublets could come from a sign pattern of period 7. Testability: exact
and cheap for L <= 11, because each class needs 3 or more crib letters. For larger L it is
underdetermined.

## 6. Chaocipher-style self-modifying alphabets (Byrne) / alphabet permuted after every letter
Two keyed alphabets are permuted after every letter by a fixed hand rule. The cipher is one-to-one
and in place, and C = P is possible. The state is fully known from position 0 when the alphabets
are keyword-mixed, so a test is an exact full decryption per keyword pair. Byrne's rule was not
published until 2010, which makes it unlikely for 1989. Simpler "shift the cipher alphabet by one
after every word or letter" rules are covered by progressive keys. Testability: about 10^7 Python
decryptions for clue x dictionary pairs. Hours.

## 7. Berlin-Clock (Mengenlehreuhr) time-coded key, base 5 / base 12
The Weltzeituhr/Mengenlehreuhr shows time as lamps in base 5: 5-hour blocks, hours, 5-minute
blocks, minutes. A key stream taken from lamp counts or patterns of a sequence of times (for
example consecutive minutes from a start time, a "clock" that advances one tick per letter) is
both "change the base" and a direct link to the BERLIN CLOCK crib. The stream is quasi-periodic
with periods 5/4/11/4 and 60/24, and at the crib positions it would look random to linear
complexity tests that use a single alphabet. Testability: cheap once a reading convention is
fixed. The convention space (which lamps, start time, step size) is moderately large, so it needs
a planted test.

## 8. Four-square / two-square digraphic substitution with keyword squares (6x6 or 26+ cells)
These are digraphic and position-preserving, so they count as one-to-one in length. Unlike
Playfair they can map a letter to itself. 5x5 squares are excluded because all 26 letters occur,
which leaves 6x6 squares (letters plus digits) as the hand option. The crib pairs contain no
repeated bigram at either alignment, so a generic digraphic system is untestable. A
keyword-squares test works with the known square structure (row/column consistency), but Sanborn's
"plain B ... out pops N" wording argues against digraphic. Testability: hours, low prior.

## 9. "Pilot's card" letter -> digit -> letter (DRYAD-style numeral encipherment)
DRYAD rows map letters to digits 0-9 through a random row alphabet. A digit key is added and the
result is mapped back, which changes the language base to decimal while staying in place. The
26 -> 10 round trip needs homophones (several letters per digit), so the ciphertext could not be
decrypted uniquely unless the key stage restores the information. In its plain form it is
excluded: the decryption map C -> P would need to be a function per key value, and with digit keys
that is already covered by the closed "digit key x dictionary alphabet" sweep. Testability: mostly
covered. Low priority.

## 10. Mod-27/mod-28 Vigenère with the plaintext including word spaces (the '?' as 27th symbol)
This is a sub-case of 3 where the extra symbol is a word space in the plaintext. The cribs are
contiguous (no space inside EASTNORTHEAST or BERLINCLOCK), which restricts it to compound words.
It is covered by T3 when the extra cell sits in any of the 27 slots.

## Ranking for cheap testing
1 (interrupted keys), 2 (Wheatstone 27/26) and 3 (base-27 groups) were chosen. All three are
exact, cheap, plausible for 1989, and missed by every linear sweep over Z26. Idea 4 was added as
a free alphabet-independent side check. Ideas 5 and 7 are the next cheapest.

## Status after testing (numbers in fe_results.out)
- 1 Interrupted keys (fe_interrupt.py): 0 exact fits across 11,077 distinct models x 842k alphabets (chance 0.03). The
  alphabet-free version is uninformative (K4 falls inside the shuffled range).
- 2 Wheatstone 27/26 (fe_wheatstone.py): 0 of 21.9M structured outer rings, with the inner alphabet arbitrary. An arbitrary
  outer ring is underdetermined (about 10^7 or more solutions estimated).
- 3 Base-27 groups (fe_group27.py): 0 of 22.7M embeddings x 4 groups, periodic L <= 23 and autokey lags 1-26.
- 4 Independent alphabets (fe_indep.py): excluded for every period dividing 42 (1, 2, 3, 6, 7, 14, 21, 42), plus 4, 5, 9,
  10, 15, 17, 25. This holds for any alphabets.
- 5 Periodic affine / Vigenere-Beaufort mixing (fe_affine.py): excluded for L <= 10 over all alphabets tested.
- Next cheapest: 7 (Berlin-clock base-5 key streams) and 6 (self-modifying alphabets).
