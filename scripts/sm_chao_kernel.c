/* Chaocipher scoring kernel (same algorithm as sm_chao.py).
   For every (clue a, dict d, assign, off, dir): assign 0: L=clue[a], R=rot(dict[d],off); assign 1: L=dict[d], R=rot(clue[a],off).
   dir 0: ct read from left disk; dir 1: ct read from right disk. Decrypt positions 0..73, score 24 cribs. */
#include <stdint.h>
#include <string.h>
static const int SIGL[26] = {0,2,3,4,5,6,7,8,9,10,11,12,13,1,14,15,16,17,18,19,20,21,22,23,24,25};
static const int SIGR[26] = {1,2,4,5,6,7,8,9,10,11,12,13,14,3,15,16,17,18,19,20,21,22,23,24,25,26};
static inline void score_one(const uint8_t *L0, const uint8_t *R0, int sw, const uint8_t *ct, const uint8_t *need,
                             int *tot, int *b1) {
    uint8_t L[52], R[52], T[26];
    memcpy(L, L0, 26); memcpy(R, R0, 26);
    int t_ = 0, b_ = 0;
    for (int t = 0; t < 74; t++) {
        uint8_t c = ct[t]; int i = 0;
        uint8_t *A = sw ? R : L, *B = sw ? L : R;
        while (A[i] != c) i++;
        uint8_t p = B[i];
        if (need[t] != 255) { if (p == need[t]) { t_++; if (t < 34) b_++; } }
        memcpy(L + 26, L, 26); memcpy(R + 26, R, 26);
        for (int j = 0; j < 26; j++) T[j] = L[i + SIGL[j]];
        memcpy(L, T, 26);
        for (int j = 0; j < 26; j++) T[j] = R[i + SIGR[j]];
        memcpy(R, T, 26);
    }
    *tot = t_; *b1 = b_;
}
/* returns number of records written (capped at maxrec). rec layout: a, d, assign, off, dir, tot, b1 (int32 each) */
int sweep(const uint8_t *clue, int na, const uint8_t *dict, int nd, const uint8_t *ct, const uint8_t *need,
          int thr, int64_t *hist_tot, int64_t *hist_b1, int32_t *rec, int maxrec) {
    int nr = 0; uint8_t Rr[26];
    for (int d = 0; d < nd; d++) for (int a = 0; a < na; a++) for (int as = 0; as < 2; as++)
    for (int off = 0; off < 26; off++) {
        const uint8_t *Ls = as ? dict + 26 * d : clue + 26 * a;
        const uint8_t *Rs = as ? clue + 26 * a : dict + 26 * d;
        for (int j = 0; j < 26; j++) Rr[j] = Rs[(j + off) % 26];
        for (int sw = 0; sw < 2; sw++) {
            int tot, b1; score_one(Ls, Rr, sw, ct, need, &tot, &b1);
            hist_tot[tot]++; hist_b1[b1]++;
            if (tot >= thr && nr < maxrec) {
                int32_t *r = rec + 7 * nr++; r[0]=a; r[1]=d; r[2]=as; r[3]=off; r[4]=sw; r[5]=tot; r[6]=b1;
            }
        }
    }
    return nr;
}
