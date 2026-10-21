/* Best partial crib score for an additive clock stream.
   key[i] ≡ sigma * pattern[i] + b  (mod 26), sigma in {+1,-1}.
   Score is the maximum, over b, of the number of the 24 crib positions that agree.
   hist[s] counts (key, pattern, sigma) triples whose score is s.
   best[5] = {score, key, pattern, sigma_code, b} for the first maximum
   (sigma_code 0 => +1, 1 => -1). */
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static void score_one(const uint8_t *k, const uint8_t *p, int *m1, int *b1, int *m2, int *b2) {
    uint8_t c1[26] = {0}, c2[26] = {0};
    int u1 = 0, a1 = 0, u2 = 0, a2 = 0;
    for (int i = 0; i < 24; i++) {
        int d = (int)k[i] - (int)p[i];
        int e = (int)k[i] + (int)p[i];
        if (d < 0) d += 26;
        if (e >= 26) e -= 26;
        int v1 = ++c1[d];
        int v2 = ++c2[e];
        if (v1 > u1) { u1 = v1; a1 = d; }
        if (v2 > u2) { u2 = v2; a2 = e; }
    }
    *m1 = u1; *b1 = a1; *m2 = u2; *b2 = a2;
}

void clock_partial(const uint8_t *keys, int nk, const uint8_t *pat, int np,
                   uint64_t *hist, int32_t *best) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(4);
#endif
    memset(hist, 0, 25 * sizeof(uint64_t));
    best[0] = -1; best[1] = -1; best[2] = -1; best[3] = -1; best[4] = -1;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        uint64_t h[25] = {0};
        int32_t lb[5] = {-1, -1, -1, -1, -1};
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int k = 0; k < nk; k++) {
            const uint8_t *kv = keys + 24 * k;
            for (int s = 0; s < np; s++) {
                int m1, b1, m2, b2;
                score_one(kv, pat + 24 * s, &m1, &b1, &m2, &b2);
                h[m1]++; h[m2]++;
                if (m1 > lb[0]) { lb[0] = m1; lb[1] = k; lb[2] = s; lb[3] = 0; lb[4] = b1; }
                if (m2 > lb[0]) { lb[0] = m2; lb[1] = k; lb[2] = s; lb[3] = 1; lb[4] = b2; }
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            for (int i = 0; i < 25; i++) hist[i] += h[i];
            if (lb[0] > best[0]) memcpy(best, lb, 5 * sizeof(int32_t));
        }
    }
}
