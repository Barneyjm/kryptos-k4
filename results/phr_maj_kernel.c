/* Majority crib score for interrupted-key models.
   Each model is a partition of 24 crib positions. A class scores its mode
   (the size of the largest set sharing one key value). Singleton positions
   are free and always score 1, passed in as nfree. Score is out of 24.
   A score of 24 is exactly an exact fit (every class is unanimous).

   members: position indices (0..23)
   cls_off: length nclasses+1, slices of members
   model_off: length nmodels+1, slices of classes (size >= 2 only)
   nfree: per model, number of singleton positions */
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int score_model(const uint8_t *key, const uint8_t *members,
                       const int32_t *cls_off, int c0, int c1, int nfree) {
    int sc = nfree;
    for (int c = c0; c < c1; c++) {
        uint8_t cnt[26] = {0};
        int mode = 0;
        for (int j = cls_off[c]; j < cls_off[c + 1]; j++) {
            int v = ++cnt[key[members[j]]];
            if (v > mode) mode = v;
        }
        sc += mode;
    }
    return sc;
}

void maj_partial(const uint8_t *keys, int nk,
                 const uint8_t *members, const int32_t *cls_off, const int32_t *model_off,
                 const uint8_t *nfree, int nmodels,
                 uint64_t *hist, int32_t *best,
                 int32_t *ex, int max_ex, int *nexact) {
    memset(hist, 0, 25 * sizeof(uint64_t));
    best[0] = -1; best[1] = -1; best[2] = -1;
    if (nexact) *nexact = 0;
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        uint64_t h[25] = {0};
        int32_t lb0 = -1, lb1 = -1, lb2 = -1;
        int32_t loc[128];
        int nloc = 0;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int k = 0; k < nk; k++) {
            const uint8_t *kv = keys + 24 * k;
            for (int m = 0; m < nmodels; m++) {
                int sc = score_model(kv, members, cls_off, model_off[m], model_off[m + 1], nfree[m]);
                h[sc]++;
                if (sc > lb0) { lb0 = sc; lb1 = k; lb2 = m; }
                if (sc == 24 && nloc < 64) { loc[2 * nloc] = k; loc[2 * nloc + 1] = m; nloc++; }
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        {
            for (int i = 0; i < 25; i++) hist[i] += h[i];
            if (lb0 > best[0]) { best[0] = lb0; best[1] = lb1; best[2] = lb2; }
            if (ex && nexact) {
                for (int i = 0; i < nloc && *nexact < max_ex; i++) {
                    ex[2 * (*nexact)] = loc[2 * i];
                    ex[2 * (*nexact) + 1] = loc[2 * i + 1];
                    (*nexact)++;
                }
            }
        }
    }
}
