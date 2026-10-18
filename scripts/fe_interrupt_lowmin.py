"""T1 secondary: weaker interrupted-key models (5-7 checks, e.g. word restarts with unknown first offset) over the
116 clue alphabets only (Quagmire I-IV), where the trial count is small enough for a low threshold. Shuffled calibration."""
import numpy as np, collections, k4, fe_common as F, fe_interrupt as T
CL = F.clue_alphabets()
for MIN in (6,):
    rng = np.random.default_rng(5)
    for label, ct in [("K4", k4.CT)] + [(f"shuffled{t}", "".join(rng.permutation(list(k4.CT)))) for t in range(5)]:
        eng = F.Engine(T.gen_models(ct), MIN)
        h = T.scan(ct, CL, eng) + T.scan(ct, CL, eng, quag4=True)
        nvec = len(CL) * 10 + len(CL) ** 2 * 2
        print(f"MIN={MIN} {label}: {len(eng)} models, {len(h)} hits, chance ~{T.expected_fp(eng, nvec):.2f}", flush=True)
        if label == "K4":
            for x in h[:20]: print("    ", x)
