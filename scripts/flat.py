import numpy as np, random, k4, fingerprint as F
from period import pic
random.seed(11)
A = np.array([ord(x)-65 for x in k4.CT])
stat = lambda s: np.mean([pic(np.array([ord(x)-65 for x in s]), L) for L in range(2, 25)])
K = stat(k4.CT)
extra = {
 "Vig KA period 26-50": lambda: F.vig(F.eng(), F.rkey(random.randint(26, 50))),
 "Vig KA period 97 (OTP)": lambda: F.vig(F.eng(), F.rkey(97)),
 "Vig KA key=KA-shuffled blocks": lambda: F.vig(F.eng(), "".join(random.sample(F.AZ, 26)) * 4),
 "Vig key = running key, then 2x columnar": lambda: F.col(F.col(F.vig(F.eng(), F.eng()), random.randint(5, 12)), random.randint(5, 12)),
 "letters dealt evenly (balanced multiset, shuffled)": lambda: "".join(random.sample((F.AZ * 4)[:97], 97)),
}
models = dict(F.MODELS, **extra)
print(f"K4 mean coincidence L=2..24 = {K:.4f}")
for m, f in models.items():
    x = np.array([stat(f()) for _ in range(3000)])
    print(f"{m:48s} mean {x.mean():.4f}  P(<= K4) = {(x <= K).mean():.4f}")
