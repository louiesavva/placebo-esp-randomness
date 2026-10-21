"""
Iversen, Longcor, Mosteller, Gilbert & Youtz (1971), "Bias and runs in dice
throwing and recording: A few million throws", Psychometrika 36, 1-19.

Rebuilds their battery of 67 standard scores per data set and runs it on
perfectly random (simulated) odd/even throws, to see how often the battery
flags |z| > 2 under ideal randomness, and how unusual their published
counts of extreme scores were.

Requires: Python 3 and numpy  (pip install numpy)
Run:      python iversen_check.py              (about 2-4 minutes)
          python iversen_check.py 20000 10000  (as reported in the paper; about 25 minutes)
"""
import sys, time
import numpy as np

N_RAND = int(sys.argv[1]) if len(sys.argv) > 1 else 1000   # 100-block data sets
N_DICE = int(sys.argv[2]) if len(sys.argv) > 2 else 500    # dice sets (100+30+31+58 blocks)
SEED = 2026
rng = np.random.default_rng(SEED)

n = 20000          # throws per block
p1 = p2 = 0.5      # P(even), P(odd)
L = 10             # run lengths 1..10

# ---- Theoretical moments: Iversen et al. Appendix 1 (after Mood, 1940) ----
def ff2(x):                      # factorial power x^(2) = x(x-1)
    return x * (x - 1)

def E(i):                        # mean number of even runs of length exactly i
    return p1**i * p2 * ((n - i - 1) * p2 + 2)

def C(i, j):                     # variance (i == j) or covariance
    if i == j:
        m = n - 2*i - 2
        return p1**(2*i) * p2**2 * (ff2(m)*p2**2 + 6*m*p2 + 6) - E(i)**2 + E(i)
    m = n - i - j - 2
    return p1**(i+j) * p2**2 * (ff2(m)*p2**2 + 6*m*p2 + 6) - E(i)*E(j)

Ev  = np.array([E(i) for i in range(1, L+1)])
Sig = np.array([[C(i, j) for j in range(1, L+1)] for i in range(1, L+1)])
K = 10                           # "long runs" = length 10 or more (their Table 2)
E_long = p1**K * ((n - K) * p2 + 1)

# ---- Simulate blocks of fair odd/even throws -------------------------------
def blocks(nb):
    """Return (evens per block, counts of even runs of length 1..10, long runs >=10)."""
    out_e, out_r, out_l = [], [], []
    for start in range(0, nb, 500):
        m = min(500, nb - start)
        x = rng.integers(0, 2, size=(m, n), dtype=np.int8)        # 1 = even
        pad = np.zeros((m, n + 2), dtype=np.int8); pad[:, 1:-1] = x
        d = np.diff(pad, axis=1)
        rs, cs = np.nonzero(d == 1)                               # run starts
        _,  ce = np.nonzero(d == -1)                              # run ends
        cat = np.minimum(ce - cs, 11)                             # 1..10, 11 = longer
        cnt = np.bincount(rs*12 + cat, minlength=m*12).reshape(m, 12)
        out_e.append(x.sum(1)); out_r.append(cnt[:, 1:11]); out_l.append(cnt[:, 10] + cnt[:, 11])
    return (np.concatenate(out_e).astype(float),
            np.concatenate(out_r).astype(float),
            np.concatenate(out_l).astype(float))

# ---- The 67 standard scores ------------------------------------------------
def battery(evens, runs, longr):
    b = len(evens); z = []
    z.append((evens.mean() - n/2) / (evens.std(ddof=1)/np.sqrt(b)))       # Table 1
    z.append((longr.mean() - E_long) / (longr.std(ddof=1)/np.sqrt(b)))    # Table 2
    sd = np.sqrt(np.diag(Sig))
    z += list((runs.mean(0) - Ev) / (sd/np.sqrt(b)))                        # Table 5
    s2 = runs.var(0, ddof=1)
    z += list(((b-1)*s2/np.diag(Sig) - (b-1)) / np.sqrt(2*(b-1)))           # Table 6
    S = np.cov(runs.T)
    for i in range(L):                                                      # Table 7
        for j in range(i+1, L):
            z.append(np.sqrt(b-1)*(S[i, j]-Sig[i, j]) / np.sqrt(Sig[i, j]**2 + Sig[i, i]*Sig[j, j]))
    return np.array(z)

def run_battery(b):
    return battery(*blocks(b))

# ---- 1. Check the formulas against simulation ------------------------------
print("1. Checking the theoretical formulas against 5,000 simulated blocks")
e, r, lg = blocks(5000)
print("   mean runs, length 1-3   theory", np.round(Ev[:3], 1), " simulated", np.round(r.mean(0)[:3], 1))
print("   variance, length 1-3    theory", np.round(np.diag(Sig)[:3], 1), " simulated", np.round(r.var(0)[:3], 1))
print("   long runs (>=10) mean   theory", round(E_long, 2), " simulated", round(lg.mean(), 2), " (their Table 2: 9.76)")

# ---- 2. Calibration and the RAND / pseudo-random comparison -----------------
print(f"\n2. Running {N_RAND} simulated 100-block data sets (like their RAND and pseudo-random series)")
t0 = time.time(); Z = []
for k in range(N_RAND):
    Z.append(run_battery(100))
    if (k+1) % 1000 == 0: print(f"   {k+1}/{N_RAND}  ({time.time()-t0:.0f}s)")
Z = np.array(Z)
ext = (np.abs(Z) > 2).sum(1)
cov_var = Z[:, 22:67].var(1)
nominal = 0.0455
print(f"   Share of scores beyond +/-2: {(np.abs(Z) > 2).mean():.4f}   (nominal {nominal})")
print(f"   Covariance scores only:      {(np.abs(Z[:, 22:]) > 2).mean():.4f}")
print(f"   Mean extremes per 67 scores: {ext.mean():.2f}   (nominal {67*nominal:.2f})")
def prop(hits, n):
    """Proportion with an exact 95% (Clopper-Pearson) interval."""
    from scipy.stats import beta
    lo = beta.ppf(0.025, hits, n - hits + 1) if hits else 0.0
    hi = beta.ppf(0.975, hits + 1, n - hits) if hits < n else 1.0
    return f"{hits/n:.4f}  ({hits} of {n}; 95% CI {lo:.4f} to {hi:.4f})"
print(f"   P(>= 5 extremes)  [their pseudo-random had 5]: {prop(int((ext >= 5).sum()), len(ext))}")
print(f"   P(>= 10 extremes) [their RAND had 10]:         {prop(int((ext >= 10).sum()), len(ext))}")
print(f"   P(covariance-score variance >= 1.97) [their RAND]: {prop(int((cov_var >= 1.97).sum()), len(ext))}")

# ---- 3. The dice: brands A, B, C, X with their real block counts -----------
print(f"\n3. Running {N_DICE} simulated dice sets (100 + 30 + 31 + 58 blocks)")
# Iversen et al. counted 267 dice scores, not 4 x 67 = 268: their summary covers "all the
# standard scores that we show in this publication", one fewer than the full battery. So we
# also count extremes with one score, chosen at random, left out of each simulated set.
t0 = time.time(); D = []; D267 = []
for k in range(N_DICE):
    flags = np.concatenate([np.abs(run_battery(b)) > 2 for b in (100, 30, 31, 58)])
    D.append(int(flags.sum()))
    D267.append(int(np.delete(flags, rng.integers(len(flags))).sum()))
    if (k+1) % 500 == 0: print(f"   {k+1}/{N_DICE}  ({time.time()-t0:.0f}s)")
D = np.array(D); D267 = np.array(D267)
print(f"   Mean extremes per 268 scores: {D.mean():.2f}; per 267 scores: {D267.mean():.2f}")
print(f"   P(>= 16 extremes of 268): {prop(int((D >= 16).sum()), len(D))}")
print(f"   P(>= 16 extremes of 267) [their dice had 16 of 267]: {prop(int((D267 >= 16).sum()), len(D267))}")

print("\nSummary: compare with Iversen et al.'s ratios of 1.3 (dice), 3 (RAND) and 1.7 (pseudo-random).")
