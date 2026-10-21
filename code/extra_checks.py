"""
extra_checks.py - follow-up checks on the random.org main series, reported in the
section "Does the standard error vary?" and in the description of Table 1.

    python code/extra_checks.py

Writes results/extra_checks.md. Needs numpy and scipy; run after extract.py.

1. Fifteen thresholds. Counts of |CR| >= t at t = 2.0, 2.1, ..., 3.0, 3.5, 4.0, 4.5, 5.0
   overlap, so the chance of the largest |z| among them reaching a given size is found by
   simulating series that follow the binomial exactly (multinomial draws, 20,000 series).
2. Kurtosis and central peak of the difference M - A and the total M + A. A standard
   error that varied from one experiment to the next (a scale mixture) would raise the
   kurtosis of both scores alike, by 3 x (relative variance of the variance).
3. Upper and lower tails of the difference separately.
4. The +6 displacement fit of the random.org digit series, with its Bonferroni threshold.
"""
import os
import numpy as np
from scipy.stats import binom, chi2, norm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'results')
load = lambda f: np.load(os.path.join(DATA, f))

m = load('random_org_main_morning.npz')['morning'].astype(np.int64)
a = load('random_org_main_afternoon.npz')['afternoon'].astype(np.int64)
N = len(m)
D = m - a                 # difference, mean 0
T = m + a - 250           # total, centred
pm = binom.pmf(np.arange(501), 500, 0.5)
md = ['# Extra checks on the random.org main series', '', f'{N:,} experiments.', '']

# ---- 1. fifteen thresholds --------------------------------------------------
ts = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0]
lims = [int(np.ceil(t * np.sqrt(125) - 1e-9)) for t in ts]      # |M - A| needed
absp = np.r_[pm[250], pm[251:] + pm[249::-1]]                   # P(|D| = k), k = 0..250
cnt_abs = np.bincount(np.abs(D), minlength=251)
obs = np.array([cnt_abs[L:].sum() for L in lims])
exp = np.array([N * absp[L:].sum() for L in lims])
z = (obs - exp) / np.sqrt(exp * (1 - exp / N))
edges = [0] + lims + [251]
cells = np.array([absp[edges[i]:edges[i + 1]].sum() for i in range(len(edges) - 1)])
rng = np.random.default_rng(2026)
S = rng.multinomial(N, cells / cells.sum(), size=20000)
sim_tails = np.cumsum(S[:, ::-1], axis=1)[:, ::-1][:, 1:]
sim_z = (sim_tails - exp) / np.sqrt(exp * (1 - exp / N))
zmax = np.abs(z).max()
p15 = (np.abs(sim_z).max(1) >= zmax - 1e-9).mean()
md += ['## 1. Fifteen thresholds', '', '| t | needs abs(M - A) >= | observed | exact | obs/exact | z |', '|---|---|---|---|---|---|']
md += [f'| {t} | {L} | {o:,} | {e:,.1f} | {o / e:.3f} | {zz:+.2f} |' for t, L, o, e, zz in zip(ts, lims, obs, exp, z)]
md += ['', f'Largest |z| = {zmax:.2f}. In 20,000 simulated binomial series, a largest |z| at least this big '
       f'occurred in {p15:.1%} of them.', '']

# ---- 2. kurtosis and peak ----------------------------------------------------
se_k = np.sqrt(24 / N)


def kurt(x):
    v = x.var()
    return v, (x.astype(float) ** 4).mean() / v ** 2 - 3


md += ['## 2. Kurtosis and central peak', '',
       f'Standard error of excess kurtosis = sqrt(24/N) = {se_k:.6f}.', '',
       '| score | variance (expected) | excess kurtosis | binomial value | z |', '|---|---|---|---|---|']
rows = {}
for name, x, n, c in (('difference M - A', D, 500, 0), ('total M + A', T, 500, 0),
                      ('morning', m - 125, 250, 0), ('afternoon', a - 125, 250, 0)):
    v, k = kurt(x)
    zk = (k + 2 / n) / se_k
    rows[name] = (k, zk)
    md.append(f'| {name} | {v:.4f} ({n / 4:g}) | {k:+.5f} | {-2 / n:+.5f} | {zk:+.2f} |')
kT = rows['total M + A'][0] + 2 / 500
upper = kT + 1.96 * se_k
md += ['', 'If the standard error varied from experiment to experiment, with relative variance of the variance v, '
       'the excess kurtosis of both scores would rise by about 3v.',
       f'From the total score, 3v <= {upper:.5f} (upper 95% limit), so v <= {upper / 3:.6f}: the variance '
       f'varies by at most {np.sqrt(upper / 3):.1%} and the standard error by at most about {np.sqrt(upper / 3) / 2:.1%}.', '']
md += ['| score | cells | observed | expected | z |', '|---|---|---|---|---|']
for name, x in (('difference M - A', D), ('total M + A', T)):
    for k in (0, 1, 2):
        o = int((np.abs(x) <= k).sum())
        e = N * pm[250 - k:251 + k].sum()
        md.append(f'| {name} | abs <= {k} | {o:,} | {e:,.0f} | {(o - e) / np.sqrt(e * (1 - e / N)):+.2f} |')
z1 = (int((np.abs(D) <= 1).sum()) - N * pm[249:252].sum()) / np.sqrt(N * pm[249:252].sum() * (1 - pm[249:252].sum()))
md += ['', f'Nominal two-sided p for the excess within one hit of M = A (z = {z1:.2f}): {2 * norm.sf(abs(z1)):.4f}. '
       'This window was chosen after looking at the data, among overlapping comparisons, so the nominal p overstates the evidence.', '']

# ---- 3. tails separately ------------------------------------------------------
md += ['## 3. Upper and lower tails of the difference', '',
       '| t | morning ahead (M - A >= L) | afternoon ahead (A - M >= L) | expected each |', '|---|---|---|---|']
for t, L in ((3.0, 34), (3.5, 40), (4.0, 45)):
    up = int((D >= L).sum()); lo = int((D <= -L).sum()); e = N * pm[250 + L:].sum()
    md.append(f'| {t} | {up:,} | {lo:,} | {e:,.0f} |')
L = 40; lo = int((D <= -L).sum()); e = N * pm[250 + L:].sum(); zlo = (lo - e) / np.sqrt(e * (1 - e / N))
md += ['', f'Afternoon-ahead count at t = 3.5: z = {zlo:.2f}, nominal two-sided p = {2 * norm.sf(abs(zlo)):.4f} (chosen after looking at the data).', '']

# ---- 4. the +6 displacement ---------------------------------------------------
x = load('random_org_digits.npz')['plus6'].astype(np.int64)
n, p = 268, 0.1
c = np.bincount(x, minlength=n + 1)[:n + 1].astype(float)
e = len(x) * binom.pmf(np.arange(n + 1), n, p)
keep = e >= 5; i0 = np.argmax(keep); i1 = len(keep) - np.argmax(keep[::-1])
O = np.r_[c[:i0].sum(), c[i0:i1], c[i1:].sum()]; E = np.r_[e[:i0].sum(), e[i0:i1], e[i1:].sum()]
X2 = float(((O - E) ** 2 / E).sum()); df = len(O) - 1; pv = chi2.sf(X2, df)
half = len(x) // 2
ph = []
for part in (x[:half], x[half:]):
    c = np.bincount(part, minlength=n + 1)[:n + 1].astype(float)
    e = len(part) * binom.pmf(np.arange(n + 1), n, p)
    keep = e >= 5; i0 = np.argmax(keep); i1 = len(keep) - np.argmax(keep[::-1])
    O = np.r_[c[:i0].sum(), c[i0:i1], c[i1:].sum()]; E = np.r_[e[:i0].sum(), e[i0:i1], e[i1:].sum()]
    ph.append(chi2.sf(((O - E) ** 2 / E).sum(), len(O) - 1))
md += ['## 4. random.org digits, +6 displacement', '',
       f'chi-square {X2:.1f} on {df} df, p = {pv:.5f}; Bonferroni threshold for 35 fits = {0.05 / 35:.5f}.',
       f'Mean {x.mean():.4f} (expected {n * p}), variance {x.var():.4f} (expected {n * p * (1 - p):.1f}).',
       f'First and second halves of the series: p = {ph[0]:.3f} and {ph[1]:.3f}.', '']

with open(os.path.join(OUT, 'extra_checks.md'), 'w') as fh:
    fh.write('\n'.join(md) + '\n')
print('\n'.join(md))
