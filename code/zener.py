"""
zener.py - the random.org Zener-card series: 473,594 placebo experiments of 625 trials
(25 runs of 25), each trial a match between two random symbols out of five (p = 1/5).

    python code/zener.py

Writes results/zener.md. Needs numpy and scipy; run after extract.py.

With n = 625 and p = 1/5 the standard deviation is exactly 10, so CR = (X - 125)/10 moves
in steps of exactly 0.1. Thresholds on a 0.1 grid therefore always fall on a possible
score, and the plain normal approximation always leaves out about half of the probability
of the boundary score: unlike the 500-trial design, its error does not change sign.
"""
import os
import numpy as np
from scipy.stats import binom, chi2, norm, beta

HERE = os.path.dirname(os.path.abspath(__file__))
X = np.load(os.path.join(HERE, '..', 'data', 'random_org_zener.npz'))['hits'].astype(np.int64)
OUT = os.path.join(HERE, '..', 'results', 'zener.md')
N, n, p = len(X), 625, 0.2
mu, sd = n * p, np.sqrt(n * p * (1 - p))
pmf = binom.pmf(np.arange(n + 1), n, p)


def cp(k, N, level=0.95):
    """Clopper-Pearson interval for a proportion."""
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else beta.ppf(a, k, N - k + 1)
    hi = 1.0 if k == N else beta.ppf(1 - a, k + 1, N - k)
    return lo, hi


md = ['# random.org Zener-card series', '',
      f'{N:,} experiments of {n} trials, hit probability 1/5 (two random symbols out of five match). '
      f'CR = (X - {mu:g})/{sd:g}.', '']

# ---- moments and fit -------------------------------------------------------------
d = X - mu
var = X.var()
z_var = (var / (n * p * (1 - p)) - 1) / np.sqrt(2 / N)
skew = (d ** 3).mean() / var ** 1.5
kurt = (d ** 4).mean() / var ** 2 - 3
skew_b = (1 - 2 * p) / sd
kurt_b = (1 - 6 * p * (1 - p)) / (n * p * (1 - p))
c = np.bincount(X, minlength=n + 1)[:n + 1].astype(float)
e = N * pmf
keep = e >= 5; i0 = np.argmax(keep); i1 = len(keep) - np.argmax(keep[::-1])
O = np.r_[c[:i0].sum(), c[i0:i1], c[i1:].sum()]; E = np.r_[e[:i0].sum(), e[i0:i1], e[i1:].sum()]
X2 = float(((O - E) ** 2 / E).sum()); df = len(O) - 1
md += ['## Moments and fit', '',
       '| statistic | observed | binomial | z |', '|---|---|---|---|',
       f'| mean | {X.mean():.4f} | {mu:g} | {(X.mean() - mu) / (sd / np.sqrt(N)):+.2f} |',
       f'| variance | {var:.4f} | {n * p * (1 - p):g} | {z_var:+.2f} |',
       f'| skewness | {skew:.4f} | {skew_b:.4f} | {(skew - skew_b) / np.sqrt(6 / N):+.2f} |',
       f'| excess kurtosis | {kurt:.5f} | {kurt_b:.5f} | {(kurt - kurt_b) / np.sqrt(24 / N):+.2f} |', '',
       f'Chi-square goodness of fit to Binomial({n}, 1/5): {X2:.1f} on {df} df, p = {chi2.sf(X2, df):.3f} '
       '(cells with expected count below 5 pooled into the two tails).', '']

# ---- two-sided counts ------------------------------------------------------------
md += ['## Experiments with |CR| >= t', '',
       '| t | needs abs(X - 125) >= | observed | exact | normal | normal, corrected | obs/exact | obs/normal | exact/normal | (obs - exact)/SE |',
       '|---|---|---|---|---|---|---|---|---|---|']
for t in [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5]:
    k = int(np.ceil(t * sd - 1e-9))
    obs = int((np.abs(d) >= k).sum())
    q = binom.sf(mu + k - 1, n, p) + binom.cdf(mu - k, n, p)
    ex, nm, cc = N * q, N * 2 * norm.sf(t), N * 2 * norm.sf((k - 0.5) / sd)
    md.append(f'| {t:.1f} | {k} | {obs:,} | {ex:,.1f} | {nm:,.1f} | {cc:,.1f} | {obs / ex:.3f} | {obs / nm:.3f} | '
              f'{ex / nm:.3f} | {(obs - ex) / np.sqrt(ex * (1 - q)):+.2f} |')
md += ['', 'The exact count exceeds the normal prediction at every threshold, because each threshold '
       'coincides with a possible score.', '']

# ---- one-sided tails ---------------------------------------------------------------
md += ['## Upper and lower tails separately', '',
       '| t | upper: observed | exact | normal | lower: observed | exact | normal |', '|---|---|---|---|---|---|---|']
for t in (2.0, 2.5, 3.0, 3.5):
    k = int(np.ceil(t * sd - 1e-9))
    up, lo = int((d >= k).sum()), int((d <= -k).sum())
    md.append(f'| {t:.1f} | {up:,} | {N * binom.sf(mu + k - 1, n, p):,.1f} | {N * norm.sf(t):,.1f} | '
              f'{lo:,} | {N * binom.cdf(mu - k, n, p):,.1f} | {N * norm.sf(t):,.1f} |')
md += ['', 'With p = 1/5 the binomial is skewed to the right. Together with the lattice effect, this makes the '
       'upper tail (more hits than expected, the direction an ESP test looks for) much heavier than the normal '
       'curve, while from t = 2.5 upwards the lower tail is lighter.', '']

# ---- true size of nominal one-sided tests (as in Table 5 of the paper) ---------------
md += ['## True size of nominal one-sided z tests, 625 trials', '',
       '| nominal | smallest qualifying X | exact size | observed proportion (95% interval) |', '|---|---|---|---|']
for z, lab in ((1.645, '5%'), (1.96, '2.5%'), (2.326, '1%')):
    k = int(np.ceil(mu + z * sd - 1e-9))
    obs = int((X >= k).sum()); lo, hi = cp(obs, N)
    md.append(f'| z >= {z} ({lab}) | {k} | {binom.sf(k - 1, n, p):.2%} | {obs / N:.2%} ({lo:.2%} to {hi:.2%}) |')
md += ['']

with open(OUT, 'w') as fh:
    fh.write('\n'.join(md) + '\n')
print('\n'.join(md))
