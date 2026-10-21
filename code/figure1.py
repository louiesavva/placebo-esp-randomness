"""
figure1.py - Figure 1 of "Anomalous Significance Revisited", in black and white.

    python code/figure1.py

Percentage error of the normal and continuity-corrected normal predictions of the
number of experiments with |CR| >= t, relative to the exact binomial (n = 500, p = .5),
for thresholds t = 1.8, 1.9, ..., 5.2. Computed directly from the distributions, so
no data are needed. Writes results/figure1.png (300 dpi) and results/figure1.pdf (vector).
"""
import os
import numpy as np
from scipy.stats import binom, norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
n, p = 500, 0.5
mu, sd = n * p, np.sqrt(n * p * (1 - p))
ts = np.round(np.arange(1.8, 5.2001, 0.1), 1)
normal_err, cc_err = [], []
for t in ts:
    k = int(np.ceil(t * sd - 1e-9))                     # |X - 250| needed
    exact = 2 * binom.sf(mu + k - 1, n, p)               # P(X >= 250+k) + P(X <= 250-k)
    plain = 2 * norm.sf(t)
    cc = 2 * norm.sf((k - 0.5) / sd)
    normal_err.append(100 * (plain / exact - 1))
    cc_err.append(100 * (cc / exact - 1))

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'Liberation Serif', 'DejaVu Serif'],
                     'font.size': 10})
fig, ax = plt.subplots(figsize=(6.5, 3.6))
for lo in (2.5, 3.3, 4.2, 5.0):
    ax.axvspan(lo, lo + 0.1, color='0.88', lw=0)
ax.axhline(0, color='0.5', lw=0.6)
ax.plot(ts, normal_err, '-o', color='black', lw=1.2, ms=3.5, label='Normal approximation')
ax.plot(ts, cc_err, '--o', color='0.35', lw=1.2, ms=3.5, mfc='white', label='Continuity-corrected normal')
ax.set_xlabel('Critical ratio threshold $t$')
ax.set_ylabel('Error relative to exact binomial (%)')
ax.set_xlim(1.75, 5.25)
ax.set_xticks(np.arange(2.0, 5.01, 0.5))
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:+.0f}' if v else '0'))
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, loc='upper left')
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'figure1.png'), dpi=300)
fig.savefig(os.path.join(OUT, 'figure1.pdf'))
print('t=2.5: %.1f%%, t=2.6: %.1f%%' % (normal_err[7], normal_err[8]))
