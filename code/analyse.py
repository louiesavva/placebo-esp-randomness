"""
analyse.py - every number in the papers, computed from data/.

Usage:   python code/analyse.py          (about a minute)
Writes:  results/*.csv and results/summary.md
Needs:   numpy, scipy
"""
import os
import numpy as np
from scipy.stats import binom, norm, chi2

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, '..', 'data')
RES = os.path.join(HERE, '..', 'results')
os.makedirs(RES, exist_ok=True)
md = []                                    # lines of results/summary.md


def load(name):
    return dict(np.load(os.path.join(DATA, name)))


def write_csv(name, header, rows):
    with open(os.path.join(RES, name), 'w') as fh:
        fh.write(','.join(header) + '\n')
        for r in rows:
            fh.write(','.join(str(v) for v in r) + '\n')


def tail_counts(x, n, p, t):
    """Two-sided count of |CR| >= t for hit counts x ~ Bin(n, p), with expectations."""
    mu, sd = n * p, np.sqrt(n * p * (1 - p))
    hi = int(np.ceil(mu + t * sd - 1e-9))          # smallest upper count reaching t
    lo = int(np.floor(mu - t * sd + 1e-9))         # largest lower count reaching t
    N = len(x)
    obs = int(((x >= hi) | (x <= lo)).sum())
    p_exact = binom.sf(hi - 1, n, p) + (binom.cdf(lo, n, p) if lo >= 0 else 0.0)
    p_norm = 2 * norm.sf(t)
    p_cc = norm.sf((hi - 0.5 - mu) / sd) + (norm.cdf((lo + 0.5 - mu) / sd) if lo >= 0 else 0.0)
    return dict(obs=obs, exact=N * p_exact, normal=N * p_norm, corrected=N * p_cc,
                p_exact=p_exact, hi=hi, lo=lo, z=(obs - N * p_exact) / np.sqrt(N * p_exact * (1 - p_exact)))


def gof(x, n, p):
    """Chi-square goodness of fit of hit counts to Binomial(n, p), pooling cells with expectation < 5."""
    N = len(x)
    cnt = np.bincount(x, minlength=n + 1)[:n + 1].astype(float)
    e = N * binom.pmf(np.arange(n + 1), n, p)
    keep = e >= 5
    i0 = np.argmax(keep); i1 = len(keep) - np.argmax(keep[::-1])
    O = np.r_[cnt[:i0].sum(), cnt[i0:i1], cnt[i1:].sum()]
    E = np.r_[e[:i0].sum(), e[i0:i1], e[i1:].sum()]
    m = E > 0
    X2 = float(((O[m] - E[m]) ** 2 / E[m]).sum()); df = int(m.sum() - 1)
    return X2, df, float(chi2.sf(X2, df))


# =============================================================================
# 1. random.org main series: Spencer-Brown's morning-afternoon difference
# =============================================================================
m = load('random_org_main_morning.npz')['morning'].astype(np.int64)
a = load('random_org_main_afternoon.npz')['afternoon'].astype(np.int64)
N = len(m)
D = m - a + 250            # morning minus afternoon, shifted: exactly Binomial(500, 1/2) under chance
T = m + a                  # total hits: also Binomial(500, 1/2)

md += ['# Results', '',
       f'## 1. random.org main series ({N:,} experiments of 250 morning + 250 afternoon trials)', '',
       'Critical ratio = (morning hits - afternoon hits) / sqrt(125). Under chance, morning - afternoon + 250 ~ Binomial(500, 1/2) exactly.', '',
       '| t | observed | exact binomial | normal | continuity-corrected | (obs-exact)/SE |', '|---|---|---|---|---|---|']
rows = []
for t in [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0]:
    r = tail_counts(D, 500, 0.5, t)
    rows.append([t, r['obs'], round(r['exact'], 1), round(r['normal'], 1), round(r['corrected'], 1), round(r['z'], 2)])
    if t in (2.0, 2.2, 2.5, 2.6, 3.0, 4.0, 5.0):
        md.append(f"| {t} | {r['obs']:,} | {r['exact']:,.0f} | {r['normal']:,.0f} | {r['corrected']:,.0f} | {r['z']:+.2f} |")
write_csv('table1_random_org_difference.csv', ['t', 'observed', 'exact', 'normal', 'corrected', 'z_vs_exact'], rows)

rows = []
for t in [2.0, 2.2, 2.5, 2.6, 3.0, 4.0, 5.0]:
    r = tail_counts(T, 500, 0.5, t)
    rows.append([t, r['obs'], round(r['exact'], 1), round(r['normal'], 1), round(r['z'], 2)])
write_csv('random_org_total_score.csv', ['t', 'observed', 'exact', 'normal', 'z_vs_exact'], rows)

# error of the approximations by threshold (Figure data)
rows = []
for t in np.round(np.arange(1.8, 5.21, 0.1), 1):
    d = int(np.ceil(t * np.sqrt(125) - 1e-9))
    pe = 2 * binom.sf(250 + d - 1, 500, 0.5)
    rows.append([t, round(100 * (2 * norm.sf(t) / pe - 1), 2), round(100 * (2 * norm.sf((d - 0.5) / np.sqrt(125)) / pe - 1), 2)])
write_csv('figure_error_by_threshold.csv', ['t', 'normal_error_pct', 'corrected_error_pct'], rows)

# nominal 1.96 test at n = 500
d196 = int(np.ceil(1.96 * np.sqrt(125) - 1e-9))
size196 = 2 * binom.sf(250 + d196 - 1, 500, 0.5)
md += ['', f'Nominal |CR| >= 1.96 test: requires |deviation| >= {d196}; exact size {size196:.4f}.']

# far tail (Gilmore's "seventh decimal place")
md += ['', '| threshold | exact two-sided probability | observed / expected | 95% range if exact |', '|---|---|---|---|']
rows = []
for t in (3.0, 4.0, 5.0):
    r = tail_counts(D, 500, 0.5, t)
    from scipy.stats import poisson
    lo95, hi95 = poisson.ppf([0.025, 0.975], r['exact'])
    rows.append([t, f"{r['p_exact']:.3g}", round(r['obs'] / r['exact'], 3), round(lo95 / r['exact'], 2), round(hi95 / r['exact'], 2)])
    md.append(f"| {t} | {r['p_exact']:.2g} | {r['obs'] / r['exact']:.3f} | {lo95 / r['exact']:.2f}-{hi95 / r['exact']:.2f} |")
write_csv('far_tail.csv', ['t', 'p_exact', 'obs_over_exact', 'range95_low', 'range95_high'], rows)

# =============================================================================
# 2. Gilmore's "standard error of the standard error"
# =============================================================================
x = (D - 250).astype(float)
mu, var = x.mean(), x.var()
exk = ((x - mu) ** 4).mean() / var ** 2 - 3
se_k = np.sqrt(24 / N)
p0 = binom.pmf(250, 500, 0.5)
mode_obs = int((D == 250).sum()); mode_exp = N * p0
X2, df, pg = gof(D, 500, 0.5)
days = N // 16777
Kd = x.reshape(days, 16777)
dv = Kd.var(axis=1, ddof=1)
mu4 = 125 ** 2 * (3 - 2 / 500)
var_s2 = (mu4 - 125 ** 2 * (16777 - 3) / (16777 - 1)) / 16777
Q = float(((dv - 125) ** 2 / var_s2).sum())
dm = Kd.mean(axis=1); Qm = float(((dm / np.sqrt(125 / 16777)) ** 2).sum())
cma = np.corrcoef(m, a)[0, 1]; cma2 = np.corrcoef((m - 125.0) ** 2, (a - 125.0) ** 2)[0, 1]
kur_m = (((m - m.mean()) ** 4).mean() / m.var() ** 2 - 3); kur_a = (((a - a.mean()) ** 4).mean() / a.var() ** 2 - 3)
g = [('mean of difference', f'{mu:+.5f}', '0', f'{mu / np.sqrt(125 / N):+.2f}'),
     ('variance of difference', f'{var:.4f}', '125', f'{(var - 125) / (125 * np.sqrt((2 - 2 / 500) / N)):+.2f}'),
     ('excess kurtosis of difference', f'{exk:+.5f}', f'{-2 / 500:+.4f} (binomial); 0 (normal)', f'{(exk + 2 / 500) / se_k:+.2f} vs binomial; {exk / se_k:+.2f} vs normal'),
     ('excess kurtosis, morning', f'{kur_m:+.5f}', f'{-2 / 250:+.4f}', f'{(kur_m + 2 / 250) / se_k:+.2f}'),
     ('excess kurtosis, afternoon', f'{kur_a:+.5f}', f'{-2 / 250:+.4f}', f'{(kur_a + 2 / 250) / se_k:+.2f}'),
     ('experiments scoring exactly 0', f'{mode_obs:,}', f'{mode_exp:,.0f}', f'{(mode_obs - mode_exp) / np.sqrt(mode_exp * (1 - p0)):+.2f}'),
     ('chi-square, whole distribution', f'{X2:.1f} on {df} df', '', f'p = {pg:.3f}'),
     (f'spread of daily variances ({days} days)', f'SD {dv.std(ddof=1):.3f}', f'SD {np.sqrt(var_s2):.3f}', f'dispersion {Q:.0f} on {days} df, p = {chi2.sf(Q, days):.2f}'),
     ('spread of daily means', f'{Qm:.0f}', f'{days} df', f'p = {chi2.sf(Qm, days):.2f}'),
     ('corr(morning, afternoon)', f'{cma:+.6f}', '0', f'{cma * np.sqrt(N):+.2f}'),
     ('corr(morning^2, afternoon^2)', f'{cma2:+.6f}', '0', f'{cma2 * np.sqrt(N):+.2f}')]
write_csv('gilmore_check.csv', ['quantity', 'observed', 'expected', 'test'], [[q, f'"{o}"', f'"{e}"', f'"{z}"'] for q, o, e, z in g])
md += ['', '## 2. Gilmore: does the standard error itself vary?', '', '| quantity | observed | expected | test |', '|---|---|---|---|']
md += [f'| {q} | {o} | {e} | {z} |' for q, o, e, z in g]

# =============================================================================
# 3. Every series against the exact binomial
# =============================================================================
series = []
series += [('random.org (Example 1B)', 'morning-afternoon', 'difference', D, 500, 0.5),
           ('random.org (Example 1B)', 'morning-afternoon', 'total', T, 500, 0.5),
           ('random.org (Example 1B)', 'morning-afternoon', 'morning', m, 250, 0.5),
           ('random.org (Example 1B)', 'morning-afternoon', 'afternoon', a, 250, 0.5)]
pr = load('prng.npz'); pm, pa = pr['morning'].astype(np.int64), pr['afternoon'].astype(np.int64)
series += [('R Mersenne Twister (Example 1A)', 'morning-afternoon', 'difference', pm - pa + 250, 500, 0.5),
           ('R Mersenne Twister (Example 1A)', 'morning-afternoon', 'total', pm + pa, 500, 0.5),
           ('R Mersenne Twister (Example 1A)', 'morning-afternoon', 'morning', pm, 250, 0.5),
           ('R Mersenne Twister (Example 1A)', 'morning-afternoon', 'afternoon', pa, 250, 0.5)]
ek = load('entropy_key_digits.npz')
series += [('Entropy Key (Example 2)', 'digit matching', 'straight', ek['straight'].astype(np.int64), 500, 0.1),
           ('Entropy Key (Example 2)', 'digit matching', '+1 displacement', ek['plus1'].astype(np.int64), 490, 0.1),
           ('Entropy Key (Example 2)', 'digit matching', '-1 displacement', ek['minus1'].astype(np.int64), 490, 0.1)]
rd = load('random_org_digits.npz')
series += [('random.org (Example 3)', 'digit matching', 'straight', rd['straight'].astype(np.int64), 268, 0.1)]
for k in range(1, 8):
    series += [('random.org (Example 3)', 'digit matching', f'+{k} displacement', rd[f'plus{k}'].astype(np.int64), 268, 0.1),
               ('random.org (Example 3)', 'digit matching', f'-{k} displacement', rd[f'minus{k}'].astype(np.int64), 268, 0.1)]
r125 = load('random_org_125.npz'); a1, a2 = r125['morning'].astype(np.int64), r125['afternoon'].astype(np.int64)
series += [('random.org, 125 + 125', 'morning-afternoon', 'total', a1 + a2, 250, 0.5),
           ('random.org, 125 + 125', 'morning-afternoon', 'morning', a1, 125, 0.5),
           ('random.org, 125 + 125', 'morning-afternoon', 'afternoon', a2, 125, 0.5)]
ot = load('other_runs.npz')
for key, lab in (('mcr', 'mcr.RData run (source not recorded)'), ('aao1rng', 'aao1rng run (source not recorded)')):
    o1, o2 = ot[f'{key}_morning'].astype(np.int64), ot[f'{key}_afternoon'].astype(np.int64)
    series += [(lab, 'morning-afternoon', 'difference', o1 - o2 + 250, 500, 0.5),
               (lab, 'morning-afternoon', 'morning', o1, 250, 0.5),
               (lab, 'morning-afternoon', 'afternoon', o2, 250, 0.5)]
rows = []
md += ['', '## 3. Every series against the exact binomial', '',
       '| source | design | score | experiments | n | p | variance / expected | chi-square p | |CR|>=3: obs/exact | obs/normal |',
       '|---|---|---|---|---|---|---|---|---|---|']
for src, design, score, xs, n, p in series:
    X2, df, pv = gof(xs, n, p)
    r3 = tail_counts(xs, n, p, 3.0)
    vr = xs.var() / (n * p * (1 - p))
    rows.append([f'"{src}"', design, f'"{score}"', len(xs), n, p, round(vr, 4), round(X2, 1), df, round(pv, 3),
                 round(r3['obs'] / r3['exact'], 3), round(r3['obs'] / r3['normal'], 3)])
    md.append(f'| {src} | {design} | {score} | {len(xs):,} | {n} | {p} | {vr:.4f} | {pv:.3f} | {r3["obs"] / r3["exact"]:.3f} | {r3["obs"] / r3["normal"]:.3f} |')
write_csv('all_series_fit.csv', ['source', 'design', 'score', 'experiments', 'n', 'p', 'variance_ratio', 'chi_square', 'df', 'p_value',
                                 'cr3_obs_over_exact', 'cr3_obs_over_normal'], rows)

# =============================================================================
# 4. Digit matching (p = 1/10): normal versus exact
# =============================================================================
md += ['', '## 4. Digit matching (p = 1/10): two-sided counts against normal and exact', '']
rows = []
for lab, xs, n in (('Entropy Key, straight, 500 pairs', ek['straight'].astype(np.int64), 500),
                   ('random.org, straight, 268 pairs', rd['straight'].astype(np.int64), 268)):
    md += [f'**{lab}** ({len(xs):,} experiments)', '', '| t | observed | normal | exact | obs/normal | obs/exact |', '|---|---|---|---|---|---|']
    for t in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5):
        r = tail_counts(xs, n, 0.1, t)
        rows.append([f'"{lab}"', t, r['obs'], round(r['normal'], 1), round(r['exact'], 1), round(r['obs'] / r['normal'], 3), round(r['obs'] / r['exact'], 3)])
        md.append(f"| {t} | {r['obs']:,} | {r['normal']:,.0f} | {r['exact']:,.0f} | {r['obs'] / r['normal']:.2f} | {r['obs'] / r['exact']:.3f} |")
    md.append('')
    # upper tail at 3.5 SD: continuity correction
    mu_, sd_ = n * 0.1, np.sqrt(n * 0.09); hi = int(np.ceil(mu_ + 3.5 * sd_ - 1e-9))
    up_obs = int((xs >= hi).sum()); up_ex = len(xs) * binom.sf(hi - 1, n, 0.1)
    up_cc = len(xs) * norm.sf((hi - 0.5 - mu_) / sd_); up_n = len(xs) * norm.sf(3.5)
    md += [f'Upper tail at 3.5 SD: observed {up_obs:,}; exact {up_ex:,.0f}; normal {up_n:,.0f}; continuity-corrected {up_cc:,.0f}.', '']
write_csv('digit_matching.csv', ['series', 't', 'observed', 'normal', 'exact', 'obs_over_normal', 'obs_over_exact'], rows)

# =============================================================================
# 5. True size of nominal one-sided z tests in small forced-choice designs
# =============================================================================
md += ['## 5. True size of nominal one-sided tests', '', '| design | n | p | z>=1.645 (5%) | z>=1.96 (2.5%) | z>=2.326 (1%) |', '|---|---|---|---|---|---|']
rows = []
for lab, n, p in (('One Zener run', 25, 0.2), ('Four Zener runs', 100, 0.2), ('Eight Zener runs', 200, 0.2),
                  ('Four-choice, 40 trials', 40, 0.25), ('Four-choice, 100 trials', 100, 0.25),
                  ('Binary, 100 trials', 100, 0.5), ('Binary, 500 trials', 500, 0.5), ('Binary, 1,000 trials', 1000, 0.5)):
    mu_, sd_ = n * p, np.sqrt(n * p * (1 - p)); cells = []
    for z in (1.645, 1.96, 2.326):
        k = int(np.ceil(mu_ + z * sd_ - 1e-9)); cells.append((k, 100 * binom.sf(k - 1, n, p)))
    rows.append([f'"{lab}"', n, p] + [f'{c[1]:.2f}' for c in cells] + [c[0] for c in cells])
    md.append(f'| {lab} | {n} | {p} | ' + ' | '.join(f'{c[1]:.2f}% (X>={c[0]})' for c in cells) + ' |')
write_csv('small_designs.csv', ['design', 'n', 'p', 'size_5pct', 'size_2.5pct', 'size_1pct', 'k_5', 'k_2.5', 'k_1'], rows)

with open(os.path.join(RES, 'summary.md'), 'w') as fh:
    fh.write('\n'.join(md) + '\n')
print('\n'.join(md))
