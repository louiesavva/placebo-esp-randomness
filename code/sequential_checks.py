"""
sequential_checks.py - tests for dependence between successive bits in the main series,
using random.org's public archive files (see rebuild_from_archive.py).

    python code/sequential_checks.py /folder/with/random.org/text/files

Uses exactly the bits the experiments used: characters 108 to 8,388,607 of each daily file,
19,008,341,000 bits for all 2,266 files. Each file is tested on its own and the results are
summed, so no comparison crosses from one day's file to the next. Under the null hypothesis
of independent fair bits:

1. Proportion of ones. The count of ones is Binomial(N, 1/2).
2. Runs. The number of changes between successive bits is Binomial(N - files, 1/2), because
   the XOR of two independent fair bits is a fair bit independent of the rest.
3. Agreement at lags 2 to 16, 250 and 500 (lag 1 is the runs test). For lag k the number of pairs (i, i + k) that
   agree is Binomial(N - k x files, 1/2), by the same argument. Lags 250 and 500 match the
   design: bit i and bit i + 250 are trial j of the morning and of the afternoon in the same
   experiment, and lag 500 is the same trial in consecutive experiments.
4. Run lengths. Complete runs (the first and last run of each file left out) have
   P(length = k) = 2^-k; lengths 1 to 24 and 25 or more are compared by chi-square.
5. Day to day. A z score for the proportion of ones and for the runs is computed for each
   file. If days differed only by chance, the mean of the squared z scores would be close
   to 1 (a chi-square test with one degree of freedom per file).

Writes results/sequential_checks.md when all 2,266 files are present.
"""
import os, sys, tarfile, zipfile
import numpy as np
from scipy.stats import chi2, norm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
SKIP, USED = 107, 8_388_500
LAGS = list(range(2, 17)) + [250, 500]          # lag 1 is the runs test
MAXLEN = 25

files = set(open(os.path.join(DATA, 'random_org_files.txt')).read().split())


def texts(folder):
    for root, _, names in os.walk(folder):
        for name in sorted(names):
            path = os.path.join(root, name)
            if name.endswith('.txt'):
                yield name, open(path, 'rb').read()
            elif name.endswith('.zip'):
                with zipfile.ZipFile(path) as z:
                    for member in z.namelist():
                        if member.endswith('.txt'):
                            yield os.path.basename(member), z.read(member)
            elif name.endswith(('.tgz', '.tar.gz')):
                with tarfile.open(path) as t:
                    for member in t.getmembers():
                        if member.isfile() and member.name.endswith('.txt'):
                            yield os.path.basename(member.name), t.extractfile(member).read()


folder = os.path.abspath(os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '.'))
if not os.path.isdir(folder):
    sys.exit(f'No such folder: {folder}')
print(f'Searching {folder}')

seen = set()
ones = changes = 0
agree = np.zeros(len(LAGS), dtype=np.int64)
runlen = np.zeros(MAXLEN + 1, dtype=np.int64)          # index k = length k; MAXLEN = MAXLEN or more
z_ones, z_runs = [], []
for name, raw in texts(folder):
    if name not in files or name in seen:
        continue
    b = np.frombuffer(raw, dtype=np.uint8)
    if len(b) != SKIP + USED + 1 or b.min() < 48 or b.max() > 49:
        print(f'{name}: not a complete random.org text file; skipped')
        continue
    seen.add(name)
    x = (b[SKIP:SKIP + USED] - 48).astype(np.int8)
    o = int(x.sum()); ones += o
    sw = x[1:] != x[:-1]
    c = int(sw.sum()); changes += c
    for j, k in enumerate(LAGS):
        agree[j] += int((x[k:] == x[:-k]).sum())
    edges = np.flatnonzero(sw)                         # a run ends at each change
    lengths = np.diff(edges)                           # complete runs only
    runlen += np.bincount(np.minimum(lengths, MAXLEN), minlength=MAXLEN + 1)[:MAXLEN + 1]
    z_ones.append((o - USED / 2) / np.sqrt(USED / 4))
    z_runs.append((c - (USED - 1) / 2) / np.sqrt((USED - 1) / 4))
    if len(seen) % 100 == 0:
        print(f'{len(seen):,} files read')

F = len(seen)
if not F:
    sys.exit('No random.org text files found.')
N = F * USED


def zline(label, count, n):
    z = (count - n / 2) / np.sqrt(n / 4)
    return f'| {label} | {n:,} | {count / n:.7f} | {z:+.2f} | {2 * norm.sf(abs(z)):.3f} |'


md = ['# Sequential checks on the main series', '',
      f'{F:,} of 2,266 daily random.org files; {N:,} bits, exactly those the experiments used. '
      'Each file is tested separately and the counts summed.', '',
      '## Proportions under the null hypothesis of independent fair bits (expected 0.5)', '',
      '| statistic | pairs or bits | proportion | z | two-sided p |', '|---|---|---|---|---|',
      zline('ones', ones, N),
      zline('changes between successive bits (runs)', changes, N - F)]
for j, k in enumerate(LAGS):
    md.append(zline(f'agreement at lag {k}', int(agree[j]), N - k * F))
ntest = 2 + len(LAGS)
md += ['', f'{ntest} tests; Bonferroni threshold for p = .05 is {0.05 / ntest:.4f}. '
       'Lag 250 compares the same trial in the morning and afternoon halves of an experiment; '
       'lag 500 the same trial in consecutive experiments.', '']

total = runlen[1:].sum()
p = np.array([0.5 ** k for k in range(1, MAXLEN)] + [0.5 ** (MAXLEN - 1)])
exp = total * p
obs = runlen[1:]
X2 = float(((obs - exp) ** 2 / exp).sum())
md += ['## Run lengths', '',
       f'{total:,} complete runs. Chi-square against P(length = k) = 2^-k, lengths 1 to {MAXLEN - 1} and '
       f'{MAXLEN} or more: {X2:.1f} on {MAXLEN - 1} df, p = {chi2.sf(X2, MAXLEN - 1):.3f}.', '',
       '| length | observed | expected | observed / expected |', '|---|---|---|---|']
md += [f'| {k if k < MAXLEN else f"{MAXLEN}+"} | {obs[k - 1]:,} | {exp[k - 1]:,.0f} | {obs[k - 1] / exp[k - 1]:.4f} |'
       for k in (1, 2, 3, 4, 5, 10, 15, 20, MAXLEN)]
md += ['']

zo, zr = np.array(z_ones), np.array(z_runs)
md += ['## Day to day', '',
       '| per-file z score | mean | mean of z squared (expected 1) | chi-square p |', '|---|---|---|---|',
       f'| ones | {zo.mean():+.3f} | {(zo ** 2).mean():.3f} | {chi2.sf((zo ** 2).sum(), F):.3f} |',
       f'| runs | {zr.mean():+.3f} | {(zr ** 2).mean():.3f} | {chi2.sf((zr ** 2).sum(), F):.3f} |', '']

print('\n'.join(md))
if F == 2266:
    with open(os.path.join(HERE, '..', 'results', 'sequential_checks.md'), 'w') as fh:
        fh.write('\n'.join(md) + '\n')
    print('Written to results/sequential_checks.md')
