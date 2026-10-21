"""
extract.py - turn the original R workspaces (2012-2013) into compact data files.

Usage:
    python code/extract.py /path/to/folder/with/RData/files

Expects the original file names (as stored in Google Drive, folder "RData"):
    Spencp.RData, rand1.RData        random.org main morning-afternoon series (Example 1B)
    prng.RData                       R Mersenne Twister morning-afternoon series (Example 1A)
    mcrx.RData                       Entropy Key digit matching (Example 2)
    randorgorigA.RData               random.org digit matching, 268 pairs (Example 3)
    randmorn.RData                   random.org morning-afternoon, 125 + 125 trials
    mcr.RData, aao1rng.RData         two further morning-afternoon runs (source not recorded)
    rhine1.RData                     random.org Zener-card series, 625 trials, p = 1/5 (2013)
    x.RData, rled.RData, pppv1.RData October 2012 workspaces used as provenance records

Writes to data/ (counts of hits per experiment, unsigned 8-bit, compressed):
    random_org_main_morning.npz, random_org_main_afternoon.npz
                             hits out of 250 in each half, 38,016,682 experiments
    random_org_files.txt     the 2,266 daily random.org files, in order (16,777 experiments each)
    prng.npz                 morning, afternoon  (hits out of 250 each)
    entropy_key_digits.npz   straight (of 500), plus1, minus1 (of 490): digit matches, p = 1/10
    random_org_digits.npz    straight, plus1..plus7, minus1..minus7 (of 268): digit matches, p = 1/10
    random_org_125.npz       morning, afternoon (hits out of 125 each)
    other_runs.npz           mcr_morning, mcr_afternoon, aao1rng_morning, aao1rng_afternoon (of 250)
    random_org_zener.npz     hits (of 625), Zener symbols, p = 1/5, 473,594 experiments
    random_org_zener_files.txt  the 2,266 daily random.org files in the order processed (209 experiments each)
    random_org_2012-06-17_bits.npz, random_org_2012-08-04_bits.npz
                             the 8,388,500 bits used from each of two daily files, as divided into experiments
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from rdata_reader import load

SRC = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = os.path.join(os.path.dirname(__file__), '..', 'data')
os.makedirs(OUT, exist_ok=True)


def hits(cr, n, p):
    """Convert stored critical ratios back to whole-number hit counts."""
    sd = np.sqrt(n * p * (1 - p))
    x = cr * sd + n * p
    xr = np.round(x)
    bad = np.abs(x - xr) > 1e-6
    return xr.astype(np.int64), bad


def path(name):
    return os.path.join(SRC, name)


def report(label, arr):
    print(f'  {label}: {len(arr):,} experiments, mean {arr.mean():.4f}, variance {arr.var():.4f}')


# ---- random.org main series (Example 1B) --------------------------------------
print('random.org main series')
r1 = load(path('rand1.RData'), wanted={'mcr', 'acr'})
m, bad_m = hits(r1['mcr'], 250, 0.5)
a, bad_a = hits(r1['acr'], 250, 0.5)
assert not bad_m.any() and not bad_a.any(), 'rand1 copy should be undamaged'
sp = load(path('Spencp.RData'), wanted={'mdevv', 'adevv', 'tcr', 'filenames'})
for label, stored, clean in (('morning', sp['mdevv'], m - 125), ('afternoon', sp['adevv'], a - 125)):
    damaged = np.abs(stored - np.round(stored)) > 1e-6
    agree = np.all(np.round(stored[~damaged]).astype(np.int64) == clean[~damaged])
    print(f'  Spencp {label}: {damaged.sum()} damaged values; undamaged values agree with rand1: {agree}')
tcr = sp['tcr']; k = tcr * np.sqrt(125); damaged = np.abs(k - np.round(k)) > 1e-6
agree = np.all(np.round(k[~damaged]).astype(np.int64) == (m - a)[~damaged])
print(f'  Spencp tcr (morning minus afternoon): {damaged.sum()} damaged values; undamaged agree: {agree}')
np.savez_compressed(os.path.join(OUT, 'random_org_main_morning.npz'), morning=m.astype(np.uint8))
np.savez_compressed(os.path.join(OUT, 'random_org_main_afternoon.npz'), afternoon=a.astype(np.uint8))
files = sp['filenames'][2]
with open(os.path.join(OUT, 'random_org_files.txt'), 'w') as fh:
    fh.write('\n'.join(files) + '\n')
report('morning', m); report('afternoon', a)

# ---- PRNG (Example 1A) ---------------------------------------------------------
print('PRNG series')
pr = load(path('prng.RData'), wanted={'mcr', 'acr'})
pm, b1 = hits(pr['mcr'], 250, 0.5); pa, b2 = hits(pr['acr'], 250, 0.5)
assert not b1.any() and not b2.any()
np.savez_compressed(os.path.join(OUT, 'prng.npz'), morning=pm.astype(np.uint8), afternoon=pa.astype(np.uint8))
report('morning', pm)

# ---- Entropy Key digit matching (Example 2) -------------------------------------
print('Entropy Key digit matching')
ek = load(path('mcrx.RData'), wanted={'scr', 'Acr', 'Bcr'})
s, b0 = hits(ek['scr'], 500, 0.1); p1, b1 = hits(ek['Acr'], 490, 0.1); m1, b2 = hits(ek['Bcr'], 490, 0.1)
assert not (b0.any() or b1.any() or b2.any())
np.savez_compressed(os.path.join(OUT, 'entropy_key_digits.npz'), straight=s.astype(np.uint8),
                    plus1=p1.astype(np.uint8), minus1=m1.astype(np.uint8))
report('straight', s)

# ---- random.org digit matching (Example 3) --------------------------------------
print('random.org digit matching')
names = ['cr'] + [f'crp{i}' for i in range(1, 8)] + [f'crm{i}' for i in range(1, 8)]
ro = load(path('randorgorigA.RData'), wanted=set(names))
out = {}
for nm in names:
    x, b = hits(ro[nm], 268, 0.1); assert not b.any()
    key = 'straight' if nm == 'cr' else ('plus' if nm[2] == 'p' else 'minus') + nm[3:]
    out[key] = x.astype(np.uint8)
np.savez_compressed(os.path.join(OUT, 'random_org_digits.npz'), **out)
report('straight', out['straight'])

# ---- smaller morning-afternoon runs ---------------------------------------------
print('random.org 125 + 125 series')
rm = load(path('randmorn.RData'), wanted={'mcr', 'acr'})
x1, b1 = hits(rm['mcr'], 125, 0.5); x2, b2 = hits(rm['acr'], 125, 0.5)
assert not b1.any() and not b2.any()
np.savez_compressed(os.path.join(OUT, 'random_org_125.npz'), morning=x1.astype(np.uint8), afternoon=x2.astype(np.uint8))
report('morning', x1)

print('other runs')
o = {}
for fname, key in (('mcr.RData', 'mcr'), ('aao1rng.RData', 'aao1rng')):
    d = load(path(fname), wanted={'mcr', 'acr'})
    x1, b1 = hits(d['mcr'], 250, 0.5); x2, b2 = hits(d['acr'], 250, 0.5)
    assert not b1.any() and not b2.any()
    o[f'{key}_morning'] = x1.astype(np.uint8); o[f'{key}_afternoon'] = x2.astype(np.uint8)
    report(key, x1)
np.savez_compressed(os.path.join(OUT, 'other_runs.npz'), **o)

# ---- random.org Zener-card series (2013) -------------------------------------------
print('random.org Zener-card series')
zn = load(path('rhine1.RData'), wanted={'dev', 'filenames'})
z = zn['dev']; assert np.all(z == np.round(z)) and len(z) == 2266 * 209
z = z.astype(np.int64)
np.savez_compressed(os.path.join(OUT, 'random_org_zener.npz'), hits=z.astype(np.uint8))
_, codes, levels = zn['filenames']
with open(os.path.join(OUT, 'random_org_zener_files.txt'), 'w') as fh:     # processing order, dates only
    fh.write('\n'.join(os.path.basename(levels[c - 1]) for c in codes) + '\n')
report('hits of 625', z)

# ---- provenance records (2012) -----------------------------------------------------
# x.RData and rled.RData each hold the raw bits of one daily random.org file (17 June and
# 4 August 2012) as they were divided into experiments; pppv1.RData holds frequency
# tables of the whole main series.
print('provenance records')
for fname, day in (('x.RData', '2012-06-17'), ('rled.RData', '2012-08-04')):
    xr = load(path(fname), wanted={'d[0]', 'd[1]', 'd[2]'})
    bits, exp_no, half = (xr['d'][1][k] for k in range(3))
    assert np.array_equal(exp_no, np.repeat(np.arange(1, 16778), 500))
    assert np.array_equal(half, np.tile(np.repeat([1, 2], 250), 16777))
    np.savez_compressed(os.path.join(OUT, f'random_org_{day}_bits.npz'),
                        bits=np.packbits(bits.astype(np.uint8)), n_bits=len(bits))
    k = files.index(f'{day}.txt') * 16777
    blocks = bits.reshape(16777, 2, 250).sum(axis=2)
    print(f'  {day}: {len(bits):,} bits, 500 consecutive bits per experiment (250 morning, then 250 afternoon); '
          f'counts equal experiments {k:,} to {k + 16776:,} of the main series (from 0): '
          f'{np.array_equal(blocks[:, 0], m[k:k + 16777]) and np.array_equal(blocks[:, 1], a[k:k + 16777])}')
pp = load(path('pppv1.RData'), wanted={'a', 'm', 't', 'n'})
same = int(pp['n'][0]) == len(m)
for key, score in (('m', m - 125), ('a', a - 125), ('t', m - a)):
    same &= np.array_equal(np.unique(score, return_counts=True)[1], pp[key])
print(f'  October 2012 frequency tables of morning, afternoon and morning minus afternoon match the data: {same}')
print('done')
