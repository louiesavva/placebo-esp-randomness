"""
check_rdata.py - confirm that data/placebo_esp.RData (built by make_rdata.R) holds exactly
the same values as the .npz files in data/.

    python code/check_rdata.py

Reads the R file with this project's own reader (no R needed) and compares every column,
file list and bit vector with its .npz counterpart.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rdata_reader import load

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
npz = lambda f: np.load(os.path.join(DATA, f))
lines = lambda f: open(os.path.join(DATA, f)).read().split()

r = load(os.path.join(DATA, 'placebo_esp.RData'))


def cols(name):
    """Columns of an R data frame, by name."""
    _, values, attrs = r[name]
    names = next(v for t, v, _ in attrs if t == 'names')
    return dict(zip(names, values))


def unpack(f):
    z = npz(f)
    return np.unpackbits(z['bits'])[:int(z['n_bits'])]


expected = {
    'main': {'morning': npz('random_org_main_morning.npz')['morning'],
             'afternoon': npz('random_org_main_afternoon.npz')['afternoon']},
    'short': dict(npz('random_org_125.npz')),
    'prng': dict(npz('prng.npz')),
    'entropy_key': dict(npz('entropy_key_digits.npz')),
    'digits': dict(npz('random_org_digits.npz')),
    'other_mcr': {'morning': npz('other_runs.npz')['mcr_morning'], 'afternoon': npz('other_runs.npz')['mcr_afternoon']},
    'other_aao1rng': {'morning': npz('other_runs.npz')['aao1rng_morning'],
                      'afternoon': npz('other_runs.npz')['aao1rng_afternoon']},
    'zener': {'hits': npz('random_org_zener.npz')['hits']},
}
ok = True
for name, columns in expected.items():
    got = cols(name)
    for col, want in columns.items():
        same = col in got and np.array_equal(np.asarray(got[col], dtype=np.int64), want.astype(np.int64))
        ok &= same
        print(f'{name}${col}: {len(want):,} values, identical: {same}')
for name, f in (('main_files', 'random_org_files.txt'), ('zener_files', 'random_org_zener_files.txt')):
    same = list(r[name]) == lines(f)
    ok &= same
    print(f'{name}: {len(r[name]):,} names, identical: {same}')
for name, f in (('bits_2012_06_17', 'random_org_2012-06-17_bits.npz'), ('bits_2012_08_04', 'random_org_2012-08-04_bits.npz')):
    same = np.array_equal(np.asarray(r[name], dtype=np.int64), unpack(f).astype(np.int64))
    ok &= same
    print(f'{name}: {len(r[name]):,} bits, identical: {same}')
print('ALL IDENTICAL' if ok else 'MISMATCH')
sys.exit(0 if ok else 1)
