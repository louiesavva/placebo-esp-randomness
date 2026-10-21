"""
rebuild_from_archive.py - rebuild main-series experiments from random.org's public archive
and check them against data/.

    python code/rebuild_from_archive.py /folder/with/random.org/text/files

random.org distributes its daily pregenerated files at https://archive.random.org/text
(free via BitTorrent, one torrent per month). Each text file YYYY-MM-DD.txt holds
8,388,608 random bits written as the characters 0 and 1. The main series used, from each
file, characters 108 to 8,388,607 (counting from 1) in order: 8,388,500 bits, taken 500 at a
time, the first 250 of each block as the morning hits and the next 250 as the afternoon
hits. The first 107 characters and the last one were not used. This rule was established
by comparing the archive file for 17 June 2012 with the bits of that day saved in 2012.

The folder may hold the .txt files themselves or .zip files containing them. Files whose
dates are not in data/random_org_files.txt are ignored.
"""
import os, sys, zipfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
SKIP, USED, PER_FILE = 107, 8_388_500, 16_777

files = open(os.path.join(DATA, 'random_org_files.txt')).read().split()
index = {f: i for i, f in enumerate(files)}
m = np.load(os.path.join(DATA, 'random_org_main_morning.npz'))['morning']
a = np.load(os.path.join(DATA, 'random_org_main_afternoon.npz'))['afternoon']


def texts(folder):
    for name in sorted(os.listdir(folder)):
        path = os.path.join(folder, name)
        if name.endswith('.txt'):
            yield name, open(path, 'rb').read()
        elif name.endswith('.zip'):
            with zipfile.ZipFile(path) as z:
                for member in z.namelist():
                    if member.endswith('.txt'):
                        yield os.path.basename(member), z.read(member)


checked = matched = 0
for name, raw in texts(sys.argv[1] if len(sys.argv) > 1 else '.'):
    if name not in index:
        continue
    bits = np.frombuffer(raw, dtype=np.uint8) - ord('0')
    if len(bits) != SKIP + USED + 1 or bits.max() > 1:
        print(f'{name}: not a random.org text file of 8,388,608 bits; skipped')
        continue
    counts = bits[SKIP:SKIP + USED].reshape(PER_FILE, 2, 250).sum(axis=2)
    k = index[name] * PER_FILE
    same = np.array_equal(counts[:, 0], m[k:k + PER_FILE]) and np.array_equal(counts[:, 1], a[k:k + PER_FILE])
    checked += 1; matched += same
    print(f'{name}: experiments {k:,} to {k + PER_FILE - 1:,} (from 0) rebuilt; identical to data/: {same}')
print(f'{matched} of {checked} files rebuilt identically.')
sys.exit(0 if checked and matched == checked else 1)
