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

The folder is searched recursively. It may hold the daily .txt files themselves, or .zip,
.tgz or .tar.gz archives containing them (such as random.org's monthly archives). Files
whose dates are not in data/random_org_files.txt are ignored. With all 2,266 files present,
the script checks every one of the 38,016,682 experiments.
"""
import os, sys, tarfile, zipfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
SKIP, USED, PER_FILE = 107, 8_388_500, 16_777

files = open(os.path.join(DATA, 'random_org_files.txt')).read().split()
index = {f: i for i, f in enumerate(files)}
m = np.load(os.path.join(DATA, 'random_org_main_morning.npz'))['morning']
a = np.load(os.path.join(DATA, 'random_org_main_afternoon.npz'))['afternoon']


def texts(folder):
    """Yield (file name, contents) for every random.org text file under folder."""
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


seen, mismatched = set(), []
for name, raw in texts(sys.argv[1] if len(sys.argv) > 1 else '.'):
    if name not in index or name in seen:
        continue
    bits = np.frombuffer(raw, dtype=np.uint8) - ord('0')
    if len(bits) != SKIP + USED + 1 or bits.max() > 1:
        print(f'{name}: not a random.org text file of 8,388,608 bits; skipped')
        continue
    seen.add(name)
    counts = bits[SKIP:SKIP + USED].reshape(PER_FILE, 2, 250).sum(axis=2)
    k = index[name] * PER_FILE
    same_m = counts[:, 0] == m[k:k + PER_FILE]
    same_a = counts[:, 1] == a[k:k + PER_FILE]
    bad = int((~(same_m & same_a)).sum())
    if bad:
        mismatched.append((name, bad))
        print(f'{name}: {bad:,} of {PER_FILE:,} experiments differ from data/')
    elif len(seen) % 100 == 0:
        print(f'{len(seen):,} files checked, all identical so far')

ok = len(seen) - len(mismatched)
print(f'\nFiles found: {len(seen):,} of {len(files):,}. Rebuilt identically: {ok:,}. '
      f'Experiments checked: {len(seen) * PER_FILE:,} of {len(files) * PER_FILE:,}.')
missing = [f for f in files if f not in seen]
if missing and len(missing) <= 20:
    print('Not found: ' + ', '.join(missing))
elif missing:
    print(f'Not found: {len(missing):,} files, from {missing[0]} to {missing[-1]}')
for name, bad in mismatched:
    print(f'DIFFERS: {name} ({bad:,} experiments)')
sys.exit(0 if seen and not mismatched else 1)
