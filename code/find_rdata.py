"""
find_rdata.py - list every R workspace (.RData / .rda) under a folder, say what
is in it, and flag series we have not seen yet (the hardware-key data).

Also lists R scripts that mention /dev/random, urandom, entropy keys or readRandom.

Needs: Python 3 and numpy.
Run:   python find_rdata.py "C:\\path\\to\\search"      (or several folders)
       python find_rdata.py C:\\ D:\\                    (whole drives; slow but thorough)
"""
import sys, os, gzip, bz2, lzma, struct, time, datetime
import numpy as np

KNOWN = {   # first three values of each series we already have
    (-1.13842, -1.26491, 1.01193): 'PRNG series (already have: prng)',
    (0.12649, 2.40333, -0.12649): 'mcr.RData series (already have)',
    (-1.01193, 0.63246, 2.7828): 'aao1rng series (already have)',
    (1.3914, -0.25298, -0.63246): 'random.org main series (already have)',
    (-1.26491, 0.12649, 1.01193): 'random.org main series, later batch (already have: ppp1)',
    (0.98387, 0.44721, 1.52053): 'random.org 125+125 series (already have: randmorn)',
    (-0.14907, -1.78885, 0.14907): 'random.org digit matching, 500 pairs (already have)',
    (0.44795, -1.79182, -0.36651): 'random.org digit matching, 268 pairs (already have)',
    (0.31667, -1.81726, -0.38651): 'random.org digit matching, 268 pairs, own-frequency version (already have)',
}
SERIES = ('mcr', 'scr', 'cr')   # vectors that identify a series
sys.setrecursionlimit(100000)

def opener(path):
    with open(path, 'rb') as fh:
        head = fh.read(6)
    if head[:2] == b'\x1f\x8b': return gzip.open(path, 'rb')
    if head[:3] == b'BZh':      return bz2.open(path, 'rb')
    if head[:6] == b'\xfd7zXZ\x00': return lzma.open(path, 'rb')
    return open(path, 'rb')

class Reader:
    def __init__(self, f):
        self.f = f; self.refs = []; self.objects = []
    def i32(self): return struct.unpack('>i', self.f.read(4))[0]
    def length(self):
        n = self.i32()
        if n == -1: hi, lo = self.i32(), self.i32(); n = (hi << 32) + lo
        return n
    def skip(self, nbytes, keep_first=0, dtype='>f8'):
        first = np.frombuffer(self.f.read(min(nbytes, keep_first)), dtype=dtype) if keep_first else None
        left = nbytes - (min(nbytes, keep_first) if keep_first else 0)
        while left > 0:
            chunk = self.f.read(min(left, 1 << 24))
            if not chunk: break
            left -= len(chunk)
        return first
    def strvec(self):
        self.i32(); n = self.i32(); return [self.item() for _ in range(n)]
    def item(self):
        flags = self.i32(); t = flags & 0xFF
        has_attr, has_tag = flags & (1 << 9), flags & (1 << 10)
        if t == 255: return self.refs[(flags >> 8) - 1]
        if t in (254, 253, 252, 251, 250, 242, 241): return None
        if t in (249, 248, 247): v = ('ns', self.strvec()); self.refs.append(v); return v
        if t == 238:   # ALTREP
            self.item(); state = self.item(); self.item(); return state
        if t == 1: name = self.item(); self.refs.append(name); return name
        if t == 9:
            n = self.i32(); return None if n == -1 else self.f.read(n).decode('latin1')
        if t == 4:
            env = {}; self.refs.append(env); self.i32(); self.item(); self.item(); self.item(); self.item(); return env
        if t in (2, 3, 5, 6, 17, 240, 239):
            attr = self.item() if has_attr else None
            tag = self.item() if has_tag else None
            car = self.item(); cdr = self.item()
            if t == 2:
                rest = cdr if isinstance(cdr, list) else []
                return [(tag, car, attr)] + rest
            return ('code',)
        if t in (7, 8): n = self.i32(); self.f.read(n); return ('builtin',)
        if t in (10, 13, 14, 15, 24):
            n = self.length()
            dt = {10: '>i4', 13: '>i4', 14: '>f8', 15: '>c16', 24: '>u1'}[t]
            size = np.dtype(dt).itemsize
            first = self.skip(n * size, keep_first=min(n, 4) * size, dtype=dt)
            attr = self.item() if has_attr else None
            return ('vec', dt, n, first, attr)
        if t == 16:
            n = self.length(); v = [self.item() for _ in range(n)]
            if has_attr: self.item()
            return v
        if t in (19, 20):
            n = self.length(); v = [self.item() for _ in range(n)]
            attr = self.item() if has_attr else None
            return ('list', v, attr)
        if t == 25:
            if has_attr: self.item()
            return ('s4',)
        raise ValueError(f'unhandled R type {t}')

def summarise(path):
    f = opener(path)
    magic = f.read(5)
    if magic not in (b'RDX2\n', b'RDX3\n'):
        return None
    f.read(2)
    r = Reader(f)
    version = r.i32(); r.i32(); r.i32()
    if version == 3:
        n = r.i32(); f.read(n)
    top = r.item()
    out = []
    files = None
    for tag, val, attr in top:
        if isinstance(val, tuple) and val[0] == 'vec':
            _, dt, n, first, a = val
            if a and isinstance(a, list):
                for at, av, _ in a:
                    if at == 'levels' and isinstance(av, list) and av and isinstance(av[0], str):
                        files = (len(av), av[0], av[-1])
            if n >= 100000:
                out.append((tag, n, first))
    return out, files

def main(roots):
    hits = []
    for root in roots:
        for dirpath, _, names in os.walk(root):
            for name in names:
                p = os.path.join(dirpath, name)
                low = name.lower()
                try:
                    mtime = datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime('%Y-%m-%d')
                    size = os.path.getsize(p)
                except OSError:
                    continue
                if low.endswith(('.rdata', '.rda')):
                    try:
                        res = summarise(p)
                    except Exception as e:
                        print(f'\n{p}\n   could not read: {e}'); continue
                    if res is None: continue
                    vecs, files = res
                    print(f'\n{p}   ({size/1e6:.1f} MB, saved {mtime})')
                    if files: print(f'   file list: {files[0]} files, {files[1]} ... {files[2]}')
                    new = False
                    seen = False
                    for tag, n, first in vecs:
                        label = ''
                        if tag in SERIES and first is not None and len(first) >= 3 and not seen:
                            key = tuple(round(float(v), 5) for v in first[:3])
                            label = KNOWN.get(key, '*** NOT SEEN BEFORE ***')
                            new = label.startswith('***'); seen = True
                        print(f'   {tag}: {n:,} values {label}')
                    if not seen and vecs and not files:
                        new = True
                        print('   (large vectors, no random.org file list: worth a look)')
                    if new: hits.append(p)
                elif low.endswith(('.r', '.txt', '.rhistory')) and size < 5_000_000:
                    try:
                        text = open(p, 'r', errors='ignore').read()
                    except OSError:
                        continue
                    keys = [k for k in ('/dev/random', 'urandom', 'entropy', 'ekey', 'readRandom') if k in text]
                    if keys and not name[:4].isdigit():
                        print(f'\nScript {p} ({mtime}) mentions: {", ".join(keys)}')
    print('\n==== Candidates for the hardware-key series ====')
    for h in hits: print('  ', h)
    if not hits: print('   none found in these folders')

if __name__ == '__main__':
    main(sys.argv[1:] or ['.'])
