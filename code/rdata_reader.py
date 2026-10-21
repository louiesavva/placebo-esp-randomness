"""
rdata_reader.py - read numeric vectors from R workspace files (.RData) without R.

Handles R's XDR serialisation format (RDX2/RDX3), gzip/bzip2/xz/uncompressed.
Only what this project needs: numeric, integer and logical vectors, character
vectors, lists, factors (their levels), and enough of everything else (functions,
environments) to skip past it.

    from rdata_reader import load
    objs = load('Spencp.RData', wanted={'mcr', 'acr'})
    objs['mcr']   # numpy float64 array

Vectors not in `wanted` are skipped without being kept in memory.
"""
import gzip, bz2, lzma, struct, sys, threading
import numpy as np

sys.setrecursionlimit(1_000_000)


def _open(path):
    with open(path, 'rb') as fh:
        head = fh.read(6)
    if head[:2] == b'\x1f\x8b':
        return gzip.open(path, 'rb')
    if head[:3] == b'BZh':
        return bz2.open(path, 'rb')
    if head[:6] == b'\xfd7zXZ\x00':
        return lzma.open(path, 'rb')
    return open(path, 'rb')


class _Reader:
    DT = {10: '>i4', 13: '>i4', 14: '>f8', 15: '>c16', 24: '>u1'}

    def __init__(self, f, wanted):
        self.f, self.wanted, self.refs = f, wanted, []

    def i32(self):
        return struct.unpack('>i', self.f.read(4))[0]

    def length(self):
        n = self.i32()
        if n == -1:
            hi, lo = self.i32(), self.i32()
            n = (hi << 32) + lo
        return n

    def vector(self, dt, n, keep):
        size = np.dtype(dt).itemsize
        if keep:
            out = np.empty(n, dtype=dt.replace('>', '<'))
            step, i = 4_000_000, 0
            while i < n:
                m = min(step, n - i)
                out[i:i + m] = np.frombuffer(self.f.read(m * size), dtype=dt)
                i += m
            return out
        left = n * size
        while left > 0:
            chunk = self.f.read(min(left, 1 << 24))
            if not chunk:
                break
            left -= len(chunk)
        return None

    def item(self, tag=None):
        flags = self.i32()
        t = flags & 0xFF
        has_attr, has_tag = flags & (1 << 9), flags & (1 << 10)
        if t == 255:                                   # reference
            return self.refs[(flags >> 8) - 1]
        if t in (254, 253, 252, 251, 250, 242, 241):  # NULL, env markers
            return None
        if t in (249, 248, 247):                       # namespace/package/persist
            self.i32(); n = self.i32()
            v = [self.item() for _ in range(n)]
            self.refs.append(v)
            return v
        if t == 238:                                   # ALTREP
            self.item(); state = self.item(); self.item()
            return state
        if t == 1:                                     # symbol
            name = self.item(); self.refs.append(name)
            return name
        if t == 9:                                     # CHARSXP
            n = self.i32()
            return None if n == -1 else self.f.read(n).decode('latin1')
        if t == 4:                                     # environment
            env = {}; self.refs.append(env)
            self.i32(); self.item(); self.item(); self.item(); self.item()
            return env
        if t in (2, 3, 5, 6, 17, 240, 239):            # pairlists, closures, calls
            attr = self.item() if has_attr else None
            tg = self.item() if has_tag else None
            car = self.item(tag=tg if t == 2 else None)
            cdr = self.item()
            if t == 2:
                rest = cdr if isinstance(cdr, list) else []
                return [(tg, car, attr)] + rest
            return ('code',)
        if t in (7, 8):                                # builtins
            n = self.i32(); self.f.read(n)
            return ('builtin',)
        if t in self.DT:                               # atomic vectors
            n = self.length()
            keep = self.wanted is None or tag in self.wanted
            val = self.vector(self.DT[t], n, keep)
            attr = self.item() if has_attr else None
            if attr and isinstance(attr, list):
                for at, av, _ in attr:
                    if at == 'levels':
                        return ('factor', val, av)
            return val
        if t == 16:                                    # character vector
            n = self.length()
            v = [self.item() for _ in range(n)]
            if has_attr:
                self.item()
            return v
        if t in (19, 20):                              # list / expression
            n = self.length()
            v = [self.item(tag=f'{tag}[{k}]' if tag else None) for k in range(n)]
            attr = self.item() if has_attr else None
            return ('list', v, attr)
        if t == 25:                                    # S4
            if has_attr:
                self.item()
            return ('s4',)
        raise ValueError(f'unhandled R type {t}')


def load(path, wanted=None):
    """Return {name: value} for the top-level objects of an .RData file."""
    result = {}

    def run():
        f = _open(path)
        magic = f.read(5)
        if magic not in (b'RDX2\n', b'RDX3\n'):
            raise ValueError(f'{path}: not an RData file')
        f.read(2)
        r = _Reader(f, set(wanted) if wanted else None)
        version = r.i32(); r.i32(); r.i32()
        if version == 3:
            n = r.i32(); f.read(n)
        for tag, val, _ in r.item():
            result[tag] = val

    threading.stack_size(512 * 1024 * 1024)
    th = threading.Thread(target=run)
    th.start(); th.join()
    return result
