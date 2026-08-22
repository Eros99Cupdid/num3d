from num3d import Num3D
import math
import numpy as np  # untuk sampling
import random

class SolusiSet:
    """
    Mewakili himpunan solusi dari operasi non-injektif.
    Bisa berupa:
      - elemen diskrit (list Num3D)
      - generator dengan parameter (fungsi)
    """
    def __init__(self, elements=None, generator=None, params=None, rank=0, default_params=None):
        self.elements = elements if elements is not None else []
        self.generator = generator
        self.params = params if params is not None else {}
        self.rank = rank
        self.default_params = default_params if default_params is not None else {}

    def canonical(self):
        """Mengembalikan satu solusi representatif (default)."""
        if self.elements:
            return self.elements[0]
        if self.generator:
            kwargs = {}
            for key in self.params:
                if key in self.default_params:
                    kwargs[key] = self.default_params[key]
                else:# default umum
                    if key == 'k':      # parameter integer (cabang)
                        kwargs[key] = 0
                    elif key == 'u':    # parameter sudut
                        kwargs[key] = 0.0
                    elif key == 't':    # parameter garis
                        kwargs[key] = 0.0
                    else:
                        kwargs[key] = 0.0
            return self.generator(**kwargs)
        raise ValueError("SolusiSet kosong.")

    def set_default(self, **kwargs):
        """Atur parameter default (kanonik)."""
        for k in kwargs:
            if k not in self.params:
                raise ValueError(f"Parameter '{k}' tidak dikenal.")
        self.default_params.update(kwargs)

    def select(self, **kwargs):
        """Pilih satu solusi dengan parameter tertentu."""
        if self.elements:
            if kwargs:
                raise ValueError("Solusi diskrit tidak memiliki parameter.")
            return self.elements[0]
        if not self.generator:
            raise ValueError("Tidak ada generator.")
        for k in kwargs:
            if k not in self.params:
                raise ValueError(f"Parameter '{k}' tidak dikenal. Tersedia: {list(self.params.keys())}")
        return self.generator(**kwargs)

    def sample(self, n=10):
        """
        Ambil n sampel solusi. Untuk rank>0, hasilnya list Num3D.
        """
        if self.elements:
            return self.elements[:n]
        if not self.generator:
            return []

        if 'k' in self.params:
            k_values = range(0, 5)  # misal beberapa cabang integer
            n_u = max(1, n // len(k_values))
            result = []
            for k in k_values:
                for u in np.linspace(0, 2*math.pi, n_u):
                    result.append(self.select(k=k, u=u))
            return result

        if self.rank == 1:
            # Satu parameter, ambil linspace dalam domain
            key = list(self.params.keys())[0]
            if 'u' in self.params and self.params['u'] == '[0, 2π)':
                vals = np.linspace(0, 2*math.pi, n)
            else:
                vals = np.linspace(-5, 5, n)
            return [self.select(**{key: v}) for v in vals]
        
        elif self.rank == 2:
            # Dua parameter, ambil grid sederhana
            keys = list(self.params.keys())[:2]
            n1 = int(math.sqrt(n))
            n2 = n // n1
            result = []
            for v1 in np.linspace(0, 1, n1):
                for v2 in np.linspace(0, 1, n2):
                    kwargs = {keys[0]: v1, keys[1]: v2}
                    result.append(self.select(**kwargs))
            return result

        else:
            # rank >= 3 : ambil acak
            keys = list(self.params.keys())
            result = []
            for _ in range(n):
                kwargs = {}
                for key in keys:
                    if key == 'k':
                        kwargs[key] = random.randint(0, 5)
                    elif key == 'u':
                        kwargs[key] = random.uniform(0, 2*math.pi)
                    else:
                        kwargs[key] = random.uniform(-5, 5)
                result.append(self.select(**kwargs))
            return result

    # ---------- OPERATOR ARITMATIKA (dengan Num3D) ----------
    def __add__(self, other):
        return self.canonical() + other
    def __radd__(self, other):
        return other + self.canonical()
    def __sub__(self, other):
        return self.canonical() - other
    def __rsub__(self, other):
        return other - self.canonical()
    def __mul__(self, other):
        return self.canonical() * other
    def __rmul__(self, other):
        return other * self.canonical()
    def __truediv__(self, other):
        return self.canonical() / other
    def __rtruediv__(self, other):
        return other / self.canonical()
    def __pow__(self, power):
        return self.canonical() ** power
    def __rpow__(self, base):
        return base ** self.canonical()

    def __repr__(self):
        if self.elements:
            return f"SolusiSet(rank=0, elements={self.elements})"
        return f"SolusiSet(rank={self.rank}, params={self.params})"

class NumL(SolusiSet):
    """Garis solusi: rank=1."""
    def __init__(self, basis, arah):
        def generator(t):
            return basis + arah * t
        super().__init__(generator=generator, params={'t': 'R'}, rank=1)
        self.basis = basis
        self.arah = arah

class NumA(SolusiSet):
    """Bidang/luas solusi: rank=2."""
    def __init__(self, basis, v1, v2):
        def generator(u, v):
            return basis + v1 * u + v2 * v
        super().__init__(generator=generator, params={'u': 'R', 'v': 'R'}, rank=2)
        self.basis = basis
        self.v1 = v1
        self.v2 = v2

class NumV(SolusiSet):
    """Volume solusi: rank=3."""
    def __init__(self, basis, v1, v2, v3):
        def generator(u, v, w):
            return basis + v1 * u + v2 * v + v3 * w
        super().__init__(generator=generator, params={'u': 'R', 'v': 'R', 'w': 'R'}, rank=3)
        self.basis = basis
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3