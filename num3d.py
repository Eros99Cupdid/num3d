"""
NUM3D - Aljabar 3D
==========================================
"""
from typing import Union
import math

class Num3D:
    SQRT2 = math.sqrt(2)
    SQRT3 = math.sqrt(3)
    SQRT6 = SQRT2 * SQRT3

    __slots__ = ('z', 'x', 'y')
    
    def __init__(self, z=0.0, x=0.0, y=0.0):
        self.z = float(z)
        self.x = float(x)
        self.y = float(y)

    @staticmethod
    def izin(other):
        """izin Num3D."""
        from numr import NumR
        from solusi import SolusiSet
        if isinstance(other, (int, float)):
            return Num3D(other, 0, 0)
        elif isinstance(other, Num3D):
            return other
        elif isinstance(other, NumR):
            return Num3D(other.z, other.x, other.y)
        elif isinstance(other, SolusiSet):
            return other.canonical()
        else:
            raise TypeError(f"* tidak didukung antara Num3D dan {type(other)}")

    # ---------- OPERATOR ALJABAR `*` ----------
    def __add__(self, other):
        other = Num3D.izin(other)
        return Num3D(self.z + other.z, self.x + other.x, self.y + other.y)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        other = Num3D.izin(other)
        return Num3D(self.z - other.z, self.x - other.x, self.y - other.y)

    def __rsub__(self, other):
        other = Num3D.izin(other)
        return other - self

    def __neg__(self):
        return Num3D(-self.z, -self.x, -self.y)

    def __mul__(self, other):
        other = Num3D.izin(other)
        return Num3D(
            self.z * other.z - (
            self.x * other.x + self.y * other.y),
            self.z * other.x + self.x * other.z,
            self.z * other.y + self.y * other.z
        )
    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other = Num3D.izin(other)
        return self.true_div(other)

    def __rtruediv__(self, other):
        other = Num3D.izin(other)
        return other / self

    def __pow__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                return Num3D(1, 0, 0)
            if other == 1:
                return self
            if other < 0:
                return self.inv() ** (-other)
            if isinstance(other, int):
                result = Num3D(1, 0, 0)
                for _ in range(int(other)):
                    result = result * self
                return result
            # float: gunakan exp(other * ln(self))
            other_el = Num3D(float(other), 0, 0)
            return (other_el * self.ln_principal()).exp()
        other = Num3D.izin(other)
        # pangkat dengan Num3D
        return (other * self.ln_principal()).exp()

    def __xor__(self, other):
        """Cross product quaternion (Hamilton) - mengembalikan NumR."""
        from numr import NumR
        other = Num3D.izin(other)
        self = NumR(0, self.x, self.y, self.z)
        return self ^ other

    def true_div(self, other: 'Num3D') -> 'Num3D':
        """
        Pembagian sejati: self / other = X  ->  other * X = self
        Menggunakan rumus eksplisit di basis {a, b, c}.
        """
        from solusi import SolusiSet, NumL
        other = Num3D.izin(other)
        z2 = other.z
        x2, y2 = other.x, other.y
        z1, x1, y1 = self.z, self.x, self.y
        if math.isclose(z2, 0.0, abs_tol=1e-12):
        # Pembagi nol parsial -> garis (rank=1) atau volume (rank=3) jika pembagi nol total
            if abs(x2) < 1e-12 and abs(y2) < 1e-12:
                # Pembagi nol total: semua x memenuhi jika self juga nol, selain itu tidak ada solusi
                if abs(z1) < 1e-12 and abs(x1) < 1e-12 and abs(y1) < 1e-12:
                    # Solusi seluruh ruang (rank=3)
                    def gen(u, v, w):
                        return Num3D(u, v, w)
                    return SolusiSet(generator=gen, params={'u': 'R', 'v': 'R', 'w': 'R'}, rank=3)
                else:
                    raise ZeroDivisionError("Tidak ada solusi (pembagi nol total dan bukan nol).")
            else:# Garis solusi (rank=1)
                basis = self * other.inv()
                arah = Num3D(0, y2, -x2)
                arah = arah / arah.norm()
                return NumL(basis, arah)
        # Kasus z2 != 0 -> solusi unik (rank=0)
        sq = x2*x2 + y2*y2 + z2*z2  # pasti > 0
        # Hitung komponen hasil bagi
        N_z = z1*z2 + x1*x2 + y1*y2
        N_x = x1*(sq - x2*x2) - x2*(y1*y2 + z1*z2)
        N_y = y1*(sq - y2*y2) - y2*(x1*x2 + z1*z2)
        inv_z2_sq = 1.0 / (z2*sq)
        inv_sq = inv_z2_sq*z2
        z3 = N_z*inv_sq
        x3 = N_x*inv_z2_sq
        y3 = N_y*inv_z2_sq
        return Num3D(z3, x3, y3)

    def inv(self):
        norm2 = self.z * self.z + self.x * self.x + self.y * self.y
        if norm2 == 0:
            raise ZeroDivisionError("Zero element has no inverse")
        return Num3D(self.z / norm2, -self.x / norm2, -self.y / norm2)

    # ---------- LOGARITMA BASIS ----------
    def log(self, other = None) -> 'Num3D':
        """Logaritma basis `oteher` dari `self`: log_self(other) = ln(other)/ln(self)"""
        if other is None:
            ten = Num3D(10, 0, 0)
            return self.ln_principal() / ten.ln_principal() # log_{10}(self)
        other = Num3D.izin(other)
        return  other.ln_principal() / self.ln_principal() # log_{self}(other)

    # ---------- FUNGSI ANALITIK ----------
    def exp(self):
        r = math.hypot(self.x, self.y)
        e_z = math.exp(self.z)
        if r < 1e-12:
            return Num3D(e_z, 0, 0)
        cos_r = math.cos(r)
        sin_r = math.sin(r)
        factor = e_z * sin_r / r
        return Num3D(e_z * cos_r, factor * self.x, factor * self.y)

    def ln_principal(self):
        r = math.hypot(self.x, self.y)
        norm_val = self.norm()
        if norm_val == 0:
            raise ValueError("ln(0) undefined")
        theta = math.atan2(r, self.z)
        ln_r = math.log(norm_val)
        if r == 0:
            if self.z < 0:
                # ln(-|z|) = ln(|z|) + π i
                return Num3D(math.log(-self.z), math.pi, 0)
            else:
                return Num3D(ln_r, 0, 0)
        factor = theta / r
        return Num3D(ln_r, factor * self.x, factor * self.y)

    def ln(self):
        from solusi import SolusiSet
        r = math.hypot(self.x, self.y)
        if r == 0 and self.z < 0:
            # Kasus skalar negatif: hasilnya silinder (lingkaran + cabang)
            base = math.log(-self.z)
            def generator(k, u):
                return Num3D(base,
                             (2*k + 1)*math.pi * math.cos(u),
                             (2*k + 1)*math.pi * math.sin(u))
            return SolusiSet(generator=generator, params={'k': 'Z', 'u': '[0, 2π)'}, rank=2)
        else: # Kasus umum: cabang utama + kelipatan 2πk di bidang imajiner => rank=1 (parameter k)
            principal = self.ln_principal()
            def generator(k):
                # Arah imajiner dinormalisasi
                if r > 0:
                    factor = 2 * math.pi * k / r
                    return Num3D(principal.z,
                                principal.x + factor * self.x,
                                principal.y + factor * self.y)
                else:
                    return principal  # untuk z positif real, tidak ada cabang imajiner
            return SolusiSet(generator=generator, params={'k': 'Z'}, rank=1)
        
    def sqrt_principal(self):
        r = self.norm()
        if r == 0:
            return Num3D(0, 0, 0)
        theta = math.atan2(math.hypot(self.x, self.y), self.z) / 2
        sqrt_r = math.sqrt(r)
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)
        if sin_t == 0:
            return Num3D(sqrt_r * cos_t, 0, 0)
        factor = sqrt_r * sin_t / math.hypot(self.x, self.y)
        return Num3D(sqrt_r * cos_t, factor * self.x, factor * self.y)

    def sqrt(self):
        from solusi import SolusiSet
        if self.x == 0 and self.y == 0 and self.z < 0:
            mag = math.sqrt(-self.z)
            def generator(u):
                return Num3D(0, mag * math.cos(u), mag * math.sin(u))
            return SolusiSet(generator=generator, params={'u': '[0, 2π)'}, rank=1)
        else:
            principal = self.sqrt_principal()
            return SolusiSet(elements=[principal, -principal])

    def __rshift__(self, other):
        """putaran siklik Z3: i >> j = 1"""
        other = Num3D.izin(other)
        z1, x1, y1 = self.z, self.x, self.y
        z2, x2, y2 = other.z, other.x, other.y

        z3 = z1*z2 + x1*y2 + y1*x2
        x3 = z1*x2 + x1*z2 + y1*y2
        y3 = z1*y2 + y1*z2 + x1*x2

        return Num3D(z3, x3, y3)
    
    # ============================================================
    # NORMA ALJABAR
    # ============================================================
    # ---------- GEOMETRI ----------
    def norm(self):
        return math.hypot(self.x, self.y, self.z)

    def dot(self, other):
        other = Num3D.izin(other)
        return self.z * other.z + self.x * other.x + self.y * other.y

    def angle(self, other):
        other = Num3D.izin(other)
        n1 = self.norm()
        n2 = other.norm()
        if n1 == 0 or n2 == 0:
            return 0.0
        cos_theta = self.dot(other) / (n1 * n2)
        cos_theta = max(-1.0, min(1.0, cos_theta))
        return math.acos(cos_theta)

    # ---------- FUNGSI TRIGONOMETRI ----------
    def sin(self):
        r = math.hypot(self.x, self.y)
        if r < 1e-12:
            return Num3D(math.sin(self.z), 0, 0)
        factor = math.cos(self.z) * math.sinh(r) / r
        return Num3D(
            math.sin(self.z) * math.cosh(r),
            factor * self.x,
            factor * self.y
        )

    def cos(self):
        r = math.hypot(self.x, self.y)
        if r < 1e-12:
            return Num3D(math.cos(self.z), 0, 0)
        factor = -math.sin(self.z) * math.sinh(r) / r
        return Num3D(
            math.cos(self.z) * math.cosh(r),
            factor * self.x,
            factor * self.y
        )

    def tan(self):
        return self.sin() / self.cos()

    # ---------- FUNGSI HIPERBOLIK ----------
    def sinh(self):
        r = math.hypot(self.x, self.y)
        if r < 1e-12:
            return Num3D(math.sinh(self.z), 0, 0)
        factor = math.cosh(self.z) * math.sin(r) / r
        return Num3D(
            math.sinh(self.z) * math.cos(r),
            factor * self.x,
            factor * self.y
        )

    def cosh(self):
        r = math.hypot(self.x, self.y)
        if r < 1e-12:
            return Num3D(math.cosh(self.z), 0, 0)
        factor = math.sinh(self.z) * math.sin(r) / r
        return Num3D(
            math.cosh(self.z) * math.cos(r),
            factor * self.x,
            factor * self.y
        )

    def tanh(self):
        return self.sinh() / self.cosh()

    @staticmethod
    def i(self:None) -> 'Num3D':
        """Sumbu (i, j) di bidang 2D."""
        if self is None:
            return Num3D(0, 1, 0)
        cosu = self.cos()
        sinu = self.sin()
        return Num3D(
            -(cosu.x + sinu.y),
            cosu.z,
            sinu.z)

    
    @staticmethod
    def j(self:None) -> 'Num3D':
        """Sumbu (j, i) di bidang 2D."""
        if self is None:
            return Num3D(0,0,1)
        cosu = self.cos()
        sinu = self.sin()
        return Num3D(
            -(cosu.y + sinu.x),
            sinu.z,
            cosu.z)
    
    # ============================================================
    # METRIK POSISI (EUCLIDEAN)
    # ============================================================
    def euclidean_distance(P1: 'Num3D', P2: 'Num3D') -> float:
        """Jarak antar dua titik di ruang 3D."""
        return (P1 - P2).norm()

    def centroid(points: list['Num3D']) -> 'Num3D':
        """Rata-rata posisi (centroid) untuk clustering."""
        total = Num3D(0, 0, 0)
        for p in points:
            total = total + p
        return total * (1.0 / len(points))

    # ============================================================
    # METRIK ROTASI (GEODESIK)
    # ============================================================
    def geodesic_distance(R1: 'Num3D', R2: 'Num3D') -> float:
        """Sudut rotasi antara dua rotor (dalam radian)."""
        # Pastikan rotor adalah unit (norma = 1)
        dot = R1.dot(R2)
        # Clamp untuk menghindari error numerik
        dot = max(-1.0, min(1.0, dot))
        return 2 * math.acos(dot)

    def karcher_mean(rotors: list['Num3D'], max_iter: int = 10) -> 'Num3D':
        """
        Rata-rata geodesik untuk sekumpulan rotor (orientasi).
        Algoritma: iterasi log-Euclidean.
        """
        from numr import NumR
        # Inisialisasi: rata-rata aritmetika lalu normalisasi
        mean = Num3D(0, 0, 0)
        for R in rotors:
            mean = mean + R
        mean = mean / mean.norm()  # proyeksi ke bola satuan

        for _ in range(max_iter):
            # 1. Hitung rata-rata logaritma (sum of vectors)
            sum_log = Num3D(0, 0, 0)
            for R in rotors:
                # diff = R * mean^{-1} (komposisi rotasi)
                diff = NumR.n3d(R ^ mean.inv())
                # Log dari diff menghasilkan vektor axis-angle (skalar part = 0)
                log_diff = diff.ln()
                sum_log = sum_log + log_diff
            
            # 2. Rata-rata vektor axis-angle
            avg_log = sum_log * (1.0 / len(rotors))
            
            # 3. Eksponensialkan untuk mendapatkan rotor baru
            mean = NumR.n3d(avg_log.exp() ^ mean)
            mean = mean / mean.norm()  # Normalisasi ulang
            
        return mean

    # ============================================================
    # METRIK GABUNGAN (POSISI + ROTASI)
    # ============================================================
    def combined_metric(A1:tuple['Num3D','Num3D'],
                        A2:tuple['Num3D','Num3D'],
                        alpha: float = 0.1) -> float:
        """
        X1, X2 adalah tuple (posisi_Num3D, rotor_Num3D)
        alpha: bobot untuk menyeimbangkan jarak posisi dan rotasi.
        """
        pos1, rot1 = A1
        pos2, rot2 = A2
        d_pos = Num3D.euclidean_distance(pos1, pos2)
        d_rot = Num3D.geodesic_distance(rot1, rot2)
        return d_pos + alpha * d_rot
    # ============================================================
    # ROTASI VEKTOR (RODRIGUES)
    # ============================================================
    @staticmethod
    def rotate_vector(v: 'Num3D', axis: 'Num3D', angle: float) -> 'Num3D':
        from numr import NumR #mengambil ^ dari quaternion
        axis = axis / axis.norm()
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        term1 = v * cos_a
        term2 = NumR.n3d(axis ^ v) * sin_a
        term3 = axis * (axis.dot(v) * (1 - cos_a))
        return term1 + term2 + term3

    # ---------- BASIS GEOMETRI {a, b, c} ----------
    @classmethod
    def from_abc(cls, A, B, C):
        sq6 = 1/ cls.SQRT6
        sq2 = cls.SQRT3 * sq6
        sq3 = cls.SQRT2 * sq6
        z = (A + B + C) * sq3
        x = (2 * A - B - C) *sq6
        y = (B - C) * sq2
        return cls(z, x, y)

    def to_abc(self):
        sq6 = 1/ self.SQRT6
        sq2 = self.SQRT3 * sq6
        sq3 = self.SQRT2 * sq6
        A = self.z * sq3 + 2 * self.x * sq6
        B = self.z * sq3 - self.x * sq6 + self.y * sq2
        C = self.z * sq3 - self.x * sq6 - self.y * sq2
        return A, B, C

    def abc_str(self, decimals=4):
        A, B, C = self.to_abc()
        parts = []
        if abs(A) > 1e-12:
            parts.append(f"{A:.{decimals}f} a")
        if abs(B) > 1e-12:
            parts.append(f"{B:.{decimals}f} b")
        if abs(C) > 1e-12:
            parts.append(f"{C:.{decimals}f} c")
        return " + ".join(parts) if parts else "0"
    
    # ---------- BASIS COLOUR {r, g, b} ----------
    @classmethod
    def from_rgb(cls, Mer, Hij, Bir):
        A = 2*Mer - 1
        B = 2*Hij - 1
        C = 2*Bir - 1
        return cls.from_abc(A, B, C)
    
    def to_rgb(self):
        A, B, C = self.to_abc()
        Mer = (A + 1) * 0.5
        Hij = (B + 1) * 0.5
        Bir = (C + 1) * 0.5
    return Mer, Hij, Bir

    def rgb_str(self, decimals=4):
        A, B, C = self.to_rgb()
        parts = []
        if abs(A) > 1e-12:
            parts.append(f"{A:.{decimals}f} r")
        if abs(B) > 1e-12:
            parts.append(f"{B:.{decimals}f} g")
        if abs(C) > 1e-12:
            parts.append(f"{C:.{decimals}f} b")
        return " + ".join(parts) if parts else "0"

    @classmethod
    def identity(cls):
        return cls(1, 0, 0)

    def __eq__(self, other):
        other = Num3D.izin(other)
        if not isinstance(other, Num3D):
            return False
        return (abs(self.z - other.z) < 1e-12 and
                abs(self.x - other.x) < 1e-12 and
                abs(self.y - other.y) < 1e-12)

    def __repr__(self):
        return f"Num3D({self.z:.4f}, {self.x:.4f}, {self.y:.4f})"

    def __str__(self):
        parts = []
        if abs(self.z) > 1e-12:
            parts.append(f"{self.z:.3f} ")
        if abs(self.x) > 1e-12:
            parts.append(f"{self.x:.3f} i ")
        if abs(self.y) > 1e-12:
            parts.append(f"{self.y:.3f} j")
        return " + ".join(parts) if parts else "0"

    def slerp(self, other: 'Num3D', t: float) -> 'Num3D':
        """
        Interpolasi linier spherical (Slerp) antara dua rotor.
        self dan other harus memiliki norma 1 (rotor unit).
        t = 0 → self, t = 1 → other
        """
        other = Num3D.izin(other)
        # Dot product
        dot = self.dot(other)
        # Clamp untuk menghindari error numerik
        dot = max(-1.0, min(1.0, dot))
        theta = math.acos(dot)

        if abs(theta) < 1e-12:
            return self  # tidak ada rotasi

        sin_theta = math.sin(theta)
        w1 = math.sin((1 - t) * theta) / sin_theta
        w2 = math.sin(t * theta) / sin_theta

        # Kombinasi linear (penjumlahan + perkalian skalar)
        return self * w1 + other * w2

# ---------- FUNGSI BANTUAN ----------
def num3d(z=0.0, x=0.0, y=0.0):
    return Num3D(z, x, y)

@staticmethod
def i():
    return Num3D(0, 1, 0)
@staticmethod
def j():
    return Num3D(0, 0, 1)

def mer():
    return Num3D.from_rgb(1, 0, 0)
def hij():
    return Num3D.from_rgb(0, 1, 0)
def bir():
    return Num3D.from_rgb(0, 0, 1)
def a():
    return Num3D.from_abc(1, 0, 0)
def b():
    return Num3D.from_abc(0, 1, 0)
def c():
    return Num3D.from_abc(0, 0, 1)
