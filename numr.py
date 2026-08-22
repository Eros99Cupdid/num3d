import math

class NumR:
    """
    Quaternion 4D: w*I + x*i + y*j + z
    Digunakan untuk rotasi 3D yang akurat.
    """

    def __init__(self, w: float = 1.0, x: float = 0.0, y: float = 0.0, z: float = 0.0):
        self.w = float(w)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    # ============================================================
    # OPERATOR QUATERNION
    # ============================================================
    @staticmethod
    def izin(other):
        """Periksa apakah objek lain adalah NumR atau Num3D."""
        from num3d import Num3D
        if isinstance(other, NumR):
            return other
        elif isinstance(other, Num3D):
            return NumR(0, other.x, other.y, other.z)
        elif isinstance(other, (int, float)):
            return NumR(0, 0, 0, other)
        else:
            raise TypeError(f"Tipe {type(other)} tidak didukung")

    def __mul__(self, other):
        from num3d import Num3D
        other = NumR.izin(other)
        v2 = Num3D(other.z, other.x, other.y)
        v1 = Num3D(self.z, self.x, self.y)
        return v1 * v2
        
    def __rmul__(self, other):
        return self.__mul__(other)

    def conjugate(self):
        return NumR(self.w, -self.x, -self.y, -self.z)

    def norm_sq(self) -> float:
        return self.w*self.w + self.x*self.x + self.y*self.y + self.z*self.z

    def norm(self) -> float:
        return math.sqrt(self.norm_sq())

    def inverse(self) -> 'NumR':
        """Invers quaternion (untuk unit, sama dengan konjugasi)"""
        n = self.norm_sq()
        if n == 0:
            raise ZeroDivisionError("Zero quaternion tidak memiliki invers")
        return NumR(self.w/n, -self.x/n, -self.y/n, -self.z/n)

    def normalize(self) -> 'NumR':
        n = self.norm()
        if n == 0:
            return NumR(1,0,0,0)
        return NumR(self.w/n, self.x/n, self.y/n, self.z/n)

    # ============================================================
    # ROTASI VEKTOR
    # ============================================================

    def rotate(self, other):
        """
        Menerapkan rotasi ke vektor Num3D.
        Num3D dianggap sebagai (z, x, y) di mana z adalah komponen sumbu 1,
        x adalah komponen i, y adalah komponen j.
        """
        from num3d import Num3D
        v = NumR.izin(other)
        v_rot = self ^ v ^ self.inverse()
        return Num3D(v_rot.z, v_rot.x, v_rot.y)

    def rotate_inv(self, other):
        """Rotasi dengan invers quaternion: q^{-1} * v * q"""
        from num3d import Num3D
        v = NumR.izin(other)
        v_rot = self.inverse() ^ v ^ self
        return Num3D(v_rot.z, v_rot.x, v_rot.y)

    # ============================================================
    # KONSTRUKSI DARI SUMBU & SUDUT
    # ============================================================

    @classmethod
    def from_axis_angle(cls, axis, theta: float) -> 'NumR':
        """
        Membuat quaternion rotasi dari sumbu (Num3D) dan sudut (radian).
        axis tidak harus unit, akan dinormalisasi.
        """
        n = axis / axis.norm()
        half = theta / 2.0
        z = math.sin(half)
        return cls(math.cos(half), z * n.x, z * n.y, z * n.z)

    # ============================================================
    # REPRESENTASI
    # ============================================================

    def __repr__(self):
        return f"NumR({self.w:.4f}, {self.x:.4f}, {self.y:.4f}, {self.z:.4f})"

    def __str__(self):
        parts = []
        if self.w != 0: parts.append(f"{self.w:.3f} I")
        if self.x != 0: parts.append(f"{self.x:.3f} i")
        if self.y != 0: parts.append(f"{self.y:.3f} j")
        if self.z != 0: parts.append(f"{self.z:.3f}")
        return " + ".join(parts) if parts else "0"

    def __matmul__(self, other):
        return self.rotate(other)

    def __rmatmul__(self, other):
        return self.rotate_inv(other)
    
    def __add__(self, other):
        if isinstance(other, NumR):
            return NumR(self.w + other.w, self.x + other.x, self.y + other.y, self.z + other.z)
        raise TypeError(f"+ tidak didukung antara NumR dan {type(other)}")
    @staticmethod
    def n3d(other):
        """Konversi ke Num3D (hanya bagian vektor)"""
        from num3d import Num3D
        return Num3D(other.z, other.x, other.y)

    def __xor__(self, other: 'NumR') -> 'NumR':
        """Perkalian quaternion (Hamilton product) dengan operator ^."""
        from num3d import Num3D
        # Perkalian quaternion standar (Hamilton)
        if isinstance(other, Num3D):
            other = NumR(0, other.x, other.y, other.z)
            return self ^ other
        elif isinstance(other, NumR):
            w1, x1, y1, z1 = self.w, self.x, self.y, self.z
            w2, x2, y2, z2 = other.w, other.x, other.y, other.z
            return NumR(
                w1*w2 - x1*x2 - y1*y2 - z1*z2,
                w1*x2 + x1*w2 + y1*z2 - z1*y2,
                w1*y2 - x1*z2 + y1*w2 + z1*x2,
                w1*z2 + x1*y2 - y1*x2 + z1*w2)
        elif isinstance(other, (int, float)):
            # Perkalian skalar
            return self ^ NumR(other,0,0,0)
        else:
            raise TypeError(f"* tidak didukung antara NumR dan {type(other)}")

