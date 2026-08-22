"""
main.py - Tutorial penggunaan Num3D dan SolusiSet
"""
from num3d import Num3D, num3d, i, j
from solusi import SolusiSet, NumL, NumA, NumV
import math

def main():
    print("=== Tutorial Num3D ===")
    v1 = Num3D(1, 2, 3)
    v2 = num3d(0, -1, 4)
    print("v1 =", v1)
    print("v2 =", v2)

    print("\nv1 + v2 =", v1 + v2)
    print("v1 * v2 =", v1 * v2)
    print("v1 / v2 =", v1 / v2)
    print("v1 ** 2 =", v1 ** 2)
    print("exp(v1) =", v1.exp())
    print("ln_principal(v1) =", v1.ln_principal())

    print("\n=== Solusi Non-Injektif ===")
    neg = Num3D(-4, 0, 0)
    sqrt_sol = neg.sqrt()
    print("sqrt(-4) ->", sqrt_sol)
    print("  rank:", sqrt_sol.rank)
    print("  canonical:", sqrt_sol.canonical())
    print("  select(u=pi/2):", sqrt_sol.select(u=math.pi/2))
    sqrt_sol.set_default(u=math.pi)
    print("  setelah set_default(u=pi):", sqrt_sol.canonical())

    log_neg = Num3D(-2, 0, 0).ln()
    print("\nln(-2) ->", log_neg)
    print("  rank:", log_neg.rank)
    print("  canonical:", log_neg.canonical())
    print("  select(k=1, u=pi):", log_neg.select(k=1, u=math.pi))

    print("\n=== Operasi dengan SolusiSet ===")
    v3 = Num3D(1, 0, 0)
    print("sqrt(-4) + (1,0,0) =", sqrt_sol + v3)
    print("(1,0,0) + sqrt(-4) =", v3 + sqrt_sol)
    print("sqrt(-4) * 2 =", sqrt_sol * 2)

    print("\n=== Rotasi Vektor ===")
    axis = Num3D(0, 0, 1)
    vector = Num3D(0, 1, 0)
    angle = math.pi / 2
    rotated = Num3D.rotate_vector(vector, axis, angle)
    print(f"Rotasi {vector} sebesar {angle} rad mengelilingi {axis} => {rotated}")

    print("\n=== Sampling Solusi ===")
    samples = sqrt_sol.sample(n=8)
    print("Sampel sqrt(-4):")
    for s in samples:
        print("  ", s)

if __name__ == "__main__":
    main()