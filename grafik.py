"""
grafik.py - Visualisasi solusi 3D menggunakan matplotlib
"""
from num3d import Num3D
from solusi import SolusiSet, NumL, NumA, NumV
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_solution(sol: SolusiSet, n=50, title="Solusi 3D"):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    if sol.rank == 0:
        if sol.elements:
            for el in sol.elements:
                ax.scatter(el.z, el.x, el.y, color='r', s=100, label='Solusi')
        else:
            pt = sol.canonical()
            ax.scatter(pt.z, pt.x, pt.y, color='r', s=100, label='Solusi')
    else:
        samples = sol.sample(n=n)
        zs = [p.z for p in samples]
        xs = [p.x for p in samples]
        ys = [p.y for p in samples]
        ax.scatter(zs, xs, ys, c='b', s=20, alpha=0.7, label=f'rank {sol.rank}')

    ax.set_xlabel('z (1)')
    ax.set_ylabel('x (i)')
    ax.set_zlabel('y (j)')
    ax.set_title(title)
    ax.legend()
    plt.show()

def plot_points(points, title="Titik Num3D"):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    zs = [p.z for p in points]
    xs = [p.x for p in points]
    ys = [p.y for p in points]
    ax.scatter(zs, xs, ys, c='g', s=50)
    ax.set_xlabel('z')
    ax.set_ylabel('x')
    ax.set_zlabel('y')
    ax.set_title(title)
    plt.show()

if __name__ == "__main__":
    sol_sqrt = Num3D(-4, 0, 0).sqrt()
    plot_solution(sol_sqrt, n=100, title="Akar dari -4 (lingkaran imajiner)")

    sol_ln = Num3D(-2, 0, 0).ln()
    plot_solution(sol_ln, n=200, title="ln(-2) (permukaan silinder)")

    line_sol = Num3D(1,1,0) / Num3D(0,1,0)
    plot_solution(line_sol, n=50, title="Garis solusi pembagian")