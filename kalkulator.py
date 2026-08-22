"""
Kalkulator interaktif untuk Num3D
"""
from num3d import Num3D, num3d, i, j
from solusi import SolusiSet, NumL, NumA, NumV
import math

def _sqrt(x):
    if isinstance(x, Num3D):
        return x.sqrt()
    if x < 0:
        return Num3D(x, 0, 0).sqrt()
    return math.sqrt(x)

def _ln(x):
    if isinstance(x, Num3D):
        return x.ln()
    if x < 0:
        return Num3D(x, 0, 0).ln()
    if x == 0:
        raise ValueError("ln(0) undefined")
    return Num3D(math.log(x), 0, 0)  # representasi scalar sebagai Num3D

def _exp(x):
    if isinstance(x, Num3D):
        return x.exp()
    return Num3D(math.exp(x), 0, 0)

def _sin(x):
    if isinstance(x, Num3D):
        return x.sin()
    return Num3D(math.sin(x), 0, 0)

def _cos(x):
    if isinstance(x, Num3D):
        return x.cos()
    return Num3D(math.cos(x), 0, 0)

def _tan(x):
    if isinstance(x, Num3D):
        return x.tan()
    return Num3D(math.tan(x), 0, 0)

def eval_expr(expr):
    namespace = {
        'Num3D': Num3D,
        'num3d': num3d,
        'i': i(),
        'j': j(),
        'sqrt': _sqrt,
        'ln': _ln,
        'ln_principal': lambda x: x.ln_principal() if isinstance(x, Num3D) else Num3D(math.log(x), 0, 0) if x > 0 else (Num3D(math.log(-x), math.pi, 0) if x < 0 else None),
        'exp': _exp,
        'sin': _sin,
        'cos': _cos,
        'tan': _tan,
        'pi': math.pi,
        'e': math.e,
    }
    try:
        result = eval(expr, {"__builtins__": {}}, namespace)
        return result
    except Exception as e:
        return f"Error: {e}"

def main():
    print("Kalkulator Num3D")
    print("Contoh: 1 + 2*i + 3*j, sqrt(-4), ln(-2), exp(i*pi)")
    print("Ketik 'exit' untuk keluar.\n")
    while True:
        expr = input(">>> ").strip()
        if expr.lower() in ['exit', 'quit']:
            break
        if not expr:
            continue
        result = eval_expr(expr)
        if isinstance(result, SolusiSet):
            print("SolusiSet:")
            print("  rank:", result.rank)
            print("  canonical:", result.canonical())
            print("  params:", result.params)
            if result.rank > 0:
                samples = result.sample(5)
                print("  sampel:", samples)
        else:
            print(result)

if __name__ == "__main__":
    main()