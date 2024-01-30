"""
coordinate descent in multiple dimensions
"""
import inspect
from num import gss
from random import random


def coordesc(f, ss, x_rel_tol=0, y_rel_tol=1e-5, y_abs_tol=0, its=100, guess=10):
    """
    Coordinate descent in the positive real
    """
    sig = inspect.signature(f)
    paramstrs = [
        str(param)
        for param in sig.parameters.values()
        if param.kind == param.POSITIONAL_OR_KEYWORD
    ]
    n = len(paramstrs)

    if len(ss) != n:
        raise ValueError(
            f"Mismatch in argument length for searched function f({','.join(paramstrs)})"
            + " and initial guess being supplied."
        )

    print("guessing initials")
    for i in range(guess):
        args = []
        for s in ss:
            s_min, s_max = min(s), max(s)
            delta = s_max - s_min
            args.append(random() * delta + s_min)
        try:
            y = f(*args)
            break
        except ValueError:
            pass

    if i == guess - 1:
        raise ValueError(f"Maximum guesses ({guess}) exceeded.")

    def g(x, i):
        return f(*[v if i != j else x for j, v in enumerate(args)])

    def findBound(func, x_probe, x_bound, tol, record=[]):
        # we trust that the value passed in results in a valid evaluation
        x_valid = x_probe
        delta = x_bound - x_probe
        up = x_bound > x_probe

        while abs(2 * delta) > tol and (
            x_probe <= x_bound if up else x_probe >= x_bound
        ):
            try:
                record.append((x_probe, func(x_probe)))
                x_valid = x_probe
            except ValueError:
                delta *= 0.5
            finally:
                x_probe = x_valid + delta

        return x_valid

    print(f"initial value found at f{args}={y}")

    for it in range(its):
        x_delta = 0
        for i in range(n):
            record = []
            x_min, x_max = min(ss[i]), max(ss[i])
            tol = (x_max - x_min) * x_rel_tol

            x_up = findBound(lambda x: g(x, i), args[i], x_max, tol, record)
            x_lo = findBound(lambda x: g(x, i), args[i], x_min, tol, record)

            x_best = 0.5 * sum(
                gss(
                    lambda x: g(x, i),
                    x_lo,
                    x_up,
                    findMin=True,
                    x_tol=tol,
                    y_abs_tol=y_abs_tol,
                    y_rel_tol=y_rel_tol,
                )
            )

            x_delta += abs(x_best - args[i])
            args[i] = x_best

        y_star = f(*args)
        print(f"it = {it}")
        print(f"f{args} = {y_star}")
        if y is not None and abs(y_star - y) < max(y_abs_tol, abs(y_star) * y_rel_tol):
            break

        y = y_star

    if it == its - 1:
        raise ValueError(f"Maximum iterations exceeded at it={it}")
    else:
        return args, y_star


from math import sin, cos


def f(x, y):
    return (x - 3.14) ** 2 + (y - 2.72) ** 2 + sin(3 * x + 1.41) + sin(4 * y - 1.73)


if __name__ == "__main__":
    print(coordesc(f, [(0, 5), (0, 5)]))
