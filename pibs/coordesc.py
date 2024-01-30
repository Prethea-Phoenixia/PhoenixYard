"""
coordinate descent in multiple dimensions
"""
import inspect
from num import gss
from random import random
import sys, traceback


def coordesc(
    f,
    ss,
    x_rel_tol=1e-6,
    y_ref=0,
    y_rel_tol=1e-5,
    y_abs_tol=0,
    its=100,
    guess=100,
    debug=False,
):
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

    y_abs_tol = max(y_abs_tol, abs(y_ref) * y_rel_tol)

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
        except ValueError as e:
            if debug:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                errMsg = "".join(
                    traceback.format_exception(exc_type, exc_value, exc_traceback)
                )
                print(str(errMsg))
            pass

    if i == guess - 1:
        raise ValueError(f"Maximum guesses ({guess}) exceeded.")

    def g(x, i):
        return f(*[v if i != j else x for j, v in enumerate(args)])

    def findBound(func, x_probe, x_bound, tol, record=[]):
        x_valid = None
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

            # print("x_valid", x_valid)
            # print("x_probe", x_probe)
            # print("x_bound", x_bound)
            # print("delta", delta, "tol", tol)
            # input()

        return x_valid

    print(f"initial value found at f{args}={y}")
    y = None
    for it in range(its):
        x_delta = 0
        for i in range(n):
            record = []
            x_min, x_max = min(ss[i]), max(ss[i])
            tol = (x_max - x_min) * x_rel_tol

            x_up = findBound(lambda x: g(x, i), args[i], x_max, tol, record)
            x_lo = findBound(lambda x: g(x, i), args[i], x_min, tol, record)

            print(f"operating on arg {i} between {x_lo} and {x_up}")

            x_best = 0.5 * sum(
                gss(
                    lambda x: g(x, i),
                    x_lo,
                    x_up,
                    findMin=True,
                    y_abs_tol=y_abs_tol,
                    y_rel_tol=y_rel_tol,
                )
            )
            print(f"x_best found at {x_best}")

            x_delta += abs(x_best - args[i])
            args[i] = x_best

        y_star = f(*args)
        print(f"it = {it}")
        print(f"f{args} = {y_star}")

        if y is not None and abs(y_star - y_ref) < y_abs_tol:
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
