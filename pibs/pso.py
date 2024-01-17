"""
2024 1 17
Jinpeng Zhai 翟锦鹏
pso.py: Particle Swarm Optimization
"""
import inspect
from random import random  # random [0, 1)
from multiprocessing import Pool
from math import inf


def pso(
    f,
    ss,
    n=20,
    iw=0.8,
    cog=0.1,
    soc=0.1,
    y_rel_tol=0,
    y_abs_tol=1e-14,
    consecutive=100,
):
    """
    Arguments:
    f   : function to be PSO optimized. Should accept argument for each element in ss.
    ss  : search space. List of tuple, where each tuple specifies the element wise
        : (min, max)
    n   : number of particles to initialize for searching.
    iw  : inertia weight  parameter
    cog : local minima attraction parameter, "cognitive"
    soc : global minima attraction parameter, "social"
    rel_err
        : relative error by which all particles within a swarm must conform to before
        : the solution is accepted.
    abs_err
        : absolute error by which all particles within a swarm must conform to before
        : the solution is accepted.
    consecutive
        : number of iterations where the global value is not decreasing before the
        : solution is accepted.
    """

    sig = inspect.signature(f)
    paramstrs = [
        str(param)
        for param in sig.parameters.values()
        if param.kind == param.POSITIONAL_OR_KEYWORD
    ]

    if len(paramstrs) != len(ss):
        raise ValueError(
            f"Mismatch in argument length for searched function f({','.join(paramstrs)})"
            + " and search space constraint."
        )

    ps = [
        [] for _ in range(n)
    ]  # initialize positional vector of the searching particles.
    for s in ss:
        """
        generate initial position vector using random distribution,
        component by component
        """
        if (
            (isinstance(s, tuple) or isinstance(s, list))
            and len(s) == 2
            and s[0] != s[1]
        ):
            s_min, s_max = min(s), max(s)
            delta = s_max - s_min
            for p in ps:
                p.append(random() * delta + s_min)
        else:
            raise ValueError(
                f"Each element wise constraint must be of form (min, max), not {str(s)}."
            )

    pbs = [None for _ in range(n)]  # location of particle best
    fpbs = [inf for _ in range(n)]  # value of particle best
    vs = [[0 for _ in ss] for _ in range(n)]

    fgb = inf
    gb = None  # initialize global best

    itc = 0  # iterations counter
    n = 0
    with Pool() as pool:
        while True:
            fs = pool.starmap(f, ps)

            # updates global best
            minfs = min(fs)
            if minfs < fgb:
                fgb = minfs
                gb = ps[fs.index(fgb)]
                itc = 0

            maxfs = max(fs)
            if abs(maxfs - minfs) < min(abs(maxfs), abs(minfs)) * y_rel_tol:
                print("broken via rel", n)
                break
            elif abs(maxfs - minfs) < y_abs_tol:
                print("broken via abs", n)
                break
            elif itc > consecutive:
                print("broken via cons", n)
                break

            new_ps = []
            new_pbs = []
            new_fpbs = []
            new_vs = []

            for p, fp, pb, fpb, v in zip(ps, fs, pbs, fpbs, vs):
                if fp < fpb:
                    pb = p
                    fpb = fp

                new_pbs.append(pb)
                new_fpbs.append(fpb)
                new_ps.append([p_i + v_i for p_i, v_i in zip(p, v)])

                new_v = []
                for p_i, v_i, pb_i, gb_i in zip(p, v, pb, gb):
                    r_cog, r_soc = random(), random()
                    new_v.append(
                        iw * v_i
                        + r_cog * cog * (pb_i - p_i)
                        + r_soc * soc * (gb_i - p_i)
                    )

                new_vs.append(new_v)

            ps = new_ps
            pbs = new_pbs
            fpbs = new_fpbs
            vs = new_vs

            itc += 1

            n += 1

        return gb, fgb


from math import sin, cos, pi


def f(x, y):
    return (x - 3.14) ** 2 + (y - 2.72) ** 2 + sin(3 * x + 1.41) + sin(4 * y - 1.73)


if __name__ == "__main__":
    print(pso(f, [[0, 5], [0, 5]]))
