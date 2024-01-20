"""
2024 1 17
Jinpeng Zhai 翟锦鹏
pso.py: Particle Swarm Optimization
"""
import inspect
from random import random  # random [0, 1)
from multiprocessing import Pool
from math import inf, isinf


def pso(
    f,
    ss,
    n=20,
    iw=0.8,
    cog=0.1,
    soc=0.1,
    x_rel_tol=0,
    y_rel_tol=0,
    y_abs_tol=1e-14,
    y=None,
    consecutive=25,
    its=1000,
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
            pass
        else:
            raise ValueError(
                f"Each element wise constraint must be of form (min, max), not {str(s)}."
            )

    itc = 0  # iterations counter
    with Pool() as pool:
        # print("pool started")
        vps = []
        vfps = []
        while len(vps) < n:
            ps = []
            for _ in range(n):
                p = []
                for s in ss:
                    s_min, s_max = min(s), max(s)
                    delta = s_max - s_min
                    p.append(random() * delta + s_min)

                ps.append(p)

            fps = pool.starmap(f, ps)

            for p, fp in zip(ps, fps):
                if not isinf(fp):
                    vps.append(p)
                    vfps.append(fp)

        ps = vps[:n]
        fps = vfps[:n]
        pbs = [v for v in ps]
        fpbs = [v for v in fps]
        vs = [[random() * (max(s) - min(s)) + min(s) for s in ss] for _ in range(n)]

        fgb = min(fpbs)
        gb = pbs[fpbs.index(fgb)]  # initialize global best
        # print("generated valid initial guesses")

        it = 0
        while it < its:
            # updates global best
            new_ps = []
            new_vs = []

            r_cog, r_soc = random(), random()

            for p, pb, v in zip(ps, pbs, vs):
                new_v = []
                for p_i, v_i, pb_i, gb_i in zip(p, v, pb, gb):
                    new_v.append(
                        iw * v_i
                        + r_cog * cog * (pb_i - p_i)
                        + r_soc * soc * (gb_i - p_i)
                    )
                new_vs.append(new_v)
                new_ps.append([p_i + v_i for p_i, v_i in zip(p, new_v)])

            ps = new_ps
            vs = new_vs

            fps = pool.starmap(f, ps)

            new_pbs = []
            new_fpbs = []

            for p, pb, fp, fpb in zip(ps, pbs, fps, fpbs):
                if fp < fpb:
                    new_pbs.append(p)
                    new_fpbs.append(fp)
                else:
                    new_pbs.append(pb)
                    new_fpbs.append(fpb)

            pbs = new_pbs
            fpbs = new_fpbs

            minfs = min(fps)
            if minfs < fgb:
                fgb = minfs
                gb = ps[fps.index(fgb)]
                itc = 0

            print("val", minfs, fgb)

            if y is not None:
                if abs(fgb - y) < y * y_rel_tol:
                    print("broken via rel", it)
                    break
                elif abs(fgb - y) < y_abs_tol:
                    print("broken via abs", it)
                    break

            if itc > consecutive:
                print("broken via cons", it)
                break
            else:
                components = [*zip(*ps)]
                delta = 0
                for component in components:
                    delta = max(
                        (max(component) - min(component))
                        / min(abs(v) for v in component),
                        delta,
                    )
                # print("delta", delta)
                if delta < x_rel_tol:
                    print("broken via x rel tol", it)
                    break

            itc += 1
            it += 1

        return gb, fgb


from math import sin, cos, pi


def f(x, y):
    return (x - 3.14) ** 2 + (y - 2.72) ** 2 + sin(3 * x + 1.41) + sin(4 * y - 1.73)


if __name__ == "__main__":
    print(pso(f, [[0, 5], [0, 5]]))
