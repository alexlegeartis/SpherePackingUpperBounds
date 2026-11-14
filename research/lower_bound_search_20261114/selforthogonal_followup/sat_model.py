"""SAT search for a coordinate placement of a specified ternary fiber type."""

from pathlib import Path
import argparse, hashlib, itertools, json, time
import numpy as np
import pynauty
from pysat.solvers import Solver
from pysat.formula import CNF
from sympy.combinatorics import Permutation, PermutationGroup

P = Path(__file__).resolve().parent
from finite import parent_supports, shortened_supports


def source_data(G):
    n = G.shape[1]
    T = np.zeros((1, n), dtype=np.int16)
    for row in G:
        T = np.concatenate([(T + a * row) % 3 for a in range(3)])
    w = np.count_nonzero(T, axis=1)
    m = np.sum(
        (T != 0).astype(np.uint32) * (1 << np.arange(n, dtype=np.uint32)), axis=1
    )
    return [sorted(set(map(int, m[w == weight]))) for weight in (6, 9)]


def symmetries(n, six, nine):
    adjacency = {i: [] for i in range(n + len(six) + len(nine))}
    for index, support in enumerate(six + nine, n):
        for j in range(n):
            if support >> j & 1:
                adjacency[index].append(j)
                adjacency[j].append(index)
    graph = pynauty.Graph(
        len(adjacency),
        adjacency_dict=adjacency,
        vertex_coloring=[
            set(range(n)),
            set(range(n, n + len(six))),
            set(range(n + len(six), len(adjacency))),
        ],
    )
    raw = pynauty.autgrp(graph)[0]
    generators = [list(p[:n]) for p in raw]
    for permutation in generators:
        assert sorted(permutation) == list(range(n))
        for supports in (six, nine):
            assert {
                sum(1 << permutation[i] for i in range(n) if s >> i & 1)
                for s in supports
            } == set(supports)
    group = PermutationGroup([Permutation(p) for p in generators])
    order = int(group.order())
    orbits = []
    for i in range(n):
        orbits.append(sorted(group.orbit(i)))
        group = group.stabilizer(i)
    return generators, order, orbits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n", type=int, choices=(17, 18))
    parser.add_argument("index", type=int)
    parser.add_argument("--conflicts", type=int, default=20000)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument("--empty-target", action="store_true")
    parser.add_argument("--proof", action="store_true")
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()
    n = args.n
    ci = args.index
    start = time.time()
    G = np.load(P / "inputs/catalogs.npz")[str(n)][ci]
    six, nine = source_data(G)
    generators, order, orbits = symmetries(n, six, nine)
    parent, _, _ = parent_supports()
    current = shortened_supports(parent, list(range(n, 20)))
    S6 = sorted(current[6])
    S7 = sorted(current[7])
    S9 = set(current[9])
    if args.empty_target:
        S6 = []
        S7 = []
        S9 = set()
    bad6 = set(S6)
    for s in S7:
        for i in range(n):
            if s >> i & 1:
                bad6.add(s ^ (1 << i))
    path = (
        args.output_prefix
        or P / (f"work/sat{n}_{ci}" + ("_control" if args.empty_target else ""))
    ).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    clauses = []
    nextvar = n * n
    membership = {}
    seen9 = set()
    iterations = 0
    status = "UNKNOWN"
    reason = "iteration cap"

    def x(i, j):
        return 1 + i * n + j

    with Solver(name="g4", with_proof=args.proof) as solver:

        def add(clause):
            clauses.append(clause)
            solver.add_clause(clause)

        def member(F):
            nonlocal nextvar
            if F in membership:
                return membership[F]
            v = list(range(nextvar + 1, nextvar + n + 1))
            nextvar += n
            membership[F] = v
            for i in range(n):
                for j in range(n):
                    add([-x(i, j), v[i] if F >> j & 1 else -v[i]])
            return v

        for i in range(n):
            add([x(i, j) for j in range(n)])
            add([x(j, i) for j in range(n)])
            for a, b in itertools.combinations(range(n), 2):
                add([-x(i, a), -x(i, b)])
                add([-x(a, i), -x(b, i)])
        # Lexicographically least image tuple in each source-automorphism orbit.
        for i, orbit in enumerate(orbits):
            for j in orbit:
                if j == i:
                    continue
                assert j > i
                for a in range(n):
                    for b in range(a):
                        add([-x(i, a), -x(j, b)])
        for F in S6 + S7:
            v = member(F)
            for s in six:
                add([-v[i] for i in range(n) if s >> i & 1])
        print(
            "BUILT",
            n,
            ci,
            "source_supports",
            len(six),
            len(nine),
            "autorder",
            order,
            "vars",
            nextvar,
            "clauses",
            len(clauses),
            "seconds",
            time.time() - start,
            flush=True,
        )
        for iteration in range(args.iterations):
            iterations = iteration + 1
            solver.conf_budget(args.conflicts)
            sat = solver.solve_limited(expect_interrupt=True)
            if sat is None:
                reason = "conflict budget"
                break
            if not sat:
                status = "UNSAT"
                reason = "solver"
                break
            model = set(v for v in solver.get_model() if v > 0)
            permutation = [
                next(j for j in range(n) if x(i, j) in model) for i in range(n)
            ]
            assert len(set(permutation)) == n

            def image(s):
                return sum(1 << permutation[i] for i in range(n) if s >> i & 1)

            assert not any(image(s) in bad6 for s in six)
            violations = [(s, image(s)) for s in nine if image(s) in S9]
            if not violations:
                status = "SAT"
                reason = "direct support audit"
                H = np.zeros_like(G)
                for i in range(n):
                    H[:, permutation[i]] = G[:, i]
                Path(str(path) + ".candidate.txt").write_text(
                    f"{n} 8 0 0\n"
                    + "\n".join(" ".join(map(str, row)) for row in H)
                    + "\n"
                )
                break
            for s, F in violations:
                assert (s, F) not in seen9
                seen9.add((s, F))
                v = member(F)
                add([-v[i] for i in range(n) if s >> i & 1])
            if iteration % 10 == 0:
                print(
                    "ITERATION",
                    iteration,
                    "violations",
                    len(violations),
                    "vars",
                    nextvar,
                    "clauses",
                    len(clauses),
                    "seconds",
                    time.time() - start,
                    flush=True,
                )
        stats = solver.accum_stats()
        if status == "UNSAT" and args.proof:
            CNF(from_clauses=clauses).to_file(str(path) + ".cnf")
            proof = solver.get_proof()
            Path(str(path) + ".drat").write_text("\n".join(proof) + "\n")
    result = {
        "n": n,
        "catalog_index": ci,
        "empty_target_control": args.empty_target,
        "proof_emitted": status == "UNSAT" and args.proof,
        "result": status,
        "reason": reason,
        "iterations": iterations,
        "variables": nextvar,
        "clauses": len(clauses),
        "automorphism_group_order": order,
        "support_automorphism_generators": generators,
        "stabilizer_orbits": orbits,
        "seconds": time.time() - start,
        "solver_stats": stats,
    }
    Path(str(path) + ".json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in result.items()
                if k not in ("support_automorphism_generators", "stabilizer_orbits")
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
