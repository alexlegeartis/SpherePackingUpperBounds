"""Run the bounded 6D warm-start discovery experiment; never certify output."""
import argparse
from pathlib import Path
from search_exact_guided_v2 import solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=Path,
                        default=Path(__file__).resolve().parent / 'seed.json')
    parser.add_argument('--output', type=Path, default=Path('candidate_new.json'))
    parser.add_argument('--rounds', type=int, default=5)
    args = parser.parse_args()
    if args.rounds < 1:
        parser.error('--rounds must be positive')
    if args.output.exists():
        parser.error('Output already exists; choose a new path to preserve prior results')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    solve(n=6, degree=31, nodes=5, angular=12, ratio=2.25,
          subdivisions=2, rounds=args.rounds, pairwise=False,
          output=str(args.output.resolve()), exclusion=1.0, counts=True,
          seed=str(args.seed.resolve()))
    print('NUMERICAL_ONLY: complete exact verification is still required.')


if __name__ == '__main__':
    main()
