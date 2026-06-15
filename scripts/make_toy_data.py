"""Generate deterministic Level A toy data for DarkDNA-Observer."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from darkdna.toy_data import make_toy_data


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate deterministic DarkDNA-Observer toy data.")
    parser.add_argument("--out", "--outdir", dest="out", default="data/toy", type=Path)
    parser.add_argument("--seed", default=42, type=int)
    args = parser.parse_args()
    paths = make_toy_data(args.out, seed=args.seed)
    for key, path in paths.items():
        print(f"{key}\t{path}")


if __name__ == "__main__":
    main()
