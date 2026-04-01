#!/usr/bin/env python3
"""Unified Python entrypoint to reproduce ONL Table-1 results."""

from __future__ import annotations

import argparse
import json

from onlpoisson import run_table1_ours


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--images-mat", default="images.mat", help="Path to images.mat")
    parser.add_argument("--repeats", type=int, default=5, help="Number of Monte-Carlo repeats")
    parser.add_argument("--seed", type=int, default=0, help="Base RNG seed")
    args = parser.parse_args()

    result = run_table1_ours(args.images_mat, repeats=args.repeats, seed=args.seed)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
