from __future__ import annotations

import argparse
import json
from pathlib import Path

from scipy.io import loadmat


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ONLPoisson Python run entry")
    parser.add_argument("--images", default="images.mat", help="Path to MATLAB .mat image file")
    parser.add_argument("--output-dir", default="results", help="Directory for run artifacts")
    return parser


def main() -> None:
    args = _build_parser().parse_args()

    images_path = Path(args.images)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mat_data = loadmat(images_path)
    image_keys = sorted([k for k in mat_data.keys() if not k.startswith("__")])

    payload = {
        "images_path": str(images_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "num_image_keys": len(image_keys),
        "image_keys": image_keys,
        "shapes": {
            key: list(getattr(mat_data[key], "shape", ()))
            for key in image_keys
        },
        "note": ".m/.c/.mexw64/.dll are historical references and not in the Python execution path.",
    }

    out_file = output_dir / "run_info.json"
    out_file.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"Wrote: {out_file}")


if __name__ == "__main__":
    main()
