#!/usr/bin/env python3
"""
Standalone CLI: derive membrane-safe binder-target residues from a PPM3.0
oriented structure ("...out.pdb"), for use directly with the
--contact_residues flag in boltz_ph/design.py or chai_ph/design.py.

Example:
    python utils/ppm_contact_residues.py \\
        --ppm_pdb /path/to/fold_bama_mono_model_0out.pdb

    python utils/ppm_contact_residues.py \\
        --ppm_pdb /path/to/fold_bama_mono_model_0out.pdb \\
        --exclude_pore_facing --out bama_contact_residues.txt
"""

import argparse
import os
import sys

sys.path.append(os.path.dirname(__file__))
from membrane_utils import (
    classify_residues,
    get_membrane_safe_residues,
    summarize_classification,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Derive membrane-avoiding contact residues from a PPM3.0-oriented PDB."
    )
    parser.add_argument("--ppm_pdb", required=True, type=str, help="Path to a PPM3.0 '...out.pdb' file.")
    parser.add_argument(
        "--chain_id",
        default=None,
        type=str,
        help="Chain ID to analyze (default: auto-detect if the structure has a single chain).",
    )
    parser.add_argument(
        "--exclude_pore_facing",
        action="store_true",
        help="Exclude barrel-lumen/pore-facing residues; keep only residues outside the membrane.",
    )
    parser.add_argument("--out", default="", type=str, help="Optional path to also write the residue list to.")
    return parser.parse_args()


def main():
    args = parse_args()
    classification = classify_residues(args.ppm_pdb, chain_id=args.chain_id)
    print(summarize_classification(classification), file=sys.stderr)

    residues = get_membrane_safe_residues(
        args.ppm_pdb,
        chain_id=args.chain_id,
        include_pore_facing=not args.exclude_pore_facing,
    )
    residue_str = ",".join(str(r) for r in residues)
    print(residue_str)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write(residue_str + "\n")
        print(f"Wrote {len(residues)} residues to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
