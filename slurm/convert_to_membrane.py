#!/usr/bin/env python3
"""
Converts an existing (regular Protein Hunter) SLURM launch script into a
membrane-avoiding one: renames the conda env, tags the job/run name, and
inserts --ppm_pdb/--ppm_chain_id into the design.py call. Everything else
(resources, --protein_seqs, --percent_X, AF3 paths, etc.) is left untouched.

Usage:
    python slurm/convert_to_membrane.py \\
        --input /path/to/run_ph_SecE_CLIBASIA_00140_R1_NOAFVAL_mini_binders_200_6.sh \\
        --ppm_pdb fold_sece_mono_model_0out.pdb \\
        [--ppm_chain_id A] \\
        [--output slurm/run_ph_SecE_..._membrane.sh]
"""

import argparse
import os
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, type=str, help="Existing run_ph_*.sh script to convert.")
    parser.add_argument(
        "--ppm_pdb",
        required=True,
        type=str,
        help="Filename (not path) of the PPM3.0 '...out.pdb' this target, expected under ppm_structures/.",
    )
    parser.add_argument("--ppm_chain_id", default="A", type=str, help="Chain ID in --ppm_pdb (default: A).")
    parser.add_argument(
        "--output",
        default="",
        type=str,
        help="Output path. Default: alongside --input, with _membrane inserted before .sh.",
    )
    return parser.parse_args()


def default_output_path(input_path: str) -> str:
    base, ext = os.path.splitext(input_path)
    return f"{base}_membrane{ext}"


def convert(lines, ppm_pdb: str, ppm_chain_id: str):
    out = []
    warnings = []
    saw_conda_activate = False
    saw_name_flag = False
    inserted_ppm_flags = False

    for line in lines:
        stripped = line.strip()

        if any(flag in stripped for flag in ("--ppm_pdb", "--ppm_chain_id")):
            raise ValueError("Input script already contains --ppm_pdb/--ppm_chain_id -- looks like it's already a membrane script.")

        # 1. conda env rename (exact line match, so we don't clobber an
        #    already-forked env name like proteinhunter_DO)
        if stripped == "conda activate proteinhunter":
            indent = line[: len(line) - len(line.lstrip())]
            out.append(f"{indent}conda activate proteinhunter_membrane\n")
            saw_conda_activate = True
            continue

        # 2. SBATCH job name: tag with membrane_
        m = re.match(r"^(#SBATCH -J )protein_hunter_(.*)$", stripped)
        if m:
            out.append(f"{m.group(1)}protein_hunter_membrane_{m.group(2)}\n")
            continue

        # 3. --name <value>: append _membrane_avoiding to the value token only
        m = re.match(r"^(\s*--name\s+)(\S+)(.*)$", line)
        if m:
            out.append(f"{m.group(1)}{m.group(2)}_membrane_avoiding{m.group(3)}\n")
            saw_name_flag = True
            continue

        # 4. insert --ppm_pdb / --ppm_chain_id right before --alphafold_dir
        if (not inserted_ppm_flags) and stripped.startswith("--alphafold_dir"):
            indent = line[: len(line) - len(line.lstrip())]
            out.append(f"{indent}--ppm_pdb ppm_structures/{ppm_pdb} \\\n")
            out.append(f"{indent}--ppm_chain_id {ppm_chain_id} \\\n")
            inserted_ppm_flags = True

        out.append(line)

    if not saw_conda_activate:
        warnings.append("Did not find a line reading exactly 'conda activate proteinhunter' -- env was NOT renamed. Check manually.")
    if not saw_name_flag:
        warnings.append("Did not find a --name flag -- run name was NOT tagged with _membrane_avoiding. Check manually.")
    if not inserted_ppm_flags:
        warnings.append("Did not find an --alphafold_dir flag to anchor on -- --ppm_pdb/--ppm_chain_id were NOT inserted. Add them manually.")

    return out, warnings


def main():
    args = parse_args()
    output_path = args.output or default_output_path(args.input)

    with open(args.input, "r") as fh:
        lines = fh.readlines()

    try:
        new_lines, warnings = convert(lines, args.ppm_pdb, args.ppm_chain_id)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    with open(output_path, "w") as fh:
        fh.writelines(new_lines)
    os.chmod(output_path, 0o755)

    print(f"Wrote {output_path}")
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)


if __name__ == "__main__":
    main()
