"""
Membrane-aware residue classification from PPM3.0-oriented structures.

PPM3.0 (https://opm.phar.umich.edu/ppm_server) reorients a protein so the
membrane normal is the z-axis and the membrane midplane sits at z=0, and
writes:
  - a REMARK line: "REMARK      1/2 of bilayer thickness:   <float>"
  - two rings of HETATM "DUM" pseudo-atoms (atom names "N"/"O") marking the
    +/- hydrophobic-core boundary planes.

This module uses that output to classify every residue of a target chain
into one of three categories:

  - "exposed":      CA lies outside the membrane z-band. Always a valid
                     binder-design target.
  - "pore_facing":  CA lies inside the membrane z-band, but the side chain
                     (CA->CB direction) points toward the local channel/
                     bundle axis rather than away from it. Still a valid
                     target -- e.g. the aqueous lumen of a beta barrel
                     (BamA, LptD, OmpA, a porin, ...).
  - "lipid_facing": CA lies inside the membrane z-band and the side chain
                     points away from the axis, into the lipid bilayer.
                     NOT reachable by a soluble binder -- excluded.

Geometry note: a barrel (or helix-bundle) wall is a single residue layer
thick, so the backbone (CA) radial distance from the pore/bundle axis is
nearly identical for lipid-facing and pore-facing residues in the same
strand/helix -- it is the side-chain *direction* that alternates in and out,
not the backbone radius. Classification is therefore done via the CA->CB
radial-direction test, not a radial-distance threshold.
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

REMARK_THICKNESS_RE = re.compile(
    r"1/2\s+of\s+bilayer\s+thickness\s*:\s*([\-0-9.]+)", re.IGNORECASE
)

# Standard fixed-width PDB ATOM/HETATM column layout (1-indexed spec,
# expressed here as 0-indexed python slices).
_REC = slice(0, 6)
_NAME = slice(12, 16)
_ALT = 16
_RESNAME = slice(17, 20)
_CHAIN = 21
_RESSEQ = slice(22, 26)
_ICODE = 26
_X = slice(30, 38)
_Y = slice(38, 46)
_Z = slice(46, 54)


@dataclass
class ResidueAtoms:
    resname: str
    N: Optional[np.ndarray] = None
    CA: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    CB: Optional[np.ndarray] = None


def _parse_pdb_atoms(pdb_path: str) -> Tuple[Dict[str, Dict[int, ResidueAtoms]], List[np.ndarray], Optional[float]]:
    """
    Parses ATOM records (grouped by chain / residue seq number) and DUM
    HETATM records, plus the PPM bilayer-thickness REMARK.

    Returns (chains, dum_coords, half_thickness_from_remark).
    """
    chains: Dict[str, Dict[int, ResidueAtoms]] = {}
    dum_coords: List[np.ndarray] = []
    half_thickness: Optional[float] = None

    with open(pdb_path, "r") as fh:
        for line in fh:
            if line.startswith("REMARK") and half_thickness is None:
                m = REMARK_THICKNESS_RE.search(line)
                if m:
                    half_thickness = float(m.group(1))
                continue

            is_atom = line.startswith("ATOM")
            is_hetatm = line.startswith("HETATM")
            if not (is_atom or is_hetatm):
                continue
            if len(line) < 54:
                continue

            atom_name = line[_NAME].strip()
            resname = line[_RESNAME].strip()
            try:
                x = float(line[_X])
                y = float(line[_Y])
                z = float(line[_Z])
            except ValueError:
                continue
            coord = np.array([x, y, z], dtype=float)

            if is_hetatm and resname == "DUM":
                dum_coords.append(coord)
                continue

            if not is_atom:
                continue

            chain_id = line[_CHAIN].strip() or " "
            try:
                resseq = int(line[_RESSEQ])
            except ValueError:
                continue

            res_dict = chains.setdefault(chain_id, {})
            res = res_dict.setdefault(resseq, ResidueAtoms(resname=resname))

            if atom_name == "N" and res.N is None:
                res.N = coord
            elif atom_name == "CA" and res.CA is None:
                res.CA = coord
            elif atom_name == "C" and res.C is None:
                res.C = coord
            elif atom_name == "CB" and res.CB is None:
                res.CB = coord

    return chains, dum_coords, half_thickness


def reconstruct_virtual_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Places an idealized virtual CB from backbone N/CA/C coordinates.
    Used for glycine (no CB) and for any residue missing an explicit CB.
    """
    b = ca - n
    c_vec = c - ca
    a = np.cross(b, c_vec)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * c_vec + ca


def get_half_bilayer_thickness(pdb_path: str) -> float:
    """
    Returns half the bilayer (hydrophobic core) thickness in Angstroms,
    read from the PPM3.0 REMARK line, falling back to the mean |z| of the
    DUM boundary pseudo-atoms if the REMARK is missing.
    """
    _, dum_coords, half_thickness = _parse_pdb_atoms(pdb_path)
    if half_thickness is not None:
        return half_thickness
    if dum_coords:
        return float(np.mean([abs(c[2]) for c in dum_coords]))
    raise ValueError(
        f"Could not determine membrane bilayer thickness from {pdb_path}: "
        "no PPM REMARK line and no DUM pseudo-atoms found. "
        "Is this a PPM3.0-oriented '...out.pdb' file?"
    )


def _detect_chain_id(chains: Dict[str, Dict[int, ResidueAtoms]], requested: Optional[str]) -> str:
    if requested is not None:
        if requested not in chains:
            raise ValueError(
                f"Chain '{requested}' not found in structure "
                f"(available chains: {sorted(chains.keys())})"
            )
        return requested
    if len(chains) == 1:
        return next(iter(chains))
    raise ValueError(
        f"Structure has multiple chains ({sorted(chains.keys())}); "
        "specify chain_id explicitly."
    )


def classify_residues(
    pdb_path: str,
    chain_id: Optional[str] = None,
) -> Dict[int, dict]:
    """
    Classifies every residue of one chain in a PPM3.0-oriented PDB into
    "exposed", "pore_facing", or "lipid_facing".

    Returns {resnum: {"category": str, "resname": str, "z": float}}.
    """
    chains, _, half_thickness = _parse_pdb_atoms(pdb_path)
    if half_thickness is None:
        half_thickness = get_half_bilayer_thickness(pdb_path)

    resolved_chain = _detect_chain_id(chains, chain_id)
    residues = chains[resolved_chain]

    # Local channel/bundle axis: xy-centroid of CA atoms whose z falls
    # inside the membrane band.
    in_band_ca_xy = [
        res.CA[:2]
        for res in residues.values()
        if res.CA is not None and abs(res.CA[2]) <= half_thickness
    ]
    axis_xy = (
        np.mean(np.stack(in_band_ca_xy), axis=0)
        if in_band_ca_xy
        else np.array([0.0, 0.0])
    )

    classification: Dict[int, dict] = {}
    for resnum, res in sorted(residues.items()):
        if res.CA is None:
            continue
        z = float(res.CA[2])

        if abs(z) > half_thickness:
            classification[resnum] = {
                "category": "exposed",
                "resname": res.resname,
                "z": z,
            }
            continue

        cb = res.CB
        if cb is None:
            if res.N is not None and res.C is not None:
                cb = reconstruct_virtual_cb(res.N, res.CA, res.C)
            else:
                # Can't determine side-chain direction; exclude conservatively.
                classification[resnum] = {
                    "category": "lipid_facing",
                    "resname": res.resname,
                    "z": z,
                }
                continue

        r_ca = float(np.linalg.norm(res.CA[:2] - axis_xy))
        r_cb = float(np.linalg.norm(cb[:2] - axis_xy))

        category = "lipid_facing" if r_cb > r_ca else "pore_facing"
        classification[resnum] = {
            "category": category,
            "resname": res.resname,
            "z": z,
        }

    return classification


def get_membrane_safe_residues(
    pdb_path: str,
    chain_id: Optional[str] = None,
    include_pore_facing: bool = True,
) -> List[int]:
    """
    Returns the sorted list of 1-indexed residue numbers that are valid
    binder-design targets: solvent-exposed outside the membrane, and
    (unless include_pore_facing=False) pore/lumen-facing residues inside
    barrels anchored in the membrane.
    """
    classification = classify_residues(pdb_path, chain_id=chain_id)
    allowed = {"exposed"} | ({"pore_facing"} if include_pore_facing else set())
    return sorted(
        resnum
        for resnum, info in classification.items()
        if info["category"] in allowed
    )


def summarize_classification(classification: Dict[int, dict]) -> str:
    counts = {"exposed": 0, "pore_facing": 0, "lipid_facing": 0}
    for info in classification.values():
        counts[info["category"]] += 1
    total = sum(counts.values())
    return (
        f"{total} residues classified: "
        f"{counts['exposed']} exposed, "
        f"{counts['pore_facing']} pore-facing, "
        f"{counts['lipid_facing']} lipid-facing (excluded)"
    )
