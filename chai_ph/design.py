import os
import argparse
from pipeline import ProteinHunter_Chai


def parse_args():
    """Parses command-line arguments for ProteinHunter_Chai."""
    parser = argparse.ArgumentParser(
        description="ProteinHunter input options for structure prediction and design cycles."
    )

    # --- Job and Sequence Settings ---
    sequence_group = parser.add_argument_group("Sequence and Job Settings")
    sequence_group.add_argument(
        "--jobname", type=str, default="test", help="Job name for output files and folders."
    )
    sequence_group.add_argument(
        "--percent_X",
        type=int,
        default=80,
        help="Percentage of 'X' residues in the initial sequence (0, 50, or 100).",
    )
    sequence_group.add_argument(
        "--seq", type=str, default="", help="Input sequence for the binder chain (optional). If empty, sequence is randomly sampled."
    )
    sequence_group.add_argument(
        "--min_protein_length", type=int, default=100, help="Minimum length of the designed protein chain."
    )
    sequence_group.add_argument(
        "--max_protein_length", type=int, default=150, help="Maximum length of the designed protein chain."
    )
    sequence_group.add_argument(
        "--binder_seq", default="", type=str, help="Initial sequence for the binder chain."
    )
    sequence_group.add_argument(
        "--refiner_mode", action="store_true", default=False, help="Enable refiner mode."
    )
    sequence_group.add_argument(
        "--target_seq",
        type=str,
        default="",
        help="Target sequence (protein) or SMILES (ligand) for binder design (optional).",
    )
    sequence_group.add_argument(
        "--target_pdb",
        type=str,
        default=None,
        help="Optional template PDB file (target_pdb) for design.",
    )
    sequence_group.add_argument(
        "--target_pdb_chain",
        type=str,
        default=None,
        help="Chain ID in template PDB (required if target_pdb provided).",
    )
    sequence_group.add_argument(
        "--cyclic",
        action="store_true",
        default=False,
        help="Enable cyclic topology for the designed chain.",
    )

    # --- Pocket / Contact-Residue Conditioning (Protein Targets Only) ---
    contact_group = parser.add_argument_group("Pocket / Contact-Residue Conditioning")
    contact_group.add_argument(
        "--contact_residues",
        default="",
        type=str,
        help=(
            "Comma-separated 1-indexed target_seq residue numbers the binder must contact. "
            "Uses chai-lab's native pocket restraint (TokenPairPocketRestraint) to actively "
            "steer the diffusion, the same mechanism as Boltz2's pocket constraint."
        ),
    )
    contact_group.add_argument(
        "--ppm_pdb",
        default="",
        type=str,
        help=(
            "Path to a PPM3.0-oriented target structure ('...out.pdb') used to restrict "
            "--contact_residues to membrane-safe residues: solvent-exposed outside the "
            "membrane, or pore/lumen-facing inside a barrel anchored in the membrane. If "
            "--contact_residues is also given, it is intersected with the membrane-safe set "
            "(residues that fall in the lipid-facing membrane region are dropped, with a "
            "warning); if omitted, the full membrane-safe set is used as the pocket-"
            "conditioning contacts. The PPM structure must correspond to the same sequence "
            "and residue numbering as --target_seq."
        ),
    )
    contact_group.add_argument(
        "--ppm_chain_id",
        default="A",
        type=str,
        help="Chain ID to read from --ppm_pdb. Default 'A' matches single-chain PPM3.0 monomer outputs.",
    )
    contact_group.add_argument(
        "--ppm_exclude_pore_facing",
        action="store_true",
        default=False,
        help="Exclude barrel-lumen/pore-facing residues; keep only residues fully outside the membrane.",
    )
    contact_group.add_argument(
        "--pocket_distance_threshold",
        default=15.0,
        type=float,
        help="Distance (Angstroms) used in the pocket-restraint conditioning for each contact residue.",
    )
    contact_group.add_argument(
        "--no_contact_filter",
        action="store_true",
        default=False,
        help="Do not filter or restart for unbound contact residues at cycle 0, and don't gate high-ipTM saves on it.",
    )
    contact_group.add_argument("--max_contact_filter_retries", default=6, type=int)
    contact_group.add_argument(
        "--contact_cutoff",
        default=15.0,
        type=float,
        help="CA-CA distance (Angstroms) used by the post-hoc geometric contact filter/gate.",
    )

    # --- Optimization Options ---
    opt_group = parser.add_argument_group("Optimization and Folding Parameters")
    opt_group.add_argument(
        "--n_trials",
        type=int,
        default=2,
        help="Number of independent optimization trials.",
    )
    opt_group.add_argument(
        "--n_cycles",
        type=int,
        default=8,
        help="Number of folding/design optimization cycles (steps).",
    )
    opt_group.add_argument(
        "--n_recycles",
        type=int,
        default=3,
        help="Number of trunk recycles per fold step.",
    )
    opt_group.add_argument(
        "--n_diff_steps",
        type=int,
        default=200,
        help="Diffusion steps for structure sampling.",
    )
    opt_group.add_argument(
        "--hysteresis_mode",
        type=str,
        default="esm",
        choices=["templates", "esm", "partial_diffusion", "none"],
        help="Strategy for template/feature reuse to guide folding (hysteresis).",
    )
    opt_group.add_argument(
        "--repredict",
        action="store_true",
        default=True,
        help="Re-predict final best structure without templates for validation.",
    )

    # --- MPNN Options ---
    mpnn_group = parser.add_argument_group("ProteinMPNN Settings")
    mpnn_group.add_argument(
        "--omit_aa", type=str, default="", help="Amino acid types to omit from design (e.g., 'C')."
    )

    # temp and bias params
    mpnn_group.add_argument("--temperature", default=0.1, type=float)
    mpnn_group.add_argument("--alanine_bias_start", default=-0.5, type=float)
    mpnn_group.add_argument("--alanine_bias_end", default=-0.1, type=float)
    mpnn_group.add_argument("--alanine_bias", action="store_true")

    mpnn_group.add_argument(
        "--scale_temp_by_plddt",
        action="store_true",
        default=False,
        help="Scale MPNN temperature inversely by pLDDT for focused design.",
    )

    # --- Visualization and Hardware Options ---
    vis_group = parser.add_argument_group("Hardware and Visualization")
    vis_group.add_argument(
        "--show_visual",
        action="store_true",
        default=False,
        help="Show interactive py3Dmol visualization during the run.",
    )
    vis_group.add_argument(
        "--render_freq",
        type=int,
        default=100,
        help="Visualization refresh frequency (diffusion steps).",
    )
    vis_group.add_argument(
        "--plot",
        action="store_true",
        default=False,
        help="Plot cycles figs per run (requires matplotlib).",
    )
    vis_group.add_argument("--gpu_id", type=int, default=0, help="GPU ID to use.")


    af_group = parser.add_argument_group("AlphaFold Settings")
    af_group.add_argument("--use_alphafold3_validation", action="store_true", default=False)
    af_group.add_argument(
        "--alphafold_dir", default=os.path.expanduser("~/alphafold3"), type=str
    )
    af_group.add_argument("--af3_docker_name", default="alphafold3_yc", type=str)
    af_group.add_argument(
        "--af3_database_settings", default="~/alphafold3/alphafold3_data_save", type=str
    )
    af_group.add_argument(
        "--hmmer_path",
        default="~/.conda/envs/alphafold3_venv",
        type=str,
    )
    af_group.add_argument("--use_msa_for_af3", action="store_true")
    af_group.add_argument("--work_dir", default="", type=str)
    af_group.add_argument("--high_iptm_threshold", default=0.8, type=float)
    af_group.add_argument("--high_plddt_threshold", default=0.8, type=float)
    return parser.parse_args()


def main():
    """Main function to run the ProteinHunter pipeline."""
    args = parse_args()
    protein_hunter = ProteinHunter_Chai(args)
    protein_hunter.run_pipeline()


if __name__ == "__main__":
    main()