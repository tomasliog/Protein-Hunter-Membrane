#!/bin/bash
#
# Template for a membrane-avoiding Protein Hunter SLURM job.
# Copy this file per run (matching the existing run_ph_<Target>_<RunTag>.sh
# convention), fill in every __FILL_ME_...__ placeholder below, then `sbatch`
# it from the Protein-Hunter-Membrane repo root.
#
# NOTE: placeholders deliberately avoid "<<" and bare "<"/">" -- those are
# heredoc/redirection operators in bash. An unedited "<<FOO>>" placeholder
# would silently turn the rest of the script into a heredoc body instead of
# failing loudly, so don't reintroduce that style here.

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH -t 6:00:00
#SBATCH -J protein_hunter_membrane___FILL_ME_TARGET____FILL_ME_RUN_TAG__
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=tomaslio@mit.edu

set -euo pipefail

START_EPOCH=$(date +%s)
START_HUMAN=$(date)
echo "Job start time: ${START_HUMAN}"

finish() {
  local end_epoch end_human elapsed h m s
  end_epoch=$(date +%s)
  end_human=$(date)
  elapsed=$((end_epoch - START_EPOCH))
  h=$((elapsed / 3600))
  m=$(((elapsed % 3600) / 60))
  s=$((elapsed % 60))
  echo "Job end time: ${end_human}"
  echo "Total runtime: ${h}h ${m}m ${s}s"
}
trap finish EXIT

echo "===== Job info ====="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node(s): ${SLURM_JOB_NODELIST}"
echo "CPUs: ${SLURM_NTASKS}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES-}"
date
echo "===================="

# Make a logs directory for stdout/stderr if it doesn't exist
mkdir -p logs

# (Optional) ensure you start in the directory you submitted from
cd "${SLURM_SUBMIT_DIR}"

# ---- Conda activation (works whether or not conda is already initialized) ----
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: conda.sh not found at $HOME/miniconda3/etc/profile.d/conda.sh"
  exit 1
fi

conda activate proteinhunter_membrane

echo "===== Environment ====="
which python
python --version
python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
print("device_count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("device 0:", torch.cuda.get_device_name(0))
PY
echo "======================="

# Helpful GPU visibility check
nvidia-smi || true

# ---- Run ProteinHunter (membrane-avoiding) ----
# --ppm_pdb restricts (or, with no --contact_residues, fully populates) the
# pocket-conditioning contact set to membrane-safe residues: solvent-exposed
# outside the membrane, or pore/lumen-facing inside a barrel. The PPM
# structure must correspond to the exact sequence/numbering in --protein_seqs.
python boltz_ph/design.py \
  --num_designs 200 \
  --num_cycles 6 \
  --protein_seqs __FILL_ME_TARGET_SEQUENCE__ \
  --msa_mode mmseqs \
  --gpu_id 0 \
  --name __FILL_ME_TARGET____FILL_ME_RUN_TAG___membrane_avoiding \
  --percent_X 70 \
  --min_protein_length 10 \
  --max_protein_length 90 \
  --high_iptm_threshold 0.7 \
  --use_msa_for_af3 \
  --plot \
  --ppm_pdb ppm_structures/__FILL_ME_TARGET___out.pdb \
  --ppm_chain_id A \
  --alphafold_dir /home/tomaslio/models/af3 \
  --af3_docker_name /orcd/software/community/001/container_images/alphafold3/20250321/alphafold3.sif \
  --hmmer_path "$HOME/miniconda3/envs/alphafold3_runner" \
  --af3_database_settings /orcd/datasets/001/alphafold3

echo "Done."
date
