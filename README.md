# Protein-Hunter 😈 (Membrane-Aware Fork)

> This fork adds **PPM3.0-based membrane avoidance** for designing binders against outer-membrane/other membrane proteins: binders are only conditioned to contact residues that are solvent-exposed outside the membrane, or that line the aqueous lumen of a barrel anchored in the membrane — never residues embedded in the lipid bilayer. See [🧬 Membrane-Avoiding Binder Design](#-membrane-avoiding-binder-design-ppm30) below.

<p align="center">
  <img src="./protein_hunter.png" alt="Protein Hunter" width="500"/>
</p>

<p align="center" style="font-size:90%">
  <em>DAlphaBall.gcc
    <strong>Note:</strong> Logo is a ChatGPT-modified version of the Netflix animation and is for illustration only.
  </em>
</p>

> 📄 **Paper**: [Protein Hunter: exploiting structure hallucination
within diffusion for protein design](https://www.biorxiv.org/content/10.1101/2025.10.10.681530v2.full.pdf)
---


## 📝 Colab Notebook

[![Open Protein Hunter (Chai)](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yehlincho/Protein-Hunter/blob/main/protein_hunter_chai_colab.ipynb)
[![Open Protein Hunter (Boltz)](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yehlincho/Protein-Hunter/blob/main/protein_hunter_boltz_colab.ipynb)


## 🚀 Quick Start

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yehlincho/Protein-Hunter.git
   cd Protein-Hunter
   ```

2. **Run the automated setup**
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```

> ⚠️ **Note**: AlphaFold3 setup is not included. Please install it separately by following the [official instructions](https://github.com/google-deepmind/alphafold3).

The setup script will automatically:
- ✅ Create a conda environment with Python 3.10
- ✅ Install all required dependencies
- ✅ Set up a Jupyter kernel for notebooks
- ✅ Download Boltz and Chai weights
- ✅ Configure LigandMPNN and ProteinMPNN
- ✅ Optionally install PyRosetta
- ❌ AF3 must be installed separately

---

We have implemented two different AF3-style models in our Protein Hunter pipeline (more models will be added in the future):
- Boltz1/2
- Chai1


## 🧬 Membrane-Avoiding Binder Design (PPM3.0)

For membrane proteins (e.g. outer-membrane barrels like BamA, LptD, OmpA, porins), you generally don't want a soluble binder conditioned to contact residues that are actually embedded in the lipid bilayer — they're not reachable. This fork uses [PPM3.0](https://opm.phar.umich.edu/ppm_server) output (`utils/membrane_utils.py`) to classify every target residue into:

- **exposed** — Cα outside the membrane z-band (from the PPM `1/2 of bilayer thickness` boundary): always a valid target.
- **pore_facing** — Cα inside the membrane band, but the side chain (Cα→Cβ direction) points *toward* the local channel/barrel axis rather than away from it: still valid (e.g. the aqueous lumen of a beta barrel).
- **lipid_facing** — Cα inside the membrane band with the side chain pointing *away* from the axis, into the lipid: excluded.

The exposed/pore-facing distinction is a direction test, not a distance threshold: a barrel wall is one residue thick, so backbone Cα radius from the pore axis is nearly identical for lipid- and pore-facing residues in the same strand — only the side-chain direction alternates.

**1. Get a PPM3.0 structure.** Run your target (e.g. an AF3/Boltz-folded monomer) through the [PPM3.0 server](https://opm.phar.umich.edu/ppm_server) — this reorients the structure so the membrane normal is the z-axis and writes the `...out.pdb` used below.

**2. Standalone: derive a residue list.** To inspect the classification or paste a residue list directly into `--contact_residues`:
```bash
python utils/ppm_contact_residues.py --ppm_pdb /path/to/target_out.pdb
# prints a summary to stderr and a comma-separated residue list to stdout
```

**3. Pipeline integration.** Both backends accept `--ppm_pdb` directly — it restricts (or, if `--contact_residues` is omitted, fully populates) the pocket-conditioning contact set to membrane-safe residues:

```bash
# Boltz edition
python boltz_ph/design.py \
    --protein_seqs "<target sequence>" \
    --ppm_pdb /path/to/target_out.pdb \
    --name bama_membrane_avoiding \
    --num_designs 20 --num_cycles 7

# Chai edition
python chai_ph/design.py \
    --target_seq "<target sequence>" \
    --ppm_pdb /path/to/target_out.pdb \
    --jobname bama_membrane_avoiding
```

If `--contact_residues` is also given, it's intersected with the membrane-safe set (any residues that fall in the lipid-facing region are dropped, with a warning). Add `--ppm_exclude_pore_facing` to keep only fully-exposed residues (drop barrel-lumen residues too).

Both backends condition generation the same way Boltz2's native `pocket` constraint works — actively biasing the diffusion trajectory, not just filtering results afterward. `chai_ph` reaches this via chai-lab's own (previously unused in this repo) `TokenPairPocketRestraint` feature; both backends also run a post-cycle-0 geometric retry and a high-ipTM save gate against the same residue set (`--contact_cutoff`, `--no_contact_filter`, `--max_contact_filter_retries`).


## Run Code End-to-End

## 1️⃣ End-to-end structure and sequence generation  
👉 See example usage in `run_protein_hunter.py` for reference. 🐍✨  
This will take you from your initial input all the way to final designed protein structures and sequences, all in one pipeline!
The design chain is "A" and other target chains are "B", "C", etc.



> 💡 **Tips:** The original evaluation in the paper used an all-X sequence for initial design. However, to increase the diversity of generated folds, you can mix random amino acids with X residues by setting the `percent_X` parameter (e.g., `--percent_X 50` for 50% X and 50% random AAs). Adjusting this ratio helps explore a broader design space—but if you see too many floating or disconnected structures, try decreasing `percent_X` to encourage more structured designs.


⚠️ **Warning**: To run the AlphaFold3 cross-validation pipeline, you need to specify your AlphaFold3 directory, Docker name, database settings, and conda environment in the configuration. These can be set using the following arguments:
- `--alphafold_dir`: Path to your AlphaFold3 installation (default: ~/alphafold3)
- `--af3_docker_name`: Name of your AlphaFold3 Docker container
- `--af3_database_settings`: Path to AlphaFold3 database
- `--af3_hmmer_path`: Path to HMMER
- `--use_alphafold3_validation`: Add this flag to enable AlphaFold3-based validation. 

## Protein Hunter (Boltz Edition ⚡) 
To use AlphaFold3 validation, make sure your AlphaFold3 Docker is installed, specify the correct AlphaFold3 directory, and turn on `--use_alphafold3_validation`.


- **Protein-protein design with all X sequence:**  
  To design a protein-protein complex using an all-X sequence (i.e., X for every residue, encouraging de novo exploration), run:  
  ```
  python boltz_ph/design.py --num_designs 3 --num_cycles 7 --protein_seqs AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --msa_mode "mmseqs" --gpu_id 0 --name PDL1_mix_aa_all_X --percent_X 100 --min_protein_length 90 --max_protein_length 150 --high_iptm_threshold 0.7 --use_msa_for_af3 --plot
  ```
  > 💡 **Tip:** The `--percent_X 100` flag ensures all positions use the X (unknown) amino acid code.

- **Protein-protein design (mixed X and random amino acids):**  
  To design a protein-protein complex using a mix of X and random amino acids for the designable chain, run:  
  ```
  python boltz_ph/design.py --num_designs 3 --num_cycles 7 --protein_seqs AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --msa_mode "mmseqs" --gpu_id 0 --name PDL1_mix_aa --percent_X 50 --min_protein_length 90 --max_protein_length 150 --high_iptm_threshold 0.7 --use_msa_for_af3 --plot
  ```

  - **Protein-protein design (sample fewer alanine with alanine bias):**  
  To design a protein-protein complex using a mix of X and random amino acids while discouraging alanine sampling during design, run:
  ```
  python boltz_ph/design.py --num_designs 3 --num_cycles 7 --protein_seqs AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --msa_mode "mmseqs" --gpu_id 0 --name PDL1_mix_aa_alanine_bias --percent_X 50 --min_protein_length 90 --max_protein_length 150 --high_iptm_threshold 0.7 --alanine_bias --use_msa_for_af3 --plot
  ```

- **Protein-protein contact specification design:**  
  By default, contact potentials are disabled. To enable contact-based potentials and specify interface residue positions (e.g., residue positions "2,3,10" in the target chain), add `--no_potentials False` and `--contact_residues 2,3,10` to your command. For example:
  ```
  python boltz_ph/design.py --num_designs 3 --num_cycles 7 --protein_seqs AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --msa_mode "mmseqs" --gpu_id 0 --name PDL1_mix_aa --min_protein_length 90 --max_protein_length 150 --high_iptm_threshold 0.7 --use_msa_for_af3 --plot --no_potentials False --contact_residues 2,3,10
  ```
  
- **Multimer binder design:**  
  To design a binder for a multimeric protein (e.g., a dimer), separate chain sequences using `:`
  ```
  python boltz_ph/design.py --num_designs 3 --num_cycles 7 --protein_seqs AGIKVFGHPASIATRRVLIALHEKNLDFELVHVELKDGEHKKEPFLSRNPFGQVPAFEDGDLKLFESRAITQYIAHRYENQGTNLLQTDSKNISQYAIMAIGMQVEDHQFDPVASKLAFEQIFKSIYGLTTDEAVVAEEEAKLAKVLDVYEARLKEFKYLAGETFTLTDLHHIPAIQYLLGTPTKKLFTERPRVNEWVAEITKRPASEKVQ:AGIKVFGHPASIATRRVLIALHEKNLDFELVHVELKDGEHKKEPFLSRNPFGQVPAFEDGDLKLFESRAITQYIAHRYENQGTNLLQTDSKNISQYAIMAIGMQVEDHQFDPVASKLAFEQIFKSIYGLTTDEAVVAEEEAKLAKVLDVYEARLKEFKYLAGETFTLTDLHHIPAIQYLLGTPTKKLFTERPRVNEWVAEITKRPASEKVQ --msa_mode "mmseqs" --gpu_id 0 --name 1GNW_mix_aa --min_protein_length 90 --max_protein_length 150 --high_iptm_threshold 0.7 --use_msa_for_af3 --plot
  ```


- **Cyclic peptide binder design**
  ```
  python boltz_ph/design.py --num_designs 3 --num_cycles 7 --protein_seqs AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --msa_mode "mmseqs" --gpu_id 0 --name PDL1_cyclic_peptide_binder --min_protein_length 10 --max_protein_length 20 --high_iptm_threshold 0.8 --percent_X 100 --use_msa_for_af3 --plot --cyclic
  ```

- **Small molecule binder design:**  
  For designing a protein binder for a small molecule (e.g., SAM), use:  
  ```
  python boltz_ph/design.py --num_designs 5 --num_cycles 7 --ligand_ccd SAM --gpu_id 2 --name SAM_binder --min_protein_length 130 --max_protein_length 150 --high_iptm_threshold 0.7 --use_msa_for_af3 --plot
  ```

- **DNA/RNA PDB design:**  
  To design a protein binder for a nucleic acid (e.g., an RNA sequence), run:
  ```
  python boltz_ph/design.py --num_designs 5 --num_cycles 7 --nucleic_seq AGAGAGAGA --nucleic_type rna --gpu_id 0 --name RNA_bind --min_protein_length 130 --max_protein_length 150 --high_iptm_threshold 0.7 --use_msa_for_af3 --plot
  ```

- **Designs with multiple/heterogeneous target types:**  
  Want to target multiple types of molecules (e.g., a protein with a ligand and a template)? Run:
  ```
  python boltz_ph/design.py --num_designs 5 --num_cycles 7 --protein_seqs AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --msa_mode "mmseqs" --ligand_ccd SAM --gpu_id 0 --name PDL1_SAM --min_protein_length 90 --max_protein_length 150 --high_iptm_threshold 0.8 --use_msa_for_af3 --plot
  ```

## Protein Hunter (Chai Edition ☕) 

> ⚠️ **Caution:** The Chai version is under active development
- [x] Support for multiple targets
- [x] Full AlphaFold (AF3) validation

- **Unconditional protein design:**  
  Generate de novo proteins of a desired length:
  ```
  python chai_ph/design.py --jobname unconditional_design --min_protein_length 150 --max_protein_length 150 --percent_X 0 --seq "" --target_seq ACDEFGHIKLMNPQRSTVWY --n_trials 1 --n_cycles 5 --n_recycles 3 --n_diff_steps 200 --hysteresis_mode templates --repredict --omit_aa "" --temperature 0.1 --scale_temp_by_plddt --render_freq 100 --gpu_id 0 --plot

- **Protein binder design:**  
  To design a binder for a specific protein target (e.g., PDL1):
  ```
  python chai_ph/design.py --jobname PDL1_binder --min_protein_length 90 --max_protein_length 150 --percent_X 80 --seq "" --target_seq AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --n_trials 1 --n_cycles 5 --n_recycles 3 --n_diff_steps 200 --hysteresis_mode templates --repredict --omit_aa "" --temperature 0.1 --scale_temp_by_plddt --render_freq 100 --gpu_id 0 --use_msa_for_af3 --plot
  ```

- **Cyclic peptide binder design:**  
  Design a cyclic peptide binder for a specific protein target:
  ```
  python chai_ph/design.py --jobname PDL1_cyclic_binder --percent_X 80 --seq "" --min_protein_length 10 --max_protein_length 20 --cyclic --target_seq AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYAAALE --n_trials 1 --n_cycles 5 --n_recycles 3 --n_diff_steps 200 --hysteresis_mode templates --repredict --omit_aa "" --temperature 0.1 --scale_temp_by_plddt --render_freq 100 --gpu_id 0 --use_msa_for_af3 --plot
  ```

- **Small molecule (ligand) binder design:**  
  To design a protein binder for a small molecule or ligand (SMILES string as target):
  ```
  python chai_ph/design.py --jobname ligand_binder  --percent_X 0 --seq "" --min_protein_length 130 --max_protein_length 150 --target_seq O=C(NCc1cocn1)c1cnn(C)c1C(=O)Nc1ccn2cc(nc2n1)c1ccccc1 --n_trials 1 --n_cycles 5 --n_recycles 3 --n_diff_steps 200 --hysteresis_mode esm --repredict --omit_aa "" --temperature 0.01 --scale_temp_by_plddt --render_freq 100 --gpu_id 2 --plot
  ```

## 2️⃣ Refine your own designs!
🛠️ You can provide your initial designs as input and further improve their structures by iteratively redesigning and predicting them. Repeat as needed for optimal results!

See the code in `refiner.ipynb` for example usage.

For example, you can generate a design using Boltzgen, take the final output, and refine it further using the iterative pipeline. 

---




## 🎥 Trajectory Visualization
We have implemented trajectory visualization using LogMD and py2Dmol (developed by Sergey Ovchinnikov).

## ✅ Structure Validation

### Primary Evaluation: AlphaFold3
Final structures are validated using **AlphaFold3** for:
- Structure quality assessment 
- Confidence scoring
- Cross-validation against design targets

## 🎯 Successful Designs

By default, the `high_iptm_yaml` folder contains only the designs with ipTM above your chosen threshold and alanine percentage below 20%. If you want to browse **all** generated designs, please check the metrics in `summary_all_runs.csv`.

After running the pipeline with `run_protein_hunter.py`, high-confidence designs can be found in:

- `your_output_folder/high_iptm_yaml`
- `your_output_folder/high_iptm_cif`
- Metrics summary: `your_output_folder/summary_high_iptm.csv`

After running AlphaFold3, the validated (successfully predicted) structures are saved in:

`your_output_folder/03_af_pdb_success`


## 📝 To-Do List

- [ ] Add cross-validation between Boltz and Chai (both directions) without using AlphaFold3 as an option
- [ ] Implement multi-timer support for Protein Hunter Chai edition
- [ ] Specify multiple contacts on multiple targets 
- [ ] Explore other cool applications

Collaboration is always welcome! Email me. Let's chat.

## 📄 License & Citation

**License**: MIT License - See LICENSE file for details  
**Citation**: If you use Protein Hunter in your research, please cite:
```
@article{cho2025protein,
  title={Protein Hunter: exploiting structure hallucination within diffusion for protein design},
  author={Cho, Yehlin and Rangel, Griffin and Bhardwaj, Gaurav and Ovchinnikov, Sergey},
  journal={bioRxiv},
  pages={2025--10},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
```
---

## 📧 Contact & Support

**Questions or Collaboration**: yehlin@mit.edu  
**Issues**: Please report bugs and feature requests via GitHub Issues

---

## ⚠️ Important Disclaimer

> **EXPERIMENTAL SOFTWARE**: This pipeline is under active development and has **NOT been experimentally validated** in laboratory settings. We release this code to enable community contributions and collaborative development. Use at your own discretion and validate results independently.

