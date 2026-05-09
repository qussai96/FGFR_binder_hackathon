# FGFR2 Binder Hackathon

Computational protein design workflow combining **BoltzGen** and **RFdiffusion3** to generate FGFR2 binder candidates targeting the FGF1 interface.

## Dashboard

Open the dashboard here:

- [docs/fgfr2/index.html](docs/fgfr2/index.html)

The dashboard shows ranked designs and summary metrics from the latest runs with live filtering and KPI cards.

## Latest Updates

**May 9, 2026 - Pipeline Complete:**
- ✅ **Full design generation pipeline**: BoltzGen (59) + RFdiffusion3 (15) = **74 total designs**
- ✅ **Sequence generation**: Extracted sequences from all outputs
- ✅ **Scoring & evaluation**: Computed composite scores ranked by specificity-adjusted metrics
- ✅ **Dashboard updated**: [docs/fgfr2/index.html](docs/fgfr2/index.html) with interactive filtering
- ✅ **RFdiffusion3 fixed**: Discontinuous residue numbering bug resolved (contig: `A147-161,A163-185,A187-188,A190-217,A219-296,A307-362`)
- ✅ **3 RFdiffusion jobs completed** (jobs 70787-70789): 5 designs per strategy with structures

**Design breakdown:**
- **strategy1_ligand_site_blocker**: 30 (25 BoltzGen + 5 RFdiffusion3)
- **strategy2_dimer_or_activation_blocker**: 22 (17 BoltzGen + 5 RFdiffusion3)
- **strategy3_surface_explorer**: 22 (17 BoltzGen + 5 RFdiffusion3)

## Project Layout

- `fgfr2_campaign/`: Main campaign directory
  - `configs/`: Strategy YAML, target structures, hotspots, specificity filters
  - `configs/rfd3/`: RFdiffusion3 JSON design specs (auto-generated)
  - `scripts/`: Pipeline automation scripts
  - `out/`: Generated designs and logs
- `docs/fgfr2/`: Generated dashboard page and live data index

## Environment Setup

### BoltzGen (CPU/GPU)
```bash
conda activate FGFR_hack  # Python 3.13, torch 2.5.1+cu121
```

### RFdiffusion3 (GPU required)
```bash
conda activate rfd3_py312  # Python 3.12, rfd3 with foundry checkpoints
```

## Quick Start

### Generate BoltzGen Designs
```bash
cd fgfr2_campaign
bash scripts/submit_round0_slurm.sh --all-strategies --relaxed-filters
```

### Generate RFdiffusion Designs
```bash
cd fgfr2_campaign
bash scripts/submit_rfdiffusion_slurm.sh --all-strategies
```

### Update Dashboard
```bash
cd fgfr2_campaign
python3 sync_designs_local.py  # Scan outputs, compute scores
python3 generate_pages_local.py  # Update dashboard HTML/data
# Open docs/fgfr2/index.html
```

## Design Strategies

1. **strategy1_ligand_site_blocker**: Block FGF1 binding interface
2. **strategy2_dimer_or_activation_blocker**: Prevent FGFR2 dimerization/activation
3. **strategy3_surface_explorer**: Broad interface coverage with minimal specificity constraints

## Key Scripts

| Script | Purpose |
|--------|---------|
| `02_generate_strategy_configs.py` | Generate hotspots from interface and specificity filters |
| `04_convert_strategy_to_boltzgen_spec.py` | Convert YAML strategy to BoltzGen native format |
| `06_convert_strategy_to_rfd3_spec.py` | Convert YAML strategy to RFdiffusion3 JSON (handles discontinuous residues) |
| `sync_designs_local.py` | Scan outputs, extract sequences, compute scores |
| `generate_pages_local.py` | Generate interactive dashboard from design index |
| `slurm_boltzgen.sh` | Slurm wrapper for BoltzGen |
| `slurm_rfdiffusion.sh` | Slurm wrapper for RFdiffusion3 |

## Hotspots

User-provided hotspot residues: **A283, A251, A346, A281, A288, A170, A173**

Comprehensive interface hotspots auto-detected from target/ligand structures and stored in:
- `configs/fgfr2_fgf1_interface_hotspots.json`
- `configs/fgfr2_vs_fgfr1_specificity.json`

## Outputs

- **BoltzGen**: JSON designs with sequences, scores, confidence metrics
- **RFdiffusion3**: JSON + CIF structures with diffusion trajectory metadata
- **Dashboard**: Live-updated HTML with filtering, sorting, and KPI summaries

## Troubleshooting

**RFdiffusion validation error (residue not found)**:
- Cause: Contig included missing residue numbers from PDB gaps
- Solution: Rebuild specs with `06_convert_strategy_to_rfd3_spec.py` (now auto-detects and skips missing residues)

**Dashboard shows no designs**:
- Run `sync_designs_local.py` to regenerate `docs/fgfr2/data/index.json`
- Then run `generate_pages_local.py` to update HTML

**GPU/CUDA issues**:
- Ensure correct PyTorch + CUDA version: `torch 2.5.1+cu121` for BoltzGen
- Use local rfd3 (`rfd3_py312`) for RFdiffusion3 instead of apptainer if available
