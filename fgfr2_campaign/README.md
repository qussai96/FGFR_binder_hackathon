# FGFR2 Binder Design Campaign (Hackathon, In Silico)

This directory contains a local, GPU-friendly workflow for computational binder candidate design against the extracellular domain of human FGFR2 while screening to avoid FGFR1 off-target binding.

## Scientific Focus

- Target: FGFR2 extracellular ligand-binding domain.
- Starting structure: PDB 1DJS (FGF1-FGFR2 complex).
- Off-target specificity structure: PDB 1CVS (FGF2-FGFR1 complex).
- Intent: prioritize binders that block FGF ligand engagement on FGFR2 with increased specificity over FGFR1.
- Scope: in-silico hackathon workflow only.

## Strategies

1. Strategy 1: ligand-site blocker
- Uses FGFR2 residues contacting FGF1 in 1DJS.
- Prioritizes non-conserved interface positions versus FGFR1 from 1CVS.
- Binder length target: 60-120 aa.

2. Strategy 2: activation/dimerization-adjacent surface
- Uses plausible secondary surface near ligand-binding region when activation interface is uncertain from 1DJS alone.
- Down-weights conserved FGFR1-like residues when possible.
- Binder length target: 80-140 aa.

3. Strategy 3: model-free surface explorer
- Broad FGFR2 surface sampling with optional FGFR1-conserved position avoidance.
- Binder length target: 60-140 aa.

## Local Compute Assumptions

- Single NVIDIA A40 GPU (typically 48 GB VRAM).
- Local Linux filesystem is source of truth.
- No Lyceum auth.
- No S3 sync.
- No Modal by default.

## Directory Layout

- structures/: staged local target/off-target structures and inspection output.
- configs/: hotspot JSON, FGFR2-vs-FGFR1 specificity JSON, and per-strategy YAML files.
- scripts/: local shell/python pipeline scripts.
- out/: tool outputs and logs.
- designs/: standardized designs and campaign index.
- rounds/: optional run metadata.

## Round 0 Quick Start

Run:

```bash
bash fgfr2_campaign/scripts/run_round0_local.sh --dry-run
```

Then run without --dry-run once environment/tools are ready.

Round 0 executes:
1. Stage 1DJS.pdb (FGFR2 target) and 1CVS.pdb (FGFR1 off-target) into campaign structures.
2. Inspect chain/entity content.
3. Ask user to confirm FGFR2 receptor and FGF1 ligand chains.
4. Compute FGFR2-FGF1 interface hotspots.
5. Build FGFR2-vs-FGFR1 specificity filter from 1CVS.
6. Generate all strategy configs.
7. Launch only Strategy 1 with 5 designs.
8. Sync local outputs into standardized designs/index.
9. Generate local dashboard pages.

## Safety and Validation Notes

- All outputs are labeled as computational binder candidates.
- Missing metrics are tolerated and reported as warnings.
- Chain-ID guesses and empty interface results are explicitly warned.
- Binder lengths over 250 aa are flagged for review.
- FGFR1 off-target risk is tracked separately and should be used during triage.
- Experimental follow-up should include FGFR2/FGFR1 differential binding assays and signaling inhibition assays.
