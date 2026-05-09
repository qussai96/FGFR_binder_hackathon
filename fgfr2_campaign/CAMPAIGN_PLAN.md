# FGFR2 Campaign Plan (Local A40, In Silico)

## Objective

Generate and rank computational protein binder candidates against the extracellular domain of human FGFR2, prioritizing candidates that may block FGF ligand binding while reducing FGFR1 off-target binding risk.

## Target and Off-target Definition

- Primary target structure: PDB 1DJS (FGF1-FGFR2 complex).
- Target input file: fgfr2_campaign/structures/1DJS.pdb.
- Off-target specificity structure: PDB 1CVS (FGF2-FGFR1 complex).
- Off-target input file: fgfr2_campaign/structures/1CVS_offtarget_fgfr1.pdb.
- Target receptor: FGFR2 extracellular domain chain (confirmed from inspection).
- Target interface guidance: FGFR2 residues contacting FGF1 at <= 5.0 A heavy-atom distance.
- Specificity guidance: prioritize FGFR2 interface residues that are non-conserved versus mapped FGFR1 receptor positions.

## Local Compute Constraints

- Environment: local Linux cluster.
- GPU: one NVIDIA A40.
- Storage: local filesystem only.
- No Lyceum, no S3, no Modal by default.
- Initial batch size: 5-10 designs per strategy.

## Round-0 Workflow

1. Stage and inspect target/off-target structures.
2. Confirm FGFR2 receptor and FGF1 ligand chains manually.
3. Compute FGFR2-FGF1 interface hotspots.
4. Build FGFR2-vs-FGFR1 specificity filter from 1CVS.
5. Generate strategy configs from hotspots plus specificity JSON.
6. Run strategy1 only (5 designs).
7. Sync local outputs into standardized design index.
8. Generate local dashboard.

## Design Strategies

1. strategy1_ligand_site_blocker
- Purpose: bind FGFR2 surface used by FGF1.
- Specificity: prefer non-conserved FGFR2 hotspots, avoid conserved FGFR1-like hotspots.
- Binder length: 60-120 aa.

2. strategy2_dimer_or_activation_blocker
- Purpose: bind plausible activation/dimerization-adjacent surface.
- Specificity: carry forward FGFR1-conserved hotspot avoidance where available.
- Binder length: 80-140 aa.

3. strategy3_surface_explorer
- Purpose: broad model-free FGFR2 surface binding exploration.
- Specificity: optional anti-FGFR1 filtering in ranking/evaluation.
- Binder length: 60-140 aa.

## Ranking and Selection

Use sync_designs_local.py to aggregate metrics and compute composite score:

- 0.30 * normalized ipTM
- 0.20 * normalized pLDDT
- 0.20 * normalized inverse min_interaction_pae
- 0.15 * normalized inverse rmsd
- 0.10 * no_clash_bonus
- 0.05 * diversity_bonus

Missing metrics are allowed and reported; partial scores are computed from available components.

Specificity extension:
- Add FGFR1 off-target metrics when available (for example, off_target_ipTM or specificity_gap).
- Penalize candidates with strong FGFR1 predicted binding.

## Evaluation Plan

- Select top 10 unevaluated candidates maximum per pass.
- One GPU job at a time.
- Prepare FGFR2+binder and FGFR1+binder inputs for local Boltz-2 or AF-style validation.
- Merge evaluation outputs back into design index via sync_designs_local.py.

## Safety and Scope

- All outputs are computational binder candidates only.
- No efficacy or inhibitor claims.
- Experimental validation is required (target and off-target binding plus signaling assays).
