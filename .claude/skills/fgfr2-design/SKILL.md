# FGFR2 Design Skill

Use this skill for FGFR2 computational binder generation workflows.

## Mandatory Context

- This campaign targets FGFR2, not PGDH.
- Local NVIDIA A40 execution only.
- No Lyceum, no S3, and no Modal by default.
- All generation outputs must go under fgfr2_campaign/out/.

## Required Steps

1. Confirm target structure and chain IDs from 1DJS inspection outputs.
2. Use fgfr2_campaign/scripts/01_find_interface_residues.py for strategy1 hotspot derivation.
3. Use local runner scripts in fgfr2_campaign/scripts/.
4. Keep initial batches small (5-10 designs).
5. Run fgfr2_campaign/sync_designs_local.py after each generation step.

## Failure Handling

- Do not delete failed designs.
- Mark failed runs clearly with warnings or status fields.
- Preserve all artifacts for reproducibility.
