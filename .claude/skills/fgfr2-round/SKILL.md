# FGFR2 Round Skill

Use this skill to execute and track iterative FGFR2 design rounds.

## Mandatory Context

- This campaign targets FGFR2, not PGDH.
- Local NVIDIA A40 execution only.
- No Lyceum, no S3, and no Modal by default.
- Outputs must go to fgfr2_campaign/out/.

## Round Workflow

1. Confirm target chains and hotspot derivation inputs.
2. Launch one strategy at a time for controlled GPU utilization.
3. Use dry-run mode before real runs when changing configs.
4. Sync designs after generation and after evaluation.
5. Generate/refresh local dashboard pages after sync.

## Failure Handling

- Never delete failed designs.
- Mark failed or partial outputs as failed with clear notes.
- Preserve artifacts and logs for debugging and reproducibility.
