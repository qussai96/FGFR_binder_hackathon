# FGFR2 Evaluate Skill

Use this skill for local FGFR2 candidate evaluation and reranking.

## Mandatory Context

- This campaign targets FGFR2, not PGDH.
- Local NVIDIA A40 execution only.
- No Lyceum, no S3, and no Modal by default.
- All outputs must be stored under fgfr2_campaign/out/.

## Required Steps

1. Read fgfr2_campaign/designs/index.json.
2. Select top unevaluated candidates (<=10 for first pass).
3. Prepare local receptor+binder evaluation inputs.
4. Run local evaluation only if installed.
5. If evaluation tooling is missing, emit TODO commands and status warnings.
6. Run fgfr2_campaign/sync_designs_local.py after each evaluation step.

## Failure Handling

- Do not delete failed designs.
- Mark failures explicitly in evaluation metrics and warnings.
- Keep unsuccessful outputs for traceability.
