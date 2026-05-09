# FGFR2 BindCraft — Handoff guide for collaborators

You're being asked to run a de novo binder design pipeline on Colab. **~5–7 hours of GPU time, fully automated.** Free Colab gives ~12 GPU-hours/day so this fits in one session if you keep the browser tab open.

## What you need

- A Google account with **Colab access** (free tier is fine — needs a T4 GPU)
- A Google Drive with at least **10 GB** free (for AlphaFold params + design outputs)
- That's it. No code to install, no PDBs to download — everything is pulled from public repos automatically.

## Setup (60 seconds)

1. Go to <https://colab.research.google.com>
2. Click **File → Upload notebook** → upload `FGFR2_BindCraft_handoff.ipynb`
3. **Runtime → Change runtime type → Hardware accelerator: T4 GPU → Save**
4. **Runtime → Disconnect and delete runtime** (forces a clean GPU container — important)
5. Click **Connect** in the top-right (or just run a cell)

## Running it

### Cell 1 — Fix dependencies (run once)

Just run the cell. It'll install some Python packages, then **the kernel will deliberately crash with "Your session crashed"** — this is intentional. Wait ~10 seconds, the runtime auto-reconnects.

### Cell 2 — Run BindCraft

Run the cell. First time:
- Mounts Google Drive (will prompt you to authorize)
- Downloads BindCraft + AlphaFold params (~6 minutes, ~5.3 GB)
- Starts designing — outputs ~20 binders in ~3.5 hours

**Important:** keep the browser tab visible. Free Colab idle-disconnects after ~90 minutes of background. If it disconnects, just re-run Cell 2 — BindCraft resumes from saved trajectories in Drive (no work lost).

You'll see progress like:

```
[trajectory 1] FGFR2_binder_l50_s421...
✓ FGFR2_binder_l50_s421_mpnn1  ipSAE_target=0.71  Δ-ipSAE=+0.18  composite=0.55
```

When done you'll see a top-5 table.

### Cell 3 — Live status (run anytime in a separate tab)

Just shows current progress. Useful while Cell 2 is running.

### Cell 4 — Download results

Packages the top 5 PDBs + FASTA + ranking CSV into a zip and downloads to your laptop. Run after Cell 2 finishes.

## What you'll send back

After Cell 2 finishes (or if you have to stop early), Cell 4 produces `fgfr2_top5.zip` (~1 MB). Send that file back. It contains:

- `bindcraft_frontier_ranking.csv` — ranked binders with ipSAE + selectivity scores
- `top_5_for_submission.fasta` — top 5 sequences with metrics in the headers
- 5 PDB structures of the top binders

The full design folder also lives in **your** Google Drive at `/My Drive/BindCraft_FGFR2/` — you can keep it or delete it after.

## Troubleshooting

**"NO GPU"** — runtime isn't on T4. Runtime → Change runtime type → T4 GPU → Save → Disconnect and delete runtime → reconnect.

**"Your session crashed"** in Cell 1 — that's expected. It's the deliberate kernel restart. Just continue to Cell 2.

**"AF2 params: 30 min..." not progressing** — download silently stalled. Stop Cell 2 (■ button), then run this fix cell:

```python
import subprocess, os
subprocess.run(['pkill','-9','aria2c','wget','tar'], check=False)
subprocess.run(['rm','-rf','/content/bindcraft/params'], check=False)
os.makedirs('/content/bindcraft/params', exist_ok=True)
subprocess.check_call(['wget','-q','--show-progress','-O','/content/af2.tar',
    'https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar'])
subprocess.check_call(['tar','-xf','/content/af2.tar','-C','/content/bindcraft/params'])
open('/content/bindcraft/params/done.txt','w').close()
print('✓ AF2 params extracted')
```

Then re-run Cell 2.

**OOM (out of memory)** — T4 has 16 GB; BindCraft uses ~10–12 GB. If a trajectory dies, BindCraft moves on automatically. Not fatal.

**Compute Units exhausted** — free Colab gives ~12 hours/day; you should be fine for 1 run.

## What this pipeline does (TL;DR for context)

It designs small protein binders that block FGF1 from binding FGFR2 (a cancer target) while *not* binding the closely-related FGFR1 (to avoid side effects). Uses AlphaFold2 + ProteinMPNN (frontier ML), ranks by ipSAE (Dunbrack 2025) + paralog selectivity gap. Hotspots and target validated with mCSM-PPI2 + paired PISA + multi-paralog MSA.

If you want the science details, see `fgfr2_campaign/bindcraft/README.md` and `hotspots.json` in the repo: <https://github.com/qussai96/FGFR_binder_hackathon>.

## One thing that helps me

After it runs, **screenshot the final top-5 table** that prints in Cell 2 and send it along with the zip. That gives me a quick read of how the campaign went without unzipping.

Thank you 🙏
