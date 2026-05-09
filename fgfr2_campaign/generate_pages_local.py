#!/usr/bin/env python3
"""Generate local FGFR2 dashboard pages from design index.

This version is fully offline for structure display: it creates local PNG previews
from structure files and embeds those previews directly in each design card.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Bio.PDB import MMCIFParser, PDBParser


def _safe_name(name: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in name)


def _copy_structure(src_path: str, dst_dir: Path, design_id: str) -> Path | None:
    if not src_path or not os.path.isfile(src_path):
        return None

    safe_id = _safe_name(design_id)
    try:
        if src_path.endswith(".cif.gz"):
            dst_path = dst_dir / f"{safe_id}.cif"
            with gzip.open(src_path, "rb") as f_in, open(dst_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            return dst_path

        # Some outputs use .cif suffix but are still gzipped. Detect by magic bytes.
        with open(src_path, "rb") as raw_in:
          header = raw_in.read(2)
        if header == b"\x1f\x8b":
          dst_path = dst_dir / f"{safe_id}.cif"
          with gzip.open(src_path, "rb") as f_in, open(dst_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
          return dst_path

        ext = Path(src_path).suffix.lower()
        if ext not in {".cif", ".pdb"}:
            ext = ".cif"
        dst_path = dst_dir / f"{safe_id}{ext}"
        shutil.copy2(src_path, dst_path)
        return dst_path
    except Exception as exc:
        print(f"Warning: Could not copy structure for {design_id}: {exc}")
        return None


def _extract_coords(struct_path: Path) -> list[tuple[float, float, float]]:
    parser = MMCIFParser(QUIET=True) if struct_path.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure(struct_path.stem, str(struct_path))

    ca_coords: list[tuple[float, float, float]] = []
    all_coords: list[tuple[float, float, float]] = []

    for atom in structure.get_atoms():
        x, y, z = atom.get_coord()
        all_coords.append((float(x), float(y), float(z)))
        if atom.get_name() == "CA":
            ca_coords.append((float(x), float(y), float(z)))

    return ca_coords if ca_coords else all_coords


def _save_preview(struct_path: Path, preview_path: Path) -> bool:
    try:
        coords = _extract_coords(struct_path)
        if len(coords) < 2:
            return False

        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]

        fig, ax = plt.subplots(1, 1, figsize=(4.6, 3.4), dpi=140)

        # Single "cartoon-like" trace view: project to XY and color along residue index.
        cmap = plt.cm.viridis
        n = max(len(xs) - 1, 1)
        for i in range(len(xs) - 1):
          color = cmap(i / n)
          ax.plot(
            [xs[i], xs[i + 1]],
            [ys[i], ys[i + 1]],
            color=color,
            linewidth=1.2,
            alpha=0.95,
          )

        ax.scatter(xs, ys, s=2.0, c="#1c2b4a", alpha=0.35)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor("#f8fbff")
        ax.set_title("Backbone Cartoon", fontsize=9)
        for spine in ax.spines.values():
          spine.set_alpha(0.35)

        fig.tight_layout(pad=0.6)
        fig.savefig(preview_path)
        plt.close(fig)
        return True
    except Exception as exc:
        print(f"Warning: Could not create preview for {struct_path.name}: {exc}")
        return False


def generate_dashboard(campaign_dir: str, docs_root: str) -> None:
    index_path = Path(campaign_dir) / "designs" / "index.json"
    if not index_path.is_file():
        raise FileNotFoundError(f"Missing index: {index_path}. Run sync_designs_local.py first.")

    with open(index_path, "r", encoding="utf-8") as handle:
        index = json.load(handle)

    out_page_dir = Path(docs_root) / "fgfr2"
    out_data_dir = out_page_dir / "data"
    out_structures_dir = out_page_dir / "structures"
    out_previews_dir = out_page_dir / "previews"

    out_data_dir.mkdir(parents=True, exist_ok=True)
    out_structures_dir.mkdir(parents=True, exist_ok=True)
    out_previews_dir.mkdir(parents=True, exist_ok=True)

    for design in index.get("designs", []):
        design_id = design.get("design_id", "unknown_design")
        src_path = design.get("design_path") or design.get("source_path")
        copied = _copy_structure(src_path, out_structures_dir, design_id)
        safe_id = _safe_name(design_id)

        design["structure_file"] = f"{safe_id}.cif"
        design["preview_file"] = f"{safe_id}.png"

        if copied is not None:
            # Keep extension accurate if it is not CIF.
            design["structure_file"] = copied.name
            _save_preview(copied, out_previews_dir / f"{safe_id}.png")

    data_json_path = out_data_dir / "index.json"
    with open(data_json_path, "w", encoding="utf-8") as handle:
        json.dump(index, handle, indent=2)

    html = """<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>FGFR2 Campaign Dashboard</title>
  <style>
    :root {
      --bg0: #070d18;
      --bg1: #0f172a;
      --panel: #152135;
      --line: #223149;
      --ink: #e8edf7;
      --muted: #91a2be;
      --blue: #5aa2ff;
      --pink: #ff4d7d;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: Segoe UI, sans-serif;
      color: var(--ink);
      background: linear-gradient(180deg, var(--bg1), var(--bg0));
      min-height: 100vh;
    }
    .wrap { max-width: 1600px; margin: 0 auto; padding: 24px 20px 36px; }
    .hero {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 16px;
      padding: 20px;
    }
    h1 { margin: 0; font-size: 36px; }
    .subtitle { margin-top: 8px; color: var(--muted); font-size: 17px; }
    .stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px; margin-top: 16px; }
    .card { background: #19263b; border: 1px solid var(--line); border-radius: 10px; padding: 14px; }
    .card .v { font-size: 30px; font-weight: 700; color: var(--blue); }
    .card .k { margin-top: 4px; color: #9fb0cc; font-size: 13px; }

    .tabs { display: flex; gap: 12px; margin-top: 14px; border-bottom: 1px solid var(--line); overflow-x: auto; }
    .tab {
      border: 0;
      background: transparent;
      color: #9ab0d2;
      cursor: pointer;
      padding: 9px 4px 10px;
      font-weight: 700;
      border-bottom: 3px solid transparent;
      white-space: nowrap;
    }
    .tab.active { color: #ffffff; border-bottom-color: var(--pink); }
    .tab .badge {
      margin-left: 6px;
      font-size: 12px;
      padding: 1px 7px;
      border-radius: 999px;
      background: #2a374f;
      color: #c3d2ec;
    }

    .design-list { margin-top: 18px; display: grid; gap: 14px; }
    .design {
      background: #151f31;
      border: 1px solid var(--line);
      border-radius: 12px;
      overflow: hidden;
    }
    .design-head {
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 10px;
      padding: 12px 14px;
      background: #1e2b42;
      border-bottom: 1px solid var(--line);
    }
    .design-id { font-size: 17px; font-weight: 700; }
    .chips { display: flex; flex-wrap: wrap; gap: 6px; justify-content: flex-end; }
    .chip {
      padding: 3px 9px;
      border-radius: 999px;
      font-size: 12px;
      font-weight: 700;
      color: #fff;
    }
    .rank { background: var(--pink); }
    .score { background: #4f73ff; }
    .tool { background: #29c760; color: #08210f; }
    .strategy { background: #1ec8cf; color: #052428; }

    .design-body { display: grid; grid-template-columns: 1.05fr 1fr; gap: 16px; padding: 14px; }
    .media {
      border: 1px solid #2b3b56;
      border-radius: 10px;
      min-height: 250px;
      background: #0f1626;
      display: grid;
      grid-template-rows: 1fr auto;
      overflow: hidden;
    }
    .media img {
      width: 100%;
      height: 100%;
      object-fit: contain;
      background: #f5f8ff;
    }
    .media-actions {
      border-top: 1px solid #2b3b56;
      padding: 8px 10px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 10px;
      color: #9cb0d0;
      font-size: 12px;
    }
    .download-btn {
      display: inline-block;
      padding: 6px 10px;
      border-radius: 7px;
      text-decoration: none;
      color: #fff;
      background: #4f73ff;
      font-size: 12px;
      font-weight: 700;
    }

    .metrics { width: 100%; border-collapse: collapse; }
    .metrics td {
      padding: 7px 0;
      border-bottom: 1px solid #23334d;
      font-size: 14px;
    }
    .metrics td:first-child { color: #a4b7d4; }
    .metrics td:last-child { text-align: right; color: #62e29b; font-weight: 700; }
    .warn { margin-top: 9px; color: #ffcd82; font-size: 12px; line-height: 1.4; }

    .table-wrap {
      margin-top: 18px;
      background: #121a2a;
      border: 1px solid var(--line);
      border-radius: 10px;
      overflow: auto;
    }
    table { border-collapse: collapse; width: 100%; min-width: 980px; }
    th, td { padding: 9px 10px; border-bottom: 1px solid #23334d; text-align: left; }
    th { background: #1a263a; color: #abc1e1; font-size: 12px; text-transform: uppercase; letter-spacing: 0.03em; }
    tr:hover td { background: #172235; }
    .path {
      max-width: 390px;
      overflow: hidden;
      text-overflow: ellipsis;
      white-space: nowrap;
      color: #97abcc;
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      font-size: 12px;
    }
    .footer { color: #8ea4c7; font-size: 12px; margin-top: 10px; }

    @media (max-width: 1050px) {
      h1 { font-size: 30px; }
      .design-body { grid-template-columns: 1fr; }
      .media { min-height: 210px; }
    }
  </style>
</head>
<body>
  <div class=\"wrap\">
    <div class=\"hero\">
      <h1 id=\"pageTitle\">FGFR2 Campaign Dashboard</h1>
      <div class=\"subtitle\" id=\"pageSubtitle\">Loading...</div>
      <div id=\"generatedAt\" style=\"margin-top:8px;color:#8ea4c8;font-size:13px;\">Generated: -</div>
      <div class=\"stats\" id=\"stats\"></div>
      <div class=\"tabs\" id=\"tabs\"></div>
    </div>

    <div class=\"design-list\" id=\"designList\"></div>

    <div class=\"table-wrap\">
      <table>
        <thead>
          <tr>
            <th>Design</th>
            <th>Strategy</th>
            <th>Tool</th>
            <th>ipTM</th>
            <th>ipSAE</th>
            <th>pTM</th>
            <th>Min PAE</th>
            <th>RMSD</th>
            <th>Score</th>
            <th>Path</th>
          </tr>
        </thead>
        <tbody id=\"tableBody\"></tbody>
      </table>
    </div>
    <div class=\"footer\">Offline previews: docs/fgfr2/previews | Structures: docs/fgfr2/structures</div>
  </div>

  <script>
    const toNum = (v) => {
      const n = Number(v);
      return Number.isFinite(n) ? n : null;
    };

    const fmt = (v, digits = 3) => {
      const n = toNum(v);
      return n === null ? 'NA' : n.toFixed(digits);
    };

    const fmt0 = (v) => {
      const n = toNum(v);
      return n === null ? 'NA' : String(Math.round(n));
    };

    const titleize = (s) => {
      if (!s) return 'Unknown';
      return String(s).replace(/[_-]+/g, ' ').replace(/\\b\\w/g, c => c.toUpperCase());
    };

    let active = 'all';
    let designs = [];
    let filters = [];

    const titleEl = document.getElementById('pageTitle');
    const subtitleEl = document.getElementById('pageSubtitle');
    const generatedEl = document.getElementById('generatedAt');
    const statsEl = document.getElementById('stats');
    const tabsEl = document.getElementById('tabs');
    const listEl = document.getElementById('designList');
    const tableBodyEl = document.getElementById('tableBody');

    const scoreValue = (d) => {
      const s = toNum(d.specificity_adjusted_score);
      if (s !== null) return s;
      return toNum(d.composite_score);
    };

    const ipsaeValue = (d) => {
      const m = d.metrics || {};
      return toNum(m.ipSAE ?? m.target_ipSAE ?? m.ipsae);
    };

    const renderCards = (items) => {
      const top = items.slice(0, 12);
      if (!top.length) {
        listEl.innerHTML = '<div class="design"><div class="design-head"><div class="design-id">No designs in this view.</div></div></div>';
        return;
      }

      listEl.innerHTML = top.map((d, i) => {
        const m = d.metrics || {};
        const s = scoreValue(d);
        const ipsae = ipsaeValue(d);
        const warnings = (d.warnings || []).concat(((d.evaluation || {}).warnings || []));
        const note = warnings.length ? warnings.slice(0, 2).join(' | ') : 'No warnings.';
        const strategy = (d.strategy || 'n/a').replace(/_/g, ' ');
        const preview = d.preview_file ? `./previews/${d.preview_file}` : '';
        const struct = d.structure_file ? `./structures/${d.structure_file}` : '#';

        return `
          <section class="design">
            <div class="design-head">
              <div class="design-id">${d.design_id || 'unknown_design'}</div>
              <div class="chips">
                <span class="chip rank">#${i + 1}</span>
                <span class="chip score">ipSAE: ${ipsae === null ? 'NA' : ipsae.toFixed(3)}</span>
                <span class="chip tool">${d.tool || 'tool'}</span>
                <span class="chip strategy">${strategy}</span>
              </div>
            </div>
            <div class="design-body">
              <div class="media">
                ${preview ? `<img src="${preview}" alt="Structure preview for ${d.design_id}" />` : '<div style="padding:14px;color:#9db0ce;">Preview unavailable</div>'}
                <div class="media-actions">
                  <span>${d.structure_file || 'structure missing'}</span>
                  <a href="${struct}" class="download-btn" download>Download</a>
                </div>
              </div>
              <div>
                <table class="metrics">
                  <tr><td>ipTM</td><td>${fmt(m.ipTM)}</td></tr>
                  <tr><td>ipSAE</td><td>${fmt(ipsae)}</td></tr>
                  <tr><td>pTM</td><td>${fmt(m.pTM)}</td></tr>
                  <tr><td>RMSD</td><td>${fmt(m.rmsd)}</td></tr>
                  <tr><td>Min PAE</td><td>${fmt(m.min_interaction_pae)}</td></tr>
                  <tr><td>Binder Length</td><td>${fmt0(m.binder_length)}</td></tr>
                </table>
                <div class="warn">${note}</div>
              </div>
            </div>
          </section>
        `;
      }).join('');
    };

    const renderTable = (items) => {
      tableBodyEl.innerHTML = items.slice(0, 400).map((d) => {
        const m = d.metrics || {};
        const s = scoreValue(d);
        const ipsae = ipsaeValue(d);
        const path = d.design_path || d.source_path || 'NA';
        return `
          <tr>
            <td>${d.design_id || 'NA'}</td>
            <td>${(d.strategy || 'NA').replace(/_/g, ' ')}</td>
            <td>${d.tool || 'NA'}</td>
            <td>${fmt(m.ipTM)}</td>
            <td>${fmt(ipsae)}</td>
            <td>${fmt(m.pTM)}</td>
            <td>${fmt(m.min_interaction_pae)}</td>
            <td>${fmt(m.rmsd)}</td>
            <td>${s === null ? 'NA' : s.toFixed(3)}</td>
            <td><div class="path" title="${path}">${path}</div></td>
          </tr>
        `;
      }).join('');
    };

    const render = () => {
      tabsEl.innerHTML = filters
        .map((f) => `<button class="tab ${f.id === active ? 'active' : ''}" data-id="${f.id}">${f.label}<span class="badge">${f.data.length}</span></button>`)
        .join('');

      const selected = (filters.find(f => f.id === active) || { data: designs }).data;
      renderCards(selected);
      renderTable(selected);
    };

    tabsEl.addEventListener('click', (ev) => {
      const btn = ev.target.closest('button[data-id]');
      if (!btn) return;
      active = btn.dataset.id;
      render();
    });

    async function loadDashboard() {
      try {
        const res = await fetch('./data/index.json?ts=' + Date.now(), { cache: 'no-store' });
        if (!res.ok) throw new Error('Failed to load data/index.json');
        const data = await res.json();

        titleEl.textContent = `${titleize(data.campaign || 'campaign')} Campaign Dashboard`;
        subtitleEl.textContent = `${data.num_designs || 0} computational binder candidates | ranked by ipSAE`;
        generatedEl.textContent = `Generated: ${new Date(data.generated_at).toLocaleString()}`;

        designs = (data.designs || []).slice();

        const toolMap = {};
        const stratMap = {};
        designs.forEach((d) => {
          const tool = titleize(d.tool || 'unknown');
          if (!toolMap[tool]) toolMap[tool] = [];
          toolMap[tool].push(d);

          const strategy = titleize(d.strategy || 'unknown');
          if (!stratMap[strategy]) stratMap[strategy] = [];
          stratMap[strategy].push(d);
        });

        const scoreNums = designs.map(scoreValue).filter((x) => x !== null);
        const avgScore = scoreNums.length ? scoreNums.reduce((a, b) => a + b, 0) / scoreNums.length : null;
        const iptmNums = designs.map((d) => toNum((d.metrics || {}).ipTM)).filter((x) => x !== null);
        const avgIptm = iptmNums.length ? iptmNums.reduce((a, b) => a + b, 0) / iptmNums.length : null;
        const ipsaeNums = designs.map((d) => ipsaeValue(d)).filter((x) => x !== null);
        const avgIpsae = ipsaeNums.length ? ipsaeNums.reduce((a, b) => a + b, 0) / ipsaeNums.length : null;

        statsEl.innerHTML = [
          { v: designs.length, k: 'Total Designs' },
          { v: Object.keys(toolMap).length, k: 'Tools' },
          { v: Object.keys(stratMap).length, k: 'Strategies' },
          { v: avgScore === null ? 'NA' : avgScore.toFixed(3), k: 'Avg Score' },
          { v: avgIptm === null ? 'NA' : avgIptm.toFixed(2), k: 'Avg ipTM' },
          { v: avgIpsae === null ? 'NA' : avgIpsae.toFixed(2), k: 'Avg ipSAE' },
        ].map(({ v, k }) => `<div class="card"><div class="v">${v}</div><div class="k">${k}</div></div>`).join('');

        filters = [
          { id: 'all', label: 'All', data: designs },
          { id: 'boltzgen', label: 'BoltzGen', data: toolMap['Boltzgen'] || [] },
          { id: 'rfdiffusion3', label: 'RFdiffusion3', data: toolMap['Rfdiffusion3'] || [] },
          ...Object.entries(stratMap).map(([k, v]) => ({ id: k.toLowerCase().replace(/ /g, '-'), label: k, data: v })),
          { id: 'table', label: 'Table', data: designs },
        ].filter((f) => f.data.length > 0);

        active = 'all';
        render();
      } catch (err) {
        console.error(err);
        subtitleEl.textContent = 'Error loading dashboard data';
        subtitleEl.style.color = '#ff6b6b';
      }
    }

    document.addEventListener('DOMContentLoaded', loadDashboard);
  </script>
</body>
</html>
"""

    html_path = out_page_dir / "index.html"
    with open(html_path, "w", encoding="utf-8") as handle:
        handle.write(html)

    print(f"Wrote dashboard HTML: {html_path}")
    print(f"Wrote dashboard data JSON: {data_json_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate local FGFR2 dashboard from design index.")
    parser.add_argument("--campaign-dir", default="fgfr2_campaign", help="Path to campaign directory")
    parser.add_argument("--docs-root", default="docs", help="Path to docs output root")
    args = parser.parse_args()

    generate_dashboard(args.campaign_dir, args.docs_root)
    print("Dashboard generation complete!")


if __name__ == "__main__":
    main()
