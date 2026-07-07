"""Parse CellWiki cell type markdown files and query marker genes."""

import re
from pathlib import Path
from typing import Optional


def _parse_frontmatter(content: str) -> dict:
    m = re.match(r"^---\n(.*?)\n---", content, re.DOTALL)
    if not m:
        return {}
    result = {}
    for line in m.group(1).splitlines():
        if ":" in line:
            key, val = line.split(":", 1)
            key = key.strip()
            val = val.strip()
            if val.startswith("[") and val.endswith("]"):
                items = [x.strip().strip("-").strip() for x in val.strip("[]").split(",") if x.strip()]
                result[key] = items
            else:
                result[key] = val
    return result


def _parse_markers_table(content: str) -> list[dict]:
    markers = []
    lines = content.splitlines()
    in_table = False
    for line in lines:
        if line.startswith("| Marker |"):
            in_table = True
            continue
        if in_table and line.startswith("|"):
            parts = [p.strip() for p in line.split("|")[1:-1]]
            if len(parts) >= 4 and not all(c.startswith("---") for c in parts):
                markers.append({
                    "gene": parts[0],
                    "type": parts[1],
                    "evidence": parts[2],
                    "source": parts[3],
                })
    return markers


def _parse_contexts(content: str) -> dict:
    contexts = {"tissues": [], "diseases": []}
    tissue_match = re.search(r"\*\*Tissues:\*\*\s*(.+)", content)
    if tissue_match:
        contexts["tissues"] = [t.strip() for t in tissue_match.group(1).split(",")]
    disease_match = re.search(r"\*\*Diseases:\*\*\s*(.+)", content)
    if disease_match:
        contexts["diseases"] = [d.strip() for d in disease_match.group(1).split(",")]
    return contexts


def query_cell_type(name: str, wiki_dir: Optional[Path] = None) -> Optional[dict]:
    if wiki_dir is None:
        wiki_dir = Path.home() / "CellWiki" / "wiki"
    path = wiki_dir / "cell_types" / f"{name}.md"
    if not path.exists():
        return None
    content = path.read_text()
    frontmatter = _parse_frontmatter(content)
    markers = _parse_markers_table(content)
    contexts = _parse_contexts(content)
    return {
        "name": frontmatter.get("standard_name", name),
        "display_name": frontmatter.get("display_name", name),
        "cl_id": frontmatter.get("cl_id", ""),
        "parent_type": frontmatter.get("parent_type", ""),
        "markers": markers,
        "tissues": contexts["tissues"],
        "diseases": contexts["diseases"],
    }


def list_cell_types(wiki_dir: Optional[Path] = None) -> list[str]:
    if wiki_dir is None:
        wiki_dir = Path.home() / "CellWiki" / "wiki"
    ct_dir = wiki_dir / "cell_types"
    if not ct_dir.exists():
        return []
    return [p.stem for p in ct_dir.glob("*.md")]


def query_markers_for_genes(genes: list[str], wiki_dir: Optional[Path] = None) -> list[dict]:
    if wiki_dir is None:
        wiki_dir = Path.home() / "CellWiki" / "wiki"
    gene_set = {g.upper() for g in genes}
    results = []
    for ct_name in list_cell_types(wiki_dir):
        info = query_cell_type(ct_name, wiki_dir)
        if info is None:
            continue
        ct_markers = {m["gene"].upper() for m in info["markers"]}
        overlap = gene_set & ct_markers
        if not overlap:
            continue
        all_markers_set = gene_set | ct_markers
        jaccard = len(overlap) / len(all_markers_set) if all_markers_set else 0.0
        results.append({
            "cell_type": info["display_name"],
            "cl_id": info["cl_id"],
            "overlap": sorted(overlap),
            "jaccard": round(jaccard, 4),
            "markers_all": sorted(ct_markers),
            "parent_type": info["parent_type"],
        })
    return sorted(results, key=lambda x: x["jaccard"], reverse=True)
