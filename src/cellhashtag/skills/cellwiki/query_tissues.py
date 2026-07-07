"""Query CellWiki tissue directories for plausible cell types."""

from pathlib import Path
from typing import Optional


def get_tissue_cell_types(tissue: str, wiki_dir: Optional[Path] = None) -> list[str]:
    if wiki_dir is None:
        wiki_dir = Path.home() / "CellWiki" / "wiki"
    tissue_dir = wiki_dir / "tissues"
    if not tissue_dir.exists():
        return []
    cell_types = set()
    for tf in tissue_dir.glob("*.md"):
        content = tf.read_text().lower()
        if tissue.lower() in content:
            cell_types.add(tf.stem)
    return sorted(cell_types)
