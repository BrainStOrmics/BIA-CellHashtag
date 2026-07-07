"""Human-in-the-Loop interaction for annotation review."""

from pathlib import Path
from typing import TypedDict


class AnnotationReviewResult(TypedDict):
    cluster: str
    original_cell_type: str
    reviewed_cell_type: str
    original_confidence: float
    reviewed_confidence: float
    approved: bool


def review_low_confidence_annotations(
    results: list[dict],
    threshold: float = 0.5,
) -> list[AnnotationReviewResult]:
    """Interactively review low-confidence annotations.

    Returns updated results list. If no TTY available, skip silently.
    """
    import sys

    if not sys.stdin.isatty():
        return results

    low_conf = [r for r in results if r.get("confidence", 0) < threshold]
    if not low_conf:
        return results

    print("\n" + "=" * 60)
    print("🔍 HITL: Low-Confidence Annotation Review")
    print("=" * 60)
    print(f"{len(low_conf)} cluster(s) below confidence {threshold:.2f}")

    updated = []
    for r in results:
        if r.get("confidence", 0) >= threshold:
            updated.append(r)
            continue

        print(f"\n--- Cluster {r['cluster']} ---")
        print(f"  Suggested: {r['cell_type']} (confidence: {r['confidence']:.2f})")
        print(f"  Reasoning: {r.get('reasoning', 'N/A')}")
        print(f"\n  [a] Accept as-is")
        print(f"  [r] Reject → mark as Unknown")
        print(f"  [e] Edit cell type (enter custom name)")
        choice = input("\n  Choice [a/r/e]: ").strip().lower()

        if choice == "a":
            updated.append(r)
        elif choice == "r":
            updated.append({**r, "cell_type": "Unknown", "confidence": 0.0, "reasoning": "Rejected by user"})
        elif choice == "e":
            new_type = input("  New cell type: ").strip()
            if new_type:
                updated.append({**r, "cell_type": new_type, "confidence": max(r["confidence"], 0.5)})
            else:
                updated.append(r)
        else:
            updated.append(r)

    return updated


def build_review_report(results: list[dict], output_dir: Path):
    """Write HITL review summary to output_dir."""
    low = [r for r in results if r.get("confidence", 0) < 0.5]
    unknown = [r for r in results if r.get("cell_type") == "Unknown"]
    md = f"""# HITL Review Summary

Low confidence (<0.5): {len(low)} clusters
"""
    if low:
        md += "\n"
        for r in low:
            md += f"- Cluster {r['cluster']}: {r['cell_type']} ({r['confidence']:.2f}) — {r.get('reasoning', '')}\n"
    if unknown:
        md += f"\nUnknown annotations: {len(unknown)} clusters\n"
        for r in unknown:
            md += f"- Cluster {r['cluster']}\n"
    (output_dir / "hitl_review.md").write_text(md)
