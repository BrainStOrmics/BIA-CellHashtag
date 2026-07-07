"""Basic CellHashtag usage example.

Run with: python examples/basic_usage.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from cellhashtag import CellHashtagAgent

agent = CellHashtagAgent(profile="fast")

result = agent.run(
    input_path="data/example.h5ad",
    output_dir="output",
    cluster_key="leiden",
)

print(f"Status: {result['status']}")
print(f"Clusters: {result['n_clusters']}")
if result["errors"]:
    print(f"Errors: {result['errors']}")
