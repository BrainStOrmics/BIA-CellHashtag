"""Generate a mock h5ad for testing. 300 cells, 4 clusters, realistic PBMC markers."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

np.random.seed(42)

n_cells = 300

immune = ["CD3D", "CD3E", "CD4", "CD8A", "CD8B", "CD19", "MS4A1", "CD79A", "CD79B", "CD14", "FCGR3A", "CD68", "ITGAM"]
structural = ["VIM", "COL1A1", "COL3A1", "DCN", "LUM", "PECAM1", "CDH5", "VWF"]
epithelial = ["EPCAM", "KRT18", "KRT19", "MUC1", "CDH1"]
housekeeping = ["ACTB", "GAPDH", "B2M", "RPL13A", "RPS18", "EEF1A1", "TUBB"]
mt = ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ND4", "MT-CYB"]
ribo = ["RPS3", "RPS6", "RPL5", "RPL8"]
hb = ["HBA1", "HBA2", "HBB", "HBD"]
specific = immune + structural + epithelial + housekeeping + mt + ribo + hb
n_random = 150 - len(specific)
gene_names = specific + [f"GENE_{i:03d}" for i in range(n_random)]
n_genes = len(gene_names)

clusters = np.random.choice([0, 1, 2, 3], size=n_cells, p=[0.35, 0.25, 0.25, 0.15])
X = np.random.exponential(0.3, size=(n_cells, n_genes))

markers_by_cluster = {
    0: ["CD3D", "CD3E", "CD4"],
    1: ["CD19", "MS4A1", "CD79A"],
    2: ["CD14", "CD68", "FCGR3A"],
    3: ["VIM", "COL1A1", "DCN"],
}
for cid, markers in markers_by_cluster.items():
    mask = clusters == cid
    for m in markers:
        idx = gene_names.index(m)
        X[mask, idx] += np.random.exponential(3, mask.sum())

X = csr_matrix(X.astype(np.float32))

obs = pd.DataFrame({
    "Organism": ["Homo sapiens"] * n_cells,
    "tissue": ["PBMC"] * n_cells,
    "disease": ["healthy"] * n_cells,
    "leiden": pd.Categorical([str(c) for c in clusters]),
}, index=[f"cell_{i:04d}" for i in range(n_cells)])

var = pd.DataFrame(index=gene_names)
var["mt"] = [g.startswith("MT-") for g in gene_names]
var["ribo"] = [g.startswith("RPS") or g.startswith("RPL") for g in gene_names]
var["hb"] = [g.startswith("HB") for g in gene_names]
var["highly_variable"] = True

adata = sc.AnnData(X=X, obs=obs, var=var)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_pcs=20)
sc.tl.umap(adata)
sc.tl.rank_genes_groups(adata, groupby="leiden", reference="rest", method="wilcoxon")

out = Path(__file__).parent.parent / "data" / "example.h5ad"
adata.write_h5ad(out)
print(f"Created {out}: {adata.n_obs} cells x {adata.n_vars} genes, {adata.obs['leiden'].nunique()} clusters")
print(f"Cluster sizes: {adata.obs['leiden'].value_counts().to_dict()}")
