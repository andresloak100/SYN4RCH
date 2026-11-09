"""Mean expression of every gene for every HypoMap cell type -- no gene chosen by hand.

The input to D-091's blind fingerprint, and to D-095's gene similarity.

    python3 analysis/extract_hypomap_means.py all_genes.npz

Reads data/raw/hypomap/hypomap.h5ad (384,925 cells x 51,216 genes, not version
controlled) and writes, per C465 cell type, the mean of every gene. A sparse indicator
product per chunk keeps memory bounded: sums += onehot(cluster).T @ X_chunk, never
densifying the cells x genes block. Needs numpy, scipy and anndata.
"""
import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp

a = ad.read_h5ad("data/raw/hypomap/hypomap.h5ad", backed="r")
clusters = a.obs["C465_named"].astype(str).values
uniq = sorted(set(clusters))
col = {c: i for i, c in enumerate(uniq)}
k = np.array([col[c] for c in clusters])
sums = np.zeros((len(uniq), a.n_vars))
counts = np.bincount(k, minlength=len(uniq))
CH = 20000
for start in range(0, a.n_obs, CH):
    stop = min(start + CH, a.n_obs)
    X = a.X[start:stop]
    X = X if sp.issparse(X) else sp.csr_matrix(X)
    ind = sp.csr_matrix((np.ones(stop - start), (np.arange(stop - start), k[start:stop])),
                        shape=(stop - start, len(uniq)))
    sums += (ind.T @ X).toarray()
    print("  %d/%d" % (stop, a.n_obs), file=sys.stderr, flush=True)
means = sums / counts[:, None]
# var_names are Ensembl ids; the gene symbols are in feature_name
sym = np.array([str(s) for s in a.var["feature_name"]])
np.savez_compressed(sys.argv[1], means=means, clusters=np.array(uniq), genes=sym, counts=counts)
print("wrote %s: %d cell types x %d genes" % (sys.argv[1], *means.shape), file=sys.stderr)
