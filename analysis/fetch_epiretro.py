"""Inputs for D-104: Epi-retro-seq labels and the Allen whole-brain atlas pieces they attach to.

    python3 analysis/fetch_epiretro.py

Downloads, each skipped when its output exists (data/raw/ is not version controlled):

data/raw/epiretro/  Zhou et al. 2023, Nature 624:355-365, doi:10.1038/s41586-023-06823-w
                    (article and supplements CC BY 4.0)
  st2_cells.xlsx      Supplementary Table 2: one row per sequenced nucleus, with source,
                      target, sex, the target filter and the Allen subclass

data/raw/abc/       Allen Brain Cell Atlas (Yao et al. 2023; Allen data CC BY-NC 4.0)
  precomputed_stats.h5     per-cluster mean log2(CPM+1) over the whole-brain 10x data
                           (the MapMyCells reference; its n_cells field is 1 throughout)
  term_counts.csv          cells per cluster, to weight cluster means into subclass means
  merfish_cells.csv        MERFISH C57BL6J-638850 cells with cluster, CCF coordinates and
                           CCF parcellation
  genes_10x.csv            Ensembl identifier -> gene symbol for the 10x data
  genes_merfish.csv        the MERFISH panel's genes

The CCF ontology (data/raw/atlas/structure_graph_1.json) is fetched by src/syn4rch/atlas/build.py.

Prints each file's size and SHA-256, which the analysis pins.
"""
import hashlib
import os
import sys
import urllib.request

SPRINGER = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06823-w/"
            "MediaObjects/41586_2023_6823_{}_ESM.xlsx")
ABC = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/"
FILES = [
    ("data/raw/epiretro/st2_cells.xlsx", SPRINGER.format("MOESM4")),
    ("data/raw/abc/precomputed_stats.h5",
     ABC + "mapmycells/WMB-10X/20240831/precomputed_stats_ABC_revision_230821.h5"),
    ("data/raw/abc/merfish_cells.csv",
     ABC + "metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/"
           "cell_metadata_with_parcellation_annotation.csv"),
    ("data/raw/abc/term_counts.csv",
     ABC + "metadata/WMB-taxonomy/20231215/views/cluster_annotation_term_with_counts.csv"),
    ("data/raw/abc/genes_10x.csv", ABC + "metadata/WMB-10X/20231215/gene.csv"),
    ("data/raw/abc/genes_merfish.csv", ABC + "metadata/MERFISH-C57BL6J-638850/20231215/gene.csv"),
]
HEADERS = {"User-Agent": "syn4rch-analysis/0.1"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(dest, url):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        return
    part = dest + ".part"
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=300) as r, open(part, "wb") as f:
        while True:
            block = r.read(1 << 22)
            if not block:
                break
            f.write(block)
    os.replace(part, dest)


def main():
    for dest, url in FILES:
        fetch(dest, url)
        print("%10.1f MB  %s  %s" % (os.path.getsize(dest) / 1e6, sha256(dest), dest))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
