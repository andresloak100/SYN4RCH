"""Independent-animal cohort for the real 3D axon benchmark (Jiao et al. 2025): inventory and split.

    python3 -I analysis/axon_dataset_audit.py

METADATA ONLY. It reads the portal's neuron index (allowed fields), the list of traces already on
disk, the soma / label tables, and the per-animal archive sizes recorded in
reports/axon_dataset_sources.md. It never opens a trace (.swc) and never opens the regional outcome
tables; a guard on `open` enforces that. It writes:

    reports/axon_dataset_manifest.tsv          one row per neuron (fields from the CC BY tables + our flags)
    reports/axon_dataset_animals.tsv           one row per animal (part, counts, archive bytes)
    data/interim/axon_dataset/local_manifest.tsv   the same plus hemisphere and reconstruction type
    data/interim/axon_dataset/counts.json      every number quoted in reports/axon_dataset_plan.md

The split is by whole animal, and it cannot be tuned: an animal's place is decided by the SHA-256 of
SALT + its ID, inside its source label, with the quotas in `quota()`. SALT was written before the
split was first computed and is not to be changed.
"""
import builtins
import collections
import csv
import hashlib
import io
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = os.path.join(ROOT, "data/raw/jiao2025")
INDEX = os.path.join(J, "swc/portal_index/mouse.neuron.info.json")
INSPECTED = os.path.join(J, "swc/MANIFEST.tsv")
SOMA_XLSX = os.path.join(J, "all_neuron_soma_infor.xlsx")
SOURCES = os.path.join(ROOT, "reports/axon_dataset_sources.md")
OUT_LOCAL = os.path.join(ROOT, "data/interim/axon_dataset")
FORBIDDEN = ("all_neurons_projection_length_fig1.csv", "all_neurons_terminal_terminal_number_Supplementary_fig3_sort.csv")
INDEX_FIELDS = ("file", "acronym", "type_name", "hemisphere", "reconstruction_type")   # never type / class
SOMA_COLS = ("Unnamed: 0", "neurons", "location", "geno", "AP", "DV", "ML")            # never cluster / Class

SALT = "SYN4RCH-axon-cohort-2026-11-10"
MIN_NEURONS_HELD_OUT = 5          # an animal needs this many neurons to be a validation or test animal
SHARE_TEST, SHARE_VALIDATION = 0.30, 0.20
LI_SHARED = ("195037", "195038")  # animals also in Li et al. 2024: training only
# The other agent's development neurons, as recorded in reports/d109_protocol.md section 16.
DEVELOPMENT = ("201230_079 201232_040 201232_056 201894_001 202170_016 202559_043 202617_025 202620_016 "
               "202824_062 210007_142 210104_007 210796_002 211353_070 211889_096 212202_011 212202_026").split()

_open = builtins.open


def guarded_open(file, *a, **k):
    name = str(file)
    if name.endswith(".swc") or name.endswith(FORBIDDEN):
        raise PermissionError("the cohort audit must not open %s" % name)
    return _open(file, *a, **k)


builtins.open = io.open = guarded_open


def half_up(x):
    return int(x + 0.5)


def quota(m):
    """(test, validation) animals for a source label with m eligible new animals."""
    if m >= 10:
        return half_up(SHARE_TEST * m), half_up(SHARE_VALIDATION * m)
    if m >= 4:
        return max(1, half_up(SHARE_TEST * m)), max(1, half_up(SHARE_VALIDATION * m))
    if m >= 2:
        return 1, 1
    return (1, 0) if m == 1 else (0, 0)


def order_key(animal):
    return hashlib.sha256((SALT + ":" + animal).encode()).hexdigest()


def archive_bytes():
    """Per-animal BSDC archive sizes (exact, from the record's public file list) as tabulated in
    the appendix of reports/axon_dataset_sources.md, with that table's own counts for a cross-check."""
    out = {}
    for line in open(SOURCES):
        m = re.match(r"^\| (\d{5,6}) \| ([^|]+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| ([^|]+) \| ([\d,]+) \|", line)
        if m:
            out[m.group(1)] = {"label": m.group(2).strip(), "neurons": int(m.group(3)), "left": int(m.group(4)),
                               "right": int(m.group(5)), "axon_only": int(m.group(6)), "status": m.group(7).strip(),
                               "zip_bytes": int(m.group(8).replace(",", ""))}
    return out


def main():
    import numpy as np
    import pandas as pd

    idx = json.load(open(INDEX))["neuron_data"]
    N = {}
    for r in idx.values():
        nid = r["file"][:-4]
        N[nid] = {k: r[k] for k in INDEX_FIELDS}
    assert len(N) == 7180
    soma = pd.read_excel(SOMA_XLSX, usecols=lambda c: c in SOMA_COLS).rename(columns={"Unnamed: 0": "name"})
    soma["id"] = soma["neurons"].astype(str)
    soma = soma.set_index("id")
    assert set(soma.index) == set(N) and soma.index.is_unique
    norm = lambda g: "Orexin" if g == "Wt(orexin)" else g
    for nid, r in N.items():
        s = soma.loc[nid]
        r.update(animal=nid.split("_")[0], label=norm(s["geno"]), soma_label=s["location"],
                 ap=float(s["AP"]), dv=float(s["DV"]), ml=float(s["ML"]))
        assert r["label"] == r["type_name"], (nid, r["label"], r["type_name"])

    inspected = {m["neuron"] for m in csv.DictReader(open(INSPECTED), delimiter="\t") if re.match(r"^\d+_\d+$", m.get("neuron") or "")}
    dev = set(DEVELOPMENT)
    assert len(inspected) == 37 and dev <= inspected and inspected <= set(N)
    a_insp = {i.split("_")[0] for i in inspected}
    a_dev = {i.split("_")[0] for i in dev}

    A = collections.defaultdict(list)
    for nid, r in N.items():
        A[r["animal"]].append(nid)
    label = {}
    for a, ids in A.items():
        ls = {N[i]["label"] for i in ids}
        assert len(ls) == 1, (a, ls)
        label[a] = ls.pop()
    assert len(A) == 166

    # ---- the split -------------------------------------------------------------------------------
    part, reason = {}, {}
    for a in A:
        if a in a_dev:
            part[a], reason[a] = "train", "development animal (the 16 development neurons): never held out"
        elif a in a_insp:
            part[a], reason[a] = "train", "a trace of this animal was inspected before the plan: never held out"
        elif a in LI_SHARED:
            part[a], reason[a] = "train", "animal also in Li et al. 2024: never held out"
        elif len(A[a]) < MIN_NEURONS_HELD_OUT:
            part[a], reason[a] = "train", "new animal with fewer than %d neurons: never held out" % MIN_NEURONS_HELD_OUT
    eligible = collections.defaultdict(list)
    for a in A:
        if a not in part:
            eligible[label[a]].append(a)
    quotas = {}
    for L, animals in sorted(eligible.items()):
        animals.sort(key=order_key)
        nt, nv = quota(len(animals))
        quotas[L] = {"eligible": len(animals), "test": nt, "validation": nv, "train": len(animals) - nt - nv}
        for k, a in enumerate(animals):
            if k < nt:
                part[a], reason[a] = "test", "new animal; hash rank %d of %d in %s" % (k + 1, len(animals), L)
            elif k < nt + nv:
                part[a], reason[a] = "validation", "new animal; hash rank %d of %d in %s" % (k + 1, len(animals), L)
            else:
                part[a], reason[a] = "train", "new animal; hash rank %d of %d in %s" % (k + 1, len(animals), L)
    assert set(part) == set(A)
    held = {a for a in A if part[a] != "train"}
    assert not held & a_insp and not held & set(LI_SHARED)
    assert all(len(A[a]) >= MIN_NEURONS_HELD_OUT for a in held)

    def klass(nid):
        a = N[nid]["animal"]
        if nid in dev:
            return "development"
        if nid in inspected:
            return "inspected_not_development"
        return "new_neuron_of_inspected_animal" if a in a_insp else "neuron_of_new_animal"

    # ---- outputs ---------------------------------------------------------------------------------
    os.makedirs(OUT_LOCAL, exist_ok=True)
    zips = archive_bytes()
    assert set(zips) == set(A), "archive table does not list the same 166 animals"
    ids = sorted(N)
    with open(os.path.join(ROOT, "reports/axon_dataset_manifest.tsv"), "w") as f, \
            open(os.path.join(OUT_LOCAL, "local_manifest.tsv"), "w") as g:
        head = ["neuron", "animal", "source_label", "soma_label", "class", "animal_inspected", "animal_development",
                "animal_in_li_2024", "part", "part_reason"]
        f.write("\t".join(head) + "\n")
        g.write("\t".join(head + ["hemisphere", "reconstruction_type", "index_soma_acronym"]) + "\n")
        for nid in ids:
            r = N[nid]
            a = r["animal"]
            row = [nid, a, r["label"], r["soma_label"], klass(nid), str(int(a in a_insp)), str(int(a in a_dev)),
                   str(int(a in LI_SHARED)), part[a], reason[a]]
            f.write("\t".join(row) + "\n")
            g.write("\t".join(row + [r["hemisphere"], r["reconstruction_type"], r["acronym"]]) + "\n")
    with open(os.path.join(ROOT, "reports/axon_dataset_animals.tsv"), "w") as f:
        f.write("animal\tsource_label\tneurons\tpart\tstatus\tbsdc_zip_bytes\tpart_reason\n")
        for a in sorted(A):
            status = "development" if a in a_dev else "inspected" if a in a_insp else "new"
            f.write("\t".join([a, label[a], str(len(A[a])), part[a], status, str(zips[a]["zip_bytes"]), reason[a]]) + "\n")

    C = {"salt": SALT, "min_neurons_held_out": MIN_NEURONS_HELD_OUT, "share_test": SHARE_TEST,
         "share_validation": SHARE_VALIDATION, "quotas": quotas}
    cls = collections.Counter(klass(i) for i in ids)
    C["classes"] = dict(cls)
    C["animals"] = {"all": len(A), "inspected": len(a_insp), "development": len(a_dev), "new": len(set(A) - a_insp)}
    # cross-check against the appendix of the sources report (an independent tabulation)
    mism = []
    for a in A:
        z = zips[a]
        mine = (len(A[a]), sum(N[i]["hemisphere"] == "Left" for i in A[a]), sum(N[i]["hemisphere"] == "Right" for i in A[a]),
                sum(N[i]["reconstruction_type"] == "Axon_only" for i in A[a]))
        if mine != (z["neurons"], z["left"], z["right"], z["axon_only"]) or norm(z["label"]) != label[a]:
            mism.append(a)
    C["appendix_cross_check_mismatches"] = mism

    def summary(animals):
        nn = [i for a in animals for i in A[a]]
        return {"animals": len(animals), "neurons": len(nn),
                "left": sum(N[i]["hemisphere"] == "Left" for i in nn), "right": sum(N[i]["hemisphere"] == "Right" for i in nn),
                "axon_only": sum(N[i]["reconstruction_type"] == "Axon_only" for i in nn),
                "soma_labels": len({N[i]["soma_label"] for i in nn}),
                "soma_label_HY_or_tract": sum(N[i]["soma_label"] == "HY" or N[i]["soma_label"][:1].islower() for i in nn),
                "zip_bytes": sum(zips[a]["zip_bytes"] for a in animals),
                "by_label": {L: [sum(label[a] == L for a in animals), sum(len(A[a]) for a in animals if label[a] == L)]
                             for L in sorted(set(label.values()))}}

    parts = {p: sorted(a for a in A if part[a] == p) for p in ("train", "validation", "test")}
    C["parts"] = {p: summary(v) for p, v in parts.items()}
    tr = parts["train"]
    C["train_breakdown"] = {
        "development_animals": summary([a for a in tr if a in a_dev]),
        "inspected_not_development_animals": summary([a for a in tr if a in a_insp and a not in a_dev]),
        "li_shared_new_animal": summary([a for a in tr if a in LI_SHARED and a not in a_insp]),
        "small_new_animals": summary([a for a in tr if a not in a_insp and a not in LI_SHARED and len(A[a]) < MIN_NEURONS_HELD_OUT]),
        "new_animals_by_hash": summary([a for a in tr if a not in a_insp and a not in LI_SHARED and len(A[a]) >= MIN_NEURONS_HELD_OUT])}
    C["on_disk"] = {"traces": len(inspected), "development": len(dev), "inspected_not_development": len(inspected - dev),
                    "animals": len(a_insp), "in_train": sum(part[i.split("_")[0]] == "train" for i in inspected)}
    C["new_animals"] = summary(sorted(set(A) - a_insp))
    C["everything"] = summary(sorted(A))

    # donor coverage, inputs only: distance from each held-out soma to the nearest training soma
    xyz = lambda nn: np.array([[N[i]["ap"], N[i]["dv"], N[i]["ml"]] for i in nn])
    tr_ids = [i for a in tr for i in A[a]]
    T = xyz(tr_ids)
    cov = {}
    for p in ("validation", "test"):
        q_ids = [i for a in parts[p] for i in A[a]]
        Q = xyz(q_ids)
        d = np.empty(len(Q))
        for s in range(0, len(Q), 500):
            d[s:s + 500] = np.sqrt(((Q[s:s + 500, None, :] - T[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
        lab = np.array([N[i]["label"] for i in q_ids])
        same = np.empty(len(Q))
        tl = np.array([N[i]["label"] for i in tr_ids])
        for k in range(len(Q)):
            m = tl == lab[k]
            same[k] = np.sqrt(((T[m] - Q[k]) ** 2).sum(axis=1)).min() if m.any() else np.nan
        pct = lambda v: [float(np.nanpercentile(v, x)) for x in (50, 90, 100)]
        cov[p] = {"any_label_um_p50_p90_max": pct(d), "same_label_um_p50_p90_max": pct(same),
                  "same_label_missing": int(np.isnan(same).sum()),
                  "by_label_p50_any_same": {L: [float(np.median(d[lab == L])), float(np.nanmedian(same[lab == L]))]
                                            for L in sorted(set(lab))}}
    C["donor_coverage_um"] = cov
    C["manifest_sha256"] = hashlib.sha256(open(os.path.join(ROOT, "reports/axon_dataset_manifest.tsv"), "rb").read()).hexdigest()
    C["animals_sha256"] = hashlib.sha256(open(os.path.join(ROOT, "reports/axon_dataset_animals.tsv"), "rb").read()).hexdigest()
    json.dump(C, open(os.path.join(OUT_LOCAL, "counts.json"), "w"), indent=1, sort_keys=True)

    # ---- printed summary -------------------------------------------------------------------------
    p = print
    p("classes:", dict(cls))
    p("animals:", C["animals"], " appendix cross-check mismatches:", mism)
    p("quotas by source label (eligible new animals -> test / validation / train):")
    for L, q in quotas.items():
        p("  %-8s %2d -> %2d / %2d / %2d" % (L, q["eligible"], q["test"], q["validation"], q["train"]))
    for name, s in list(C["parts"].items()) + [("new animals (all)", C["new_animals"]), ("everything", C["everything"])]:
        p("%-18s animals %3d  neurons %4d  Left %4d Right %4d  Axon_only %3d  soma labels %2d (HY or tract: %d)  zip bytes %s"
          % (name, s["animals"], s["neurons"], s["left"], s["right"], s["axon_only"], s["soma_labels"],
             s["soma_label_HY_or_tract"], format(s["zip_bytes"], ",")))
    p("training pool:")
    for k, s in C["train_breakdown"].items():
        p("  %-36s animals %3d  neurons %4d  zip bytes %s" % (k, s["animals"], s["neurons"], format(s["zip_bytes"], ",")))
    p("per source label, animals / neurons:   train | validation | test")
    for L in sorted(set(label.values())):
        p("  %-8s %s" % (L, " | ".join("%2d / %4d" % tuple(C["parts"][x]["by_label"][L]) for x in ("train", "validation", "test"))))
    p("on disk:", C["on_disk"])
    for x in ("validation", "test"):
        c = cov[x]
        p("donor coverage, %s: nearest training soma, any label p50 %.0f p90 %.0f max %.0f um; same label p50 %.0f p90 %.0f max %.0f um"
          % ((x,) + tuple(c["any_label_um_p50_p90_max"]) + tuple(c["same_label_um_p50_p90_max"])))
    p("manifest sha256", C["manifest_sha256"])
    p("animals  sha256", C["animals_sha256"])


if __name__ == "__main__":
    sys.exit(main())
