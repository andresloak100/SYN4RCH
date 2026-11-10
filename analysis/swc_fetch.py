"""SWC prototype input: a small public sample of Jiao et al. 2025 single-neuron reconstructions.

    python3 -I analysis/swc_fetch.py [out_dir]        # default data/raw/jiao2025/swc

Two public routes, kept apart because their coordinate frames differ (reports/swc_access.md):

  zenodo/  The 18 Jiao SWCs of Supplementary Fig. 4 inside the Zenodo 14588094 zip (CC BY 4.0),
           read by HTTP range requests (the 918 MB archive is not downloaded). The 18 MouseLight
           files (AA*.swc) in the same folder are NOT fetched. FRAME: mirrored-left (the authors'
           fig_for_plot_neurons.ipynb sets ML' = 11400 - ML for every neuron whose soma ML > 5700).
  portal/  Single files from the public digital-brain.cn portal route (no login, no terms dialog;
           page shows only "Copyright (c) 2023 Digital Brain. All Rights Reserved."). FRAME: true
           hemisphere (Allen CCFv3 um as served). At most 20 new GETs, sequential, 3 s apart,
           User-Agent syn4rch-analysis/0.1. Three portal files already fetched by the user's own probe
           (/private/tmp/synarch_jiao_geometry_probe) are copied byte for byte when their SHA-256
           matches the recorded hash, and are fetched otherwise.

Files are written unchanged. Re-running skips any file already present with the manifest hash.
Prints bytes, full SHA-256 and source per file and writes MANIFEST.tsv."""
import io, os, sys, time, shutil, zipfile, hashlib, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data/raw/jiao2025/swc")
UA = {"User-Agent": "syn4rch-analysis/0.1"}
ZIP_URL = "https://zenodo.org/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip?download=1"
ZIP_DIR = "Supplementary Fig/Supplementary Fig4/swc/"
PORTAL = "https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hy/hy_neuron_download/%s.swc"
PROBE = "/private/tmp/synarch_jiao_geometry_probe"
LIC_Z = "CC BY 4.0 (Jiao et al. 2025, Zenodo 14588094)"
LIC_P = "no licence stated; portal footer 'Copyright (c) 2023 Digital Brain. All Rights Reserved.'; local analysis and private display with attribution only, do not redistribute"
FRAME_Z = "zenodo_mirrored_left"
FRAME_P = "portal_true_hemisphere"

# Already fetched by the user's probe (hash recorded in reports/swc_access.md).
PROBED = {"192710_002": "8873e27f89d858c36483467bfd4dd14dc5a9318072d7148d953484def72f5255",
          "195038_041": "8d6f299313ed1a015329bb21706274dd30ebc6d142eaf224042ef5b7ec3b55d5",
          "201232_056": "8a6bb67e6acb61937f96957204a06f298e9634b200d77e635b0b13781bdd78b3"}
# New portal GETs (<= 20): part B of the access report (every peptide line, both hemispheres where
# they exist) minus 192710_002, plus two Right-hemisphere neurons that also have a Zenodo copy
# (frame check on more than one Right neuron).
NEW = ["202457_066", "202813_035", "202813_009", "202825_072", "202825_065", "212248_071",
       "201221_022", "210105_001", "210104_007", "195039_009", "210073_044", "202170_016",
       "210007_142", "210951_004", "201230_079", "212202_026", "212202_011", "202824_062",
       "211461_026", "202618_143"]
assert len(NEW) <= 20
MAX_BYTES = 100e6
IDX_URL = "https://mouse.digital-brain.cn/projectome/2/srv/info/mouse/hy/mouse.neuron.info.json"
IDX_SHA = "fdf1fa2eff099bc755e5218dc6437e772f12fe0a934c01a4072d846d9b3f4da3"
IDX_CACHE = "/private/tmp/claude-501/-Users-andresarriaga-Documents-SYN4RCH/6fc8191b-130f-4b76-b6a8-69340161f99e/scratchpad/li/hy.neuron.info.json"


class Ranged(io.RawIOBase):
    """Seekable read-only view of a remote file through HTTP Range requests."""
    def __init__(s, url):
        s.pos = 0
        r = urllib.request.urlopen(urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"}), timeout=60)
        s.size = int(r.headers["Content-Range"].split("/")[1]); s.url = r.geturl()
    def seekable(s): return True
    def readable(s): return True
    def tell(s): return s.pos
    def seek(s, off, wh=0):
        s.pos = off if wh == 0 else s.pos + off if wh == 1 else s.size + off; return s.pos
    def read(s, n=-1):
        if n < 0: n = s.size - s.pos
        if n == 0 or s.pos >= s.size: return b""
        end = min(s.size, s.pos + n) - 1
        r = urllib.request.urlopen(urllib.request.Request(s.url, headers={**UA, "Range": "bytes=%d-%d" % (s.pos, end)}), timeout=120)
        b = r.read(); s.pos += len(b); return b
    def readinto(s, b):
        d = s.read(len(b)); b[:len(d)] = d; return len(d)


def sha(b): return hashlib.sha256(b).hexdigest()


def load_manifest(p):
    m = {}
    if os.path.exists(p):
        for i, line in enumerate(open(p)):
            if i == 0: continue
            f = line.rstrip("\n").split("\t")
            m[f[0]] = f
    return m


def main():
    os.makedirs(os.path.join(OUT, "zenodo"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "portal"), exist_ok=True)
    mpath = os.path.join(OUT, "MANIFEST.tsv")
    old = load_manifest(mpath)
    rows = []

    def have(rel):
        p = os.path.join(OUT, rel)
        if rel in old and os.path.exists(p):
            b = open(p, "rb").read()
            if sha(b) == old[rel][3]: return b
        return None

    # --- Zenodo (CC BY 4.0), Jiao files only ---------------------------------------------------
    f = Ranged(ZIP_URL)
    z = zipfile.ZipFile(io.BufferedReader(f, 1 << 16))
    members = [i for i in z.infolist() if i.filename.startswith(ZIP_DIR) and i.filename.endswith(".swc")]
    jiao = [i for i in members if not os.path.basename(i.filename).startswith("AA")]
    ml = [i for i in members if os.path.basename(i.filename).startswith("AA")]
    print("zip %d bytes; %s has %d SWC: %d Jiao, %d MouseLight (MouseLight not fetched)" % (f.size, ZIP_DIR, len(members), len(jiao), len(ml)))
    for i in sorted(jiao, key=lambda i: i.filename):
        base = os.path.basename(i.filename); rel = "zenodo/" + base
        b = have(rel)
        if b is None:
            b = z.read(i); open(os.path.join(OUT, rel), "wb").write(b)
        nid = "_".join(base.split("_")[:2])
        rows.append([rel, nid, str(len(b)), sha(b), ZIP_URL + " :: " + i.filename, FRAME_Z, LIC_Z, "zip member, range request"])
        print("%-48s %9d  %s  <- zip member %s" % (rel, len(b), sha(b), i.filename))

    # --- Portal (no licence; private display only) ----------------------------------------------
    total = sum(int(r[2]) for r in rows)
    first = True
    for nid in list(PROBED) + NEW:
        rel = "portal/%s.swc" % nid; url = PORTAL % nid
        b = have(rel); how = "already present (hash matches manifest)"
        if b is None and nid in PROBED:
            src = os.path.join(PROBE, nid + ".swc")
            if os.path.exists(src) and sha(open(src, "rb").read()) == PROBED[nid]:
                shutil.copyfile(src, os.path.join(OUT, rel)); b = open(os.path.join(OUT, rel), "rb").read()
                how = "copied from user probe %s (fetched 2026-11-10 from the URL)" % src
        if b is None:
            if total > MAX_BYTES: print("STOP: %.1f MB reached" % (total / 1e6)); break
            if not first: time.sleep(3)
            first = False
            r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120)
            b = r.read()
            if r.status != 200 or not b[:2].strip(): raise SystemExit("bad response for " + url)
            open(os.path.join(OUT, rel), "wb").write(b)
            how = "GET %s, Last-Modified %s" % (time.strftime("%Y-%m-%d"), r.headers.get("Last-Modified"))
        else:
            how = old.get(rel, [None] * 8)[7] if how.startswith("already") and rel in old else how
        total += len(b)
        rows.append([rel, nid, str(len(b)), sha(b), url, FRAME_P, LIC_P, how])
        print("%-48s %9d  %s  <- %s" % (rel, len(b), sha(b), url))

    # --- Portal metadata index (public JSON; hemisphere, nucleus, peptide line per neuron) --------
    rel = "portal_index/mouse.neuron.info.json"; url = IDX_URL
    os.makedirs(os.path.join(OUT, "portal_index"), exist_ok=True)
    b = have(rel); how = old[rel][7] if b is not None else None
    if b is None and os.path.exists(IDX_CACHE) and sha(open(IDX_CACHE, "rb").read()) == IDX_SHA:
        shutil.copyfile(IDX_CACHE, os.path.join(OUT, rel)); b = open(IDX_CACHE, "rb").read()
        how = "copied from session cache %s (fetched 2026-11-10 from the URL; SHA-256 matches)" % IDX_CACHE
    if b is None:
        time.sleep(3)
        r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120); b = r.read()
        open(os.path.join(OUT, rel), "wb").write(b); how = "GET %s" % time.strftime("%Y-%m-%d")
    rows.append([rel, "index", str(len(b)), sha(b), url, "metadata", LIC_P, how])
    print("%-48s %9d  %s  <- %s" % (rel, len(b), sha(b), url))

    with open(mpath, "w") as fh:
        fh.write("\t".join(["file", "neuron", "bytes", "sha256", "source", "frame", "licence", "how"]) + "\n")
        for r in rows: fh.write("\t".join(r) + "\n")
    print("files %d, bytes %d (%.1f MB); manifest %s" % (len(rows), total, total / 1e6, mpath))


if __name__ == "__main__":
    main()
