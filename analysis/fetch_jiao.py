"""D-109 input: the small region-level tables of Jiao et al. 2025 (Nat Neurosci, doi:10.1038/s41593-025-01919-0).

    python3 analysis/fetch_jiao.py [out_dir]      # default data/raw/jiao2025

Selected members of the Zenodo 14588094 zip (CC BY 4.0), read by HTTP range requests so the
918 MB archive is not downloaded, plus Nature Source Data Fig. 2 (no licence stated; identical
bytes to the Zenodo all_neuron_soma_infor.xlsx). Prints each file's size and full SHA-256.
Sources, member paths and the audit: reports/d109_jiao_audit.md."""
import io, os, sys, zipfile, hashlib, urllib.request
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data/raw/jiao2025")
os.makedirs(OUT, exist_ok=True)
URL = "https://zenodo.org/records/14588094/files/Code%20for%20Hypothalamic%20projectome.zip?download=1"
UA = {"User-Agent": "syn4rch-analysis/0.1"}
class Ranged(io.RawIOBase):
    def __init__(s, url):
        s.url = url; s.pos = 0
        r = urllib.request.urlopen(urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"}))
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
        r = urllib.request.urlopen(urllib.request.Request(s.url, headers={**UA, "Range": "bytes=%d-%d" % (s.pos, end)}))
        b = r.read(); s.pos += len(b); return b
    def readinto(s, b):
        d = s.read(len(b)); b[:len(d)] = d; return len(d)
WANT = ("all_neurons_projection_length_fig1.csv", "all_neurons_terminal_terminal_number_Supplementary_fig3_sort.csv",
        "All_soma_infor_HY_neurons.csv", "geno_location_labels.csv", "all_neuron_soma_infor.xlsx")
f = Ranged(URL); print("zip size", f.size)
z = zipfile.ZipFile(io.BufferedReader(f, 1 << 20))
done = set()
for info in z.infolist():
    base = os.path.basename(info.filename)
    if base in WANT and base not in done and info.file_size < 20e6:
        data = z.read(info); done.add(base)
        open(os.path.join(OUT, base), "wb").write(data)
        print("%-70s %9d  %s  <- %s" % (base, len(data), hashlib.sha256(data).hexdigest(), info.filename))
print("missing:", sorted(set(WANT) - done))
nat = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-01919-0/MediaObjects/41593_2025_1919_MOESM5_ESM.xlsx"
try:
    d = urllib.request.urlopen(urllib.request.Request(nat, headers=UA)).read()
    open(os.path.join(OUT, "nature_source_data_fig2.xlsx"), "wb").write(d)
    print("nature_source_data_fig2.xlsx", len(d), hashlib.sha256(d).hexdigest())
except Exception as e:
    print("Nature file failed:", e)
