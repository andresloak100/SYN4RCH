"""Target-free local measurements from a small public raw 3D microscopy volume (MouseLight imagery).

Protocol: docs/raw_brain_image_pilot_protocol.json. Descriptive data/feature pilot only: no fitted
biological accuracy, no projection outcome, no axon labels.

* Input: a whole-brain low-resolution octree root TIFF (one channel) from the official public AWS bucket
  janelia-mouselight-imagery (CC BY 4.0), its transform.txt, and the hypothalamus mesh that MouseLight
  back-registered into the same sample space. The mesh only places the crop; it is never a feature.
* Units: transform.txt origin and voxel sizes are nanometres; root-level voxel = s nm, octree level L voxel =
  s / 2**L nm. Physical coordinates are SAMPLE (imaging) space micrometres, not CCF: the CCF displacement field
  (10.3 GB) was not acquired.
* Features are computed only inside a declared crop, optionally with an exclusion mask (outcome pixels); every
  excluded voxel is removed from every statistic. Fluorescent signal is never called an axon or neuron; there
  are no gene profiles.

Usage:
  /opt/anaconda3/bin/python -I analysis/raw_brain_image_features.py pilot   (root-level TIFF; run v1 kept)
  /opt/anaconda3/bin/python -I analysis/raw_brain_image_features.py ktx     (high-resolution PVH KTX texture)
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/raw/raw_brain_image_pilot/mouselight_2018-03-26'
OUT = ROOT / 'data/interim/raw_brain_image_pilot'
PROTOCOL = ROOT / 'docs/raw_brain_image_pilot_protocol.json'
PILOT_RUN_ID = 'raw_image_pilot_mouselight_2018-03-26_v1'
SOURCE_BASE = 'https://janelia-mouselight-imagery.s3.amazonaws.com/'
PINS = {
    'default.0.tif': {'key': 'images/2018-03-26/default.0.tif', 'bytes': 106849192,
                      'sha256': '857b8b4be6c1f0641cb083b077e1b5bb84b1acc29ae8da031d6213bf0bb0f984',
                      'etag': 'fcedf4516a89853d9152905c4494abe0-13'},
    'transform.txt': {'key': 'images/2018-03-26/transform.txt', 'bytes': 112,
                      'sha256': 'bd1b56314f0df56d4b89ff2ecd8ff2dcb27cac6b611176b6b0452003c7d787f5',
                      'etag': 'daa6149bf31e60e18c3f2f5f1ba9aeab'},
    'HY_1097.obj': {'key': 'registration/2018-03-26/HortaObj/HY_1097.obj', 'bytes': 2003618,
                    'sha256': '27f092ec1de239c6be543fab7560821e2d953d520a804f5bc39d9503e5117b0a',
                    'etag': '16709d77ac0f6b3e29bc20733be25110'},
}
CROP_EDGE_UM = 1200.0                     # declared physical edge of the cubic crop
FORBIDDEN_FEATURE_WORDS = ('axon', 'dendrite', 'neuron', 'gene', 'projection', 'target', 'synapse')


class Refusal(ValueError):
    pass


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------- transform, units, coordinates

def parse_transform(text):
    vals = {}
    for line in text.strip().splitlines():
        k, v = line.split(':')
        vals[k.strip()] = float(v)
    need = {'ox', 'oy', 'oz', 'sx', 'sy', 'sz', 'nl'}
    if set(vals) != need or vals['nl'] != int(vals['nl']) or min(vals['sx'], vals['sy'], vals['sz']) <= 0:
        raise Refusal(f'transform.txt fields {sorted(vals)} invalid')
    vals['nl'] = int(vals['nl'])
    return vals


def voxel_size_um(transform, level=0):
    """Voxel size (x, y, z) in micrometres at octree level `level` (0 = root, coarsest)."""
    if not 0 <= level < transform['nl']:
        raise Refusal('octree level outside 0..nl-1')
    f = 2.0 ** level
    return np.array([transform['sx'], transform['sy'], transform['sz']]) / 1000.0 / f


def index_to_sample_um(zyx_index, transform, level=0):
    """Voxel index (z, y, x) -> sample-space (x, y, z) micrometres, Horta corner convention (not CCF)."""
    z, y, x = np.asarray(zyx_index, float).T
    s = voxel_size_um(transform, level)
    o = np.array([transform['ox'], transform['oy'], transform['oz']]) / 1000.0
    return np.stack([o[0] + x * s[0], o[1] + y * s[1], o[2] + z * s[2]], axis=-1)


def sample_um_to_index(xyz_um, transform, level=0):
    s = voxel_size_um(transform, level)
    o = np.array([transform['ox'], transform['oy'], transform['oz']]) / 1000.0
    x, y, z = ((np.asarray(xyz_um, float) - o) / s).T
    return np.stack([z, y, x], axis=-1)


def read_obj_vertices(path):
    v = [list(map(float, line.split()[1:4])) for line in Path(path).read_text().splitlines() if line.startswith('v ')]
    if not v:
        raise Refusal('mesh has no vertices')
    return np.asarray(v)


# ---------------------------------------------------------------- volume and crop

def load_tiff(path):
    import tifffile
    with tifffile.TiffFile(path) as t:
        if len(t.series) != 1:
            raise Refusal('expected one image series')
        s = t.series[0]
        meta = {'shape_zyx': list(s.shape), 'dtype': str(s.dtype), 'axes': s.axes, 'pages': len(t.pages),
                'compression': int(t.pages[0].compression)}
        vol = s.asarray()
    if vol.ndim != 3:
        raise Refusal('expected a 3D single-channel volume')
    return vol, meta


def crop_bounds(center_zyx, edge_um, spacing_zyx_um, shape_zyx):
    """Integer [start, stop) per axis for a cube of `edge_um` around a centre; refuses if it leaves the volume."""
    half = np.floor(edge_um / np.asarray(spacing_zyx_um) / 2.0 + 0.5).astype(int)      # explicit half-up rounding
    c = np.floor(np.asarray(center_zyx, float) + 0.5).astype(int)
    lo, hi = c - half, c + half
    if (lo < 0).any() or (hi > np.asarray(shape_zyx)).any() or (half < 1).any():
        raise Refusal(f'crop {lo.tolist()}..{hi.tolist()} leaves the volume {list(shape_zyx)}')
    return [(int(a), int(b)) for a, b in zip(lo, hi)]


def take(vol, bounds):
    (z0, z1), (y0, y1), (x0, x1) = bounds
    return vol[z0:z1, y0:y1, x0:x1]


# ---------------------------------------------------------------- target-free features

def _valid(crop, exclude):
    if exclude is None:
        return np.ones(crop.shape, bool)
    exclude = np.asarray(exclude, bool)
    if exclude.shape != crop.shape:
        raise Refusal('exclusion mask shape differs from the crop')
    return ~exclude


def intensity_qc(crop, exclude=None, saturation_value=None, source_units='uncalibrated fluorescence as stored'):
    """Intensity statistics in the SOURCE units and dtype actually present (no unit conversion)."""
    v = crop[_valid(crop, exclude)].astype(float)
    if v.size == 0:
        return {'status': 'NO_VALID_VOXELS', 'dtype': str(crop.dtype), 'units': source_units}
    sat = saturation_value if saturation_value is not None else float(np.iinfo(crop.dtype).max) \
        if np.issubdtype(crop.dtype, np.integer) else float(v.max())
    p = np.percentile(v, [1, 10, 50, 90, 99, 99.9])
    return {'status': 'OK', 'n_voxels': int(v.size), 'mean': float(v.mean()), 'std': float(v.std()),
            'p01': float(p[0]), 'p10': float(p[1]), 'median': float(p[2]), 'p90': float(p[3]),
            'p99': float(p[4]), 'p999': float(p[5]), 'zero_fraction': float((v == 0).mean()),
            'saturation_value': float(sat), 'saturated_fraction': float((v >= sat).mean()),
            'dynamic_range_p99_over_p10': float(p[4] / p[1]) if p[1] > 0 else None,
            'dtype': str(crop.dtype), 'units': source_units}


def foreground_sensitivity(crop, exclude=None, factors=(0.8, 1.0, 1.25)):
    """Bright-signal fraction at Otsu x factors: the spread is the segmentation uncertainty."""
    from skimage.filters import threshold_otsu
    v = crop[_valid(crop, exclude)].astype(float)
    if v.size == 0:
        return {'status': 'NO_VALID_VOXELS'}
    if np.ptp(v) == 0:
        return {'status': 'NO_SIGNAL', 'otsu': None, 'fractions': None, 'note': 'constant intensity'}
    t = float(threshold_otsu(v))
    fr = {str(f): float((v > t * f).mean()) for f in factors}
    return {'status': 'OK', 'otsu': t, 'fractions': fr, 'fraction_range': float(max(fr.values()) - min(fr.values())),
            'note': 'bright-signal fraction; NOT a cell, fibre or axon count'}


TUBE_GAP_MIN = 0.1      # sum(l2 - l1) / sum(l3) below this: no tube direction (fixed on synthetic checks only)


def _undefined_orientation(status, sigma_um, rho_um, n):
    return {'status': status, 'dominant_direction_sample_xyz': None, 'order_parameter': None, 'mean_tube_coherence': None,
            'sigma_um': float(sigma_um), 'rho_um': float(rho_um), 'n_voxels_used': int(n),
            'note': 'orientation undefined: no usable intensity structure'}


def structure_tensor_orientation(crop, spacing_zyx_um, sigma_um, rho_um, exclude=None, max_voxels=200000):
    """Local orientation from the 3D structure tensor in PHYSICAL units (anisotropic spacing handled).

    The tube-like direction is the eigenvector of the smallest eigenvalue. Returns the dominant direction (principal
    axis of the strength-weighted orientation tensor, sample-space x, y, z), its order parameter and coherence.
    A flat or fully excluded crop returns an undefined direction and order (status NO_SIGNAL / NO_VALID_VOXELS)."""
    from scipy.ndimage import gaussian_filter
    sp = np.asarray(spacing_zyx_um, float)
    valid = _valid(crop, exclude)
    if not valid.any():
        return _undefined_orientation('NO_VALID_VOXELS', sigma_um, rho_um, 0)
    base = crop.astype(float)
    if exclude is not None:                                    # excluded (outcome) voxels cannot leak via smoothing
        base = np.where(valid, base, np.median(base[valid]))
    if np.ptp(base[valid]) == 0:
        return _undefined_orientation('NO_SIGNAL', sigma_um, rho_um, int(valid.sum()))
    img = gaussian_filter(base, sigma=sigma_um / sp)
    gz, gy, gx = np.gradient(img, *sp)
    g = {'x': gx, 'y': gy, 'z': gz}
    order = ['x', 'y', 'z']
    comps = {(a, b): gaussian_filter(g[a] * g[b], sigma=rho_um / sp)
             for i, a in enumerate(order) for b in order[i:]}
    idx = np.flatnonzero(valid.ravel())
    if idx.size > max_voxels:
        idx = idx[np.linspace(0, idx.size - 1, max_voxels).astype(int)]     # deterministic subsample
    T = np.empty((idx.size, 3, 3))
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            T[:, i, j] = comps[(a, b) if (a, b) in comps else (b, a)].ravel()[idx]
    w, V = np.linalg.eigh(T)                                   # ascending eigenvalues l1 <= l2 <= l3
    lam = np.clip(w, 0, None)
    scale = float(np.max(lam)) if lam.size else 0.0
    if scale <= 0 or (lam[:, 2] - lam[:, 0]).max() <= 1e-12 * scale:
        return _undefined_orientation('NO_SIGNAL', sigma_um, rho_um, int(idx.size))
    # The tube direction (eigenvector of l1) is only defined when l1 is separated from l2. For a planar sheet
    # l1 ~ l2 and the smallest eigenvector is arbitrary within the sheet, so confidence uses the gap l2 - l1.
    tube_gap = lam[:, 1] - lam[:, 0]
    tube_ratio = float(tube_gap.sum() / max(lam[:, 2].sum(), 1e-300))
    planar_ratio = float((lam[:, 2] - lam[:, 1]).sum() / max(lam[:, 2].sum(), 1e-300))
    if tube_ratio < TUBE_GAP_MIN:
        out = _undefined_orientation('NO_TUBE_DIRECTION', sigma_um, rho_um, int(idx.size))
        out.update(tube_gap_ratio=tube_ratio, planar_ratio=planar_ratio,
                   note='no separated smallest eigenvalue (planar or isotropic texture): tube direction undefined')
        return out
    tube_dir = V[:, :, 0]
    coherence = np.where(lam[:, 1] + lam[:, 0] > 0, tube_gap / np.maximum(lam[:, 1] + lam[:, 0], 1e-300), 0.0)
    wts = tube_gap / tube_gap.sum()
    Q = np.einsum('n,ni,nj->ij', wts, tube_dir, tube_dir)
    qw, qv = np.linalg.eigh(Q)
    d = qv[:, 2] * np.sign(qv[np.argmax(np.abs(qv[:, 2])), 2])
    return {'status': 'OK', 'dominant_direction_sample_xyz': d.tolist(), 'order_parameter': float(1.5 * qw[2] - 0.5),
            'mean_tube_coherence': float(np.average(coherence, weights=wts)),
            'tube_gap_ratio': tube_ratio, 'planar_ratio': planar_ratio,
            'sigma_um': float(sigma_um), 'rho_um': float(rho_um), 'n_voxels_used': int(idx.size),
            'note': 'tube-like texture direction weighted by (l2 - l1); fluorescent signal only, not traced fibres or axons'}


def physical_log(img, spacing_zyx_um, sigma_um):
    """Laplacian of Gaussian in physical units: sum over axes of d2/dx_a2 with each axis divided by its spacing^2."""
    from scipy.ndimage import gaussian_filter
    sp = np.asarray(spacing_zyx_um, float)
    sv = sigma_um / sp
    out = np.zeros(img.shape)
    for ax in range(3):
        order = [0, 0, 0]; order[ax] = 2
        out += gaussian_filter(img, sigma=sv, order=order) / sp[ax] ** 2
    return out


def blob_scale(crop, spacing_zyx_um, radii_um, exclude=None):
    """Scale-normalised physical LoG at the brightest valid location.

    A Gaussian blob of standard deviation s (radius r = s * sqrt(3)) maximises sigma^2 * |LoG| at
    sigma = s * sqrt(2/3); each candidate radius r is therefore probed at sigma = r * sqrt(2) / 3."""
    valid = _valid(crop, exclude)
    if not valid.any():
        return {'status': 'NO_VALID_VOXELS'}
    img = crop.astype(float)
    if exclude is not None:
        img = np.where(valid, img, np.median(img[valid]))
    if np.ptp(img[valid]) == 0:
        return {'status': 'NO_SIGNAL'}
    peak = np.unravel_index(np.argmax(np.where(valid, img, -np.inf)), img.shape)
    resp = []
    for r in radii_um:
        s = r * np.sqrt(2) / 3
        resp.append(float(-s ** 2 * physical_log(img, spacing_zyx_um, s)[peak]))
    best = int(np.argmax(resp))
    return {'status': 'OK', 'radius_um': float(radii_um[best]), 'radii_probed_um': [float(r) for r in radii_um],
            'normalised_responses': resp, 'peak_index_zyx': [int(p) for p in peak],
            'note': 'scale of the brightest compact intensity object; NOT a soma or cell identification'}


def support_flags(spacing_zyx_um):
    """What sampling can and cannot support (from voxel size alone, never from outcomes).

    Sampling can show that somata-scale objects are representable; it cannot certify optical resolution, so
    single-neurite resolution is never certified by voxel size alone."""
    worst = float(np.max(spacing_zyx_um))
    return {'voxel_size_um_zyx': [float(s) for s in spacing_zyx_um],
            'soma_scale_sampling_adequate': worst <= 3.0,
            'single_neurite_resolution': 'NOT_CERTIFIED: sampling alone cannot establish optical resolution '
                                         '(point-spread function unknown)',
            'regional_texture_orientation_supported': True,
            'reason': f'largest voxel edge {worst:.2f} um; somata (~10-20 um) need several voxels per axis'}


def extract_features(crop, spacing_zyx_um, exclude=None, source_units='uncalibrated fluorescence as stored',
                     orientation_sigma_um=None):
    flags = support_flags(spacing_zyx_um)
    edge = float(np.max(spacing_zyx_um))
    sig = orientation_sigma_um if orientation_sigma_um is not None else edge
    feats = {'support': flags, 'intensity_qc': intensity_qc(crop, exclude, source_units=source_units),
             'bright_signal_fraction': foreground_sensitivity(crop, exclude),
             'orientation': structure_tensor_orientation(crop, spacing_zyx_um, sigma_um=sig, rho_um=3 * sig,
                                                         exclude=exclude)}
    if flags['soma_scale_sampling_adequate']:
        feats['compact_object_scale'] = blob_scale(crop, spacing_zyx_um, [2.0, 3.0, 4.0, 6.0, 8.0, 12.0], exclude)
    else:
        feats['compact_object_scale'] = {'status': 'NOT_SUPPORTED_AT_THIS_RESOLUTION'}
    check_names(feats)
    return feats


def check_names(obj, path=''):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if any(w in k.lower() for w in FORBIDDEN_FEATURE_WORDS):
                raise Refusal(f'forbidden feature name {path + k}')
            check_names(v, path + k + '.')


# ---------------------------------------------------------------- the descriptive pilot

def verify_pins(raw=RAW):
    for name, pin in PINS.items():
        p = Path(raw) / name
        if not p.is_file() or p.stat().st_size != pin['bytes'] or sha256(p) != pin['sha256']:
            raise Refusal(f'{name} missing or differs from its pin')


def pilot(raw=RAW, out_root=OUT, run_id=PILOT_RUN_ID):
    verify_pins(raw)
    run_dir = Path(out_root) / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    started = {'run_id': run_id, 'started_utc': datetime.now(timezone.utc).isoformat(), 'source_sha256': sha256(__file__),
               'protocol_sha256': sha256(PROTOCOL) if PROTOCOL.exists() else None,
               'inputs': {n: dict(p, url=SOURCE_BASE + p['key']) for n, p in PINS.items()},
               'classification': 'DESCRIPTIVE_DATA_FEATURE_PILOT_NOT_A_BENCHMARK',
               'environment': {'python': platform.python_version(), 'numpy': np.__version__}}
    (run_dir / 'STARTED.json').write_text(json.dumps(started, indent=2, sort_keys=True))
    tr = parse_transform((Path(raw) / 'transform.txt').read_text())
    vol, meta = load_tiff(Path(raw) / 'default.0.tif')
    spacing_xyz = voxel_size_um(tr, 0)
    spacing_zyx = spacing_xyz[::-1]
    hy = read_obj_vertices(Path(raw) / 'HY_1097.obj')
    centre_xyz = hy.mean(axis=0)                               # HY mesh centroid, sample space (placement only)
    centre_zyx = sample_um_to_index(centre_xyz, tr, 0)
    bounds = crop_bounds(centre_zyx, CROP_EDGE_UM, spacing_zyx, vol.shape)
    crop = take(vol, bounds)
    np.save(run_dir / 'crop_uint16.npy', crop)
    lo = index_to_sample_um([b[0] for b in bounds], tr, 0); hi = index_to_sample_um([b[1] for b in bounds], tr, 0)
    record = {
        'specimen': {'mouselight_image_date': '2018-03-26', 'finished_tracings_listed': False,
                     'in_hypothalamic_test_cohort': False,
                     'note': 'imaged brain with no finished or completed tracings in the public bucket listing'},
        'channel': {'file': 'default.0.tif', 'channel_index': 0,
                    'fluorophore': 'NOT RECORDED in the bucket; channel identity unverified'},
        'volume': meta, 'octree_levels': tr['nl'], 'level_used': 0,
        'voxel_size_um_xyz': spacing_xyz.tolist(),
        'transform_nm': {k: tr[k] for k in ('ox', 'oy', 'oz', 'sx', 'sy', 'sz')},
        'coordinate_frame': 'MouseLight sample (imaging) space, micrometres, Horta corner convention; NOT CCF',
        'crop': {'bounds_zyx': bounds, 'shape_zyx': list(crop.shape), 'edge_um_declared': CROP_EDGE_UM,
                 'sample_um_min_xyz': lo.tolist(), 'sample_um_max_xyz': hi.tolist(),
                 'placement': 'centred on the HY_1097 mesh centroid (MouseLight back-registered CCF mesh in sample '
                              'space); the mesh places the crop and is not a feature',
                 'hy_mesh_bbox_um_xyz': [hy.min(0).tolist(), hy.max(0).tolist()]},
        'features': extract_features(crop, spacing_zyx, source_units='raw uint16 fluorescence as stored in default.0.tif (uncalibrated)'),
        'gaps': ['CCF coordinates of crop voxels need Transform.2018-03-26.h5 (10.3 GB), not acquired',
                 'channel fluorophore identity is not stated in the bucket',
                 'root level (19 x 19 x 62 um) cannot resolve somata or single neurites',
                 'Horta voxel corner-versus-centre convention unverified (half-voxel uncertainty)'],
        'not_claimed': ['axons, neurons or projections', 'gene profiles', 'biological accuracy of any model']}
    (run_dir / 'features.json').write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False))
    files = {p.name: sha256(p) for p in sorted(run_dir.iterdir()) if p.is_file()}
    (run_dir / 'RUN_COMPLETE.json').write_text(json.dumps({'run_id': run_id, 'files_sha256': files,
                                                           'finished_utc': datetime.now(timezone.utc).isoformat()},
                                                          indent=2, sort_keys=True))
    return run_dir, record


# ---------------------------------------------------------------- KTX1 decoder (MouseLight Horta textures)

KTX_IDENTIFIER = b'\xabKTX 11\xbb\r\n\x1a\n'
GL_UNSIGNED_BYTE, GL_RG, GL_RG8 = 0x1401, 0x8227, 0x822B
KTX_ROUTE = ROOT / 'data/interim/paired_projection_search_round6/image_routes'
KTX_PIN = {'file': 'source_pilot_PVH_20170419_273377.ktx', 'bytes': 11058544,
           'sha256': '08eee2363d455e832da48275dfcd11ddb5128fd1deb4a0e422aef94271930627',
           'url': SOURCE_BASE + 'images/2017-04-19/ktx/2/7/3/3/7/7/block_8_xy_273377.ktx'}
KTX_RUN_ID = 'raw_image_ktx_pvh_20170419_273377_v3'   # v1/v2 kept; v3 uses the (l2 - l1) tube-direction rule


def read_ktx1(path):
    """Strict KTX1 reader for uncompressed little-endian RG8 3D textures; returns the first mip (Z, Y, X, C)."""
    import struct
    b = Path(path).read_bytes()
    if len(b) < 64 or b[:12] != KTX_IDENTIFIER:
        raise Refusal('not a KTX1 file')
    h = struct.unpack('<13I', b[12:64])
    names = ['endianness', 'glType', 'glTypeSize', 'glFormat', 'glInternalFormat', 'glBaseInternalFormat',
             'pixelWidth', 'pixelHeight', 'pixelDepth', 'numberOfArrayElements', 'numberOfFaces',
             'numberOfMipmapLevels', 'bytesOfKeyValueData']
    hd = dict(zip(names, h))
    if hd['endianness'] != 0x04030201:
        raise Refusal('KTX endianness is not little-endian 0x04030201')
    if (hd['glType'], hd['glTypeSize'], hd['glFormat'], hd['glInternalFormat'], hd['glBaseInternalFormat']) != \
            (GL_UNSIGNED_BYTE, 1, GL_RG, GL_RG8, GL_RG):
        raise Refusal('KTX is not uncompressed unsigned-byte RG8')
    if hd['numberOfArrayElements'] != 0 or hd['numberOfFaces'] not in (0, 1) or hd['numberOfMipmapLevels'] < 1:
        raise Refusal('KTX array/face/mip layout unsupported')
    w, ht, d = hd['pixelWidth'], hd['pixelHeight'], hd['pixelDepth']
    if min(w, ht, d) < 1:
        raise Refusal('KTX is not a 3D texture')
    kv, i, end = {}, 64, 64 + hd['bytesOfKeyValueData']
    if end > len(b):
        raise Refusal('key-value block runs past the file end')
    while i < end:
        n = struct.unpack('<I', b[i:i + 4])[0]
        raw = b[i + 4:i + 4 + n]
        k, sep, v = raw.partition(b'\0')
        if not sep:
            raise Refusal('malformed key-value pair')
        kv[k.decode()] = v.rstrip(b'\0').decode()
        i += 4 + n + (-n % 4)
    if i != end:
        raise Refusal('key-value block length mismatch')
    row = w * 2
    if row % 4:
        raise Refusal('row padding not supported for this texture width')
    size = struct.unpack('<I', b[end:end + 4])[0]
    if size != row * ht * d or end + 4 + size > len(b):
        raise Refusal(f'first mip imageSize {size} != {row * ht * d} or truncated')
    vol = np.frombuffer(b, dtype=np.uint8, count=size, offset=end + 4).reshape(d, ht, w, 2).copy()
    return vol, hd, kv


def ktx_geometry(kv, shape_zyx):
    """Voxel spacing and origin (sample-space micrometres) from the embedded texcoord transform; cross-checked."""
    import ast
    if kv.get('distance_units') != 'micrometers':
        raise Refusal('KTX distance units are not micrometers')
    M = np.array(ast.literal_eval(' '.join(kv['xyz_from_texcoord_xform'].split()).replace(' ', ',')
                                  .replace('[,', '[').replace(',]', ']')), float)
    if M.shape != (4, 4) or not np.allclose(M[:3, :3], np.diag(np.diag(M[:3, :3]))):
        raise Refusal('texcoord transform is not axis-aligned 4x4')
    extent = np.diag(M[:3, :3]); origin = M[:3, 3]
    d, h, w = shape_zyx
    spacing_xyz = extent / np.array([w, h, d])
    corners = np.array(ast.literal_eval(kv['corner_xyzs']), float)
    if not (np.allclose(corners.min(0), origin, atol=1e-3) and np.allclose(corners.max(0), origin + extent, atol=1e-3)):
        raise Refusal('corner_xyzs disagree with the texcoord transform')
    return {'spacing_um_xyz': spacing_xyz.tolist(), 'origin_um_xyz': origin.tolist(), 'extent_um_xyz': extent.tolist(),
            'voxel_centre_convention': 'texcoord (i + 0.5) / N, corners at 0 and 1 (from corner_texcoords)'}


def projection_figure(vol_zyxc, spacing_xyz, feats, path, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d, h, w, c = vol_zyxc.shape
    sx, sy, sz = spacing_xyz
    fig, axes = plt.subplots(c, 2, figsize=(11, 4.2 * c), gridspec_kw={'width_ratios': [w * sx, w * sx]})
    for ch in range(c):
        v = vol_zyxc[..., ch].astype(float)
        for k, (img, ext, lab) in enumerate(((v.max(0), [0, w * sx, h * sy, 0], 'y (um)'),
                                             (v.max(1), [0, w * sx, d * sz, 0], 'z (um)'))):
            ax = axes[ch, k]
            lo, hi = np.percentile(img, [2, 99.7])               # display contrast from the projection itself
            ax.imshow(img, cmap='gray', vmin=lo, vmax=hi, extent=ext, interpolation='nearest')
            ax.set_xlabel('x (um)'); ax.set_ylabel(lab)
            ax.plot([5, 25], [ext[2] - 5, ext[2] - 5], color='yellow', lw=3)
            ax.text(5, ext[2] - 9, '20 um', color='yellow', fontsize=9)
            ax.set_title(f'channel {ch}: max projection along {"z" if k == 0 else "y"}', fontsize=10)
        f = feats[f'channel_{ch}']
        o = f['orientation']
        txt = (f"median {f['intensity_qc']['median']:.0f}, p99 {f['intensity_qc']['p99']:.0f} (8-bit rescaled); "
               f"bright fraction range {f['bright_signal_fraction'].get('fraction_range', float('nan')):.2f}; "
               f"orientation order {o['order_parameter'] if o['order_parameter'] is None else round(o['order_parameter'], 2)}; "
               f"compact-object radius {f['compact_object_scale'].get('radius_um')} um")
        axes[ch, 0].text(0, -0.24, txt, transform=axes[ch, 0].transAxes, fontsize=8, va='top')
    fig.suptitle(title + '\nRescaled 8-bit fluorescence; channel identities unknown; sample frame, not CCF; '
                 'no neuron, gene or axon identification', fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def ktx_pilot(out_root=OUT, run_id=KTX_RUN_ID, path=None):
    path = Path(path or KTX_ROUTE / KTX_PIN['file'])
    if path.stat().st_size != KTX_PIN['bytes'] or sha256(path) != KTX_PIN['sha256']:
        raise Refusal('KTX source differs from its pin')
    run_dir = Path(out_root) / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    receipts = {n: sha256(KTX_ROUTE / n) for n in ('source_pilot_receipt.json', 'source_pilot_quality.json')
                if (KTX_ROUTE / n).exists()}
    (run_dir / 'STARTED.json').write_text(json.dumps({
        'run_id': run_id, 'started_utc': datetime.now(timezone.utc).isoformat(), 'source_sha256': sha256(__file__),
        'protocol_sha256': sha256(PROTOCOL) if PROTOCOL.exists() else None, 'input': dict(KTX_PIN),
        'route_receipts_sha256': receipts, 'classification': 'DESCRIPTIVE_DATA_FEATURE_PILOT_NOT_A_BENCHMARK'},
        indent=2, sort_keys=True))
    vol, hd, kv = read_ktx1(path)
    geo = ktx_geometry(kv, vol.shape[:3])
    sp_zyx = list(reversed(geo['spacing_um_xyz']))
    units = {ch: (f'8-bit rescaled fluorescence (gamma {kv.get(f"channel_{ch}_intensity_gamma")}, '
                  f'scale {kv.get(f"channel_{ch}_intensity_scale")}, offset {kv.get(f"channel_{ch}_intensity_offset")}); '
                  'not linear camera counts') for ch in range(vol.shape[3])}
    feats = {f'channel_{ch}': extract_features(vol[..., ch], sp_zyx, source_units=units[ch],
                                               orientation_sigma_um=1.0) for ch in range(vol.shape[3])}
    record = {
        'specimen': {'mouselight_specimen_id': kv.get('specimen_id'), 'image_date': '2017-04-19',
                     'octree_path': kv.get('octree_path'), 'multiscale_level': kv.get('multiscale_level_id'),
                     'region_selection': 'PVH atlas-mesh bounding-box centre (ROOT receipt); patch not verified to be '
                                         'soma-centred', 'cohort_note': 'brain of the MouseLight hypothalamic cohort; this is '
                                         'raw imagery only, no trace geometry opened'},
        'header': hd, 'ktx_metadata': kv, 'geometry': geo, 'shape_zyxc': list(vol.shape), 'dtype': str(vol.dtype),
        'channels': {'count': int(vol.shape[3]), 'identity': 'UNKNOWN (not stated in metadata)'},
        'coordinate_frame': 'MouseLight native sample space, micrometres; NOT CCF',
        'features': feats,
        'not_claimed': ['neurons, somata, axons or genes', 'biological accuracy', 'CCF location']}
    (run_dir / 'features.json').write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False))
    projection_figure(vol, geo['spacing_um_xyz'], feats, run_dir / 'projections.png',
                      f"MouseLight {kv.get('specimen_id')} PVH-region texture {kv.get('octree_path')}")
    files = {p.name: sha256(p) for p in sorted(run_dir.iterdir()) if p.is_file()}
    (run_dir / 'RUN_COMPLETE.json').write_text(json.dumps({'run_id': run_id, 'files_sha256': files,
                                                           'finished_utc': datetime.now(timezone.utc).isoformat()},
                                                          indent=2, sort_keys=True))
    return run_dir, record


if __name__ == '__main__':
    if sys.argv[1:] == ['pilot']:
        d, rec = pilot()
        print('pilot complete', d, json.dumps(rec['features']['support']))
    elif sys.argv[1:] == ['ktx']:
        d, rec = ktx_pilot()
        print('ktx pilot complete', d)
    else:
        sys.exit(__doc__)
