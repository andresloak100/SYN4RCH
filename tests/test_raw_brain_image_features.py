import json
import shutil

import numpy as np
import pytest

from analysis import raw_brain_image_features as F

TR = F.parse_transform('ox: 1000000\noy: 2000000\noz: 3000000\nsx: 16000\nsy: 16000\nsz: 64000\nnl: 5\n')


def line_volume(direction_xyz, spacing_zyx, shape=(40, 64, 64), width_um=3.0):
    """Bright straight tube through the centre along a PHYSICAL direction (x, y, z)."""
    d = np.asarray(direction_xyz, float); d /= np.linalg.norm(d)
    z, y, x = np.indices(shape).astype(float)
    pz, py, px = z * spacing_zyx[0], y * spacing_zyx[1], x * spacing_zyx[2]
    c = np.array([px.mean(), py.mean(), pz.mean()])
    rel = np.stack([px - c[0], py - c[1], pz - c[2]], -1)
    dist = np.linalg.norm(rel - (rel @ d)[..., None] * d, axis=-1)
    return (1000 * np.exp(-(dist / width_um) ** 2) + 10).astype(np.uint16)


def blob_volume(radius_um, spacing_zyx, shape=(48, 48, 48)):
    z, y, x = np.indices(shape).astype(float)
    sp = np.asarray(spacing_zyx)
    c = (np.asarray(shape) - 1) / 2
    r2 = ((z - c[0]) * sp[0]) ** 2 + ((y - c[1]) * sp[1]) ** 2 + ((x - c[2]) * sp[2]) ** 2
    sigma = radius_um / np.sqrt(3)
    return (2000 * np.exp(-r2 / (2 * sigma ** 2)) + 5).astype(np.uint16)


# ---------------- units and coordinates

def test_transform_units_and_levels():
    assert np.allclose(F.voxel_size_um(TR, 0), [16, 16, 64])
    assert np.allclose(F.voxel_size_um(TR, 4), [1, 1, 4])
    with pytest.raises(F.Refusal):
        F.voxel_size_um(TR, 5)
    with pytest.raises(F.Refusal):
        F.parse_transform('ox: 1\noy: 2\n')


def test_coordinate_roundtrip_exact():
    idx = np.array([[0, 0, 0], [2, 3, 4], [10.5, 1, 7]])
    xyz = F.index_to_sample_um(idx, TR, 0)
    assert np.allclose(xyz[0], [1000, 2000, 3000]) and np.allclose(xyz[1], [1000 + 64, 2000 + 48, 3000 + 128])
    assert np.allclose(F.sample_um_to_index(xyz, TR, 0), idx)


def test_crop_bounds_physical_and_refusal():
    b = F.crop_bounds([20, 32, 32], 320.0, [64, 16, 16], (40, 64, 64))
    assert [hi - lo for lo, hi in b] == [6, 20, 20]               # 320 um / spacing, rounded halves
    with pytest.raises(F.Refusal):
        F.crop_bounds([1, 32, 32], 320.0, [64, 16, 16], (40, 64, 64))


# ---------------- intensity and segmentation uncertainty

def test_intensity_qc_exact_and_exclusion():
    a = np.array([[[0, 10, 20, 65535]]], dtype=np.uint16)
    q = F.intensity_qc(a, source_units='test units')
    assert q['zero_fraction'] == 0.25 and q['saturated_fraction'] == 0.25 and q['mean'] == pytest.approx(65565 / 4)
    assert q['dtype'] == 'uint16' and q['units'] == 'test units'
    ex = np.array([[[False, False, False, True]]])
    q2 = F.intensity_qc(a, exclude=ex)
    assert q2['n_voxels'] == 3 and q2['saturated_fraction'] == 0.0
    assert F.intensity_qc(a, exclude=np.ones_like(ex))['status'] == 'NO_VALID_VOXELS'
    assert F.intensity_qc(np.array([[[1, 2]]], np.uint8))['saturation_value'] == 255


def test_foreground_fraction_monotone_in_threshold():
    rng = np.random.default_rng(0)
    v = rng.gamma(2, 100, (20, 20, 20)).astype(np.uint16)
    fr = F.foreground_sensitivity(v)['fractions']
    assert fr['0.8'] >= fr['1.0'] >= fr['1.25']


# ---------------- orientation in physical units

@pytest.mark.parametrize('direction', [(1, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 1)])
@pytest.mark.parametrize('spacing', [(1.0, 1.0, 1.0), (3.0, 1.0, 1.0)])
def test_planted_tube_direction_recovered(direction, spacing):
    vol = line_volume(direction, spacing)
    o = F.structure_tensor_orientation(vol, spacing, sigma_um=2.0, rho_um=4.0)
    d = np.asarray(direction, float) / np.linalg.norm(direction)
    assert abs(np.dot(o['dominant_direction_sample_xyz'], d)) > 0.95
    assert o['order_parameter'] > 0.5


def test_isotropic_noise_has_low_order():
    rng = np.random.default_rng(1)
    vol = rng.poisson(100, (32, 32, 32)).astype(np.uint16)
    o = F.structure_tensor_orientation(vol, (1, 1, 1), sigma_um=1.0, rho_um=3.0)
    assert o['order_parameter'] < 0.3


def test_excluded_outcome_voxels_cannot_change_features():
    vol = line_volume((1, 0, 0), (1, 1, 1))
    ex = np.zeros(vol.shape, bool); ex[:, :8, :] = True            # e.g. a projection-region slab
    f1 = F.extract_features(vol, (1, 1, 1), exclude=ex)
    v2 = vol.copy(); v2[ex] = 60000                                # change only excluded pixels
    f2 = F.extract_features(v2, (1, 1, 1), exclude=ex)
    assert json.dumps(f1, sort_keys=True) == json.dumps(f2, sort_keys=True)


# ---------------- scale

@pytest.mark.parametrize('spacing', [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (1.0, 0.5, 0.5)])
@pytest.mark.parametrize('radius', [4.0, 8.0])
def test_planted_blob_scale_recovered(spacing, radius):
    vol = blob_volume(radius, spacing)
    b = F.blob_scale(vol, spacing, [2.0, 4.0, 8.0, 16.0])
    assert b['status'] == 'OK' and b['radius_um'] == radius


@pytest.mark.parametrize('spacing', [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (1.0, 0.5, 0.5), (0.5, 1.0, 2.0)])
def test_physical_log_amplitude_matches_analytic(spacing):
    s, sigma, A = 6.0, 3.0, 1000.0
    shape = tuple(int(np.ceil(80 / sp)) | 1 for sp in spacing)
    z, y, x = np.indices(shape).astype(float)
    c = (np.asarray(shape) - 1) / 2
    r2 = sum(((g - cc) * sp) ** 2 for g, cc, sp in zip((z, y, x), c, spacing))
    img = A * np.exp(-r2 / (2 * s ** 2))
    got = F.physical_log(img, spacing, sigma)[tuple(int(cc) for cc in c)]
    want = -3 * A * s ** 3 / (s ** 2 + sigma ** 2) ** 2.5            # Laplacian of the smoothed Gaussian at its centre
    assert got == pytest.approx(want, rel=0.03)


def test_orientation_undefined_for_flat_and_all_masked():
    flat = np.full((16, 16, 16), 500, np.uint16)
    o = F.structure_tensor_orientation(flat, (1, 1, 1), 1.0, 3.0)
    assert o['status'] == 'NO_SIGNAL' and o['dominant_direction_sample_xyz'] is None and o['order_parameter'] is None
    vol = line_volume((1, 0, 0), (1, 1, 1), shape=(16, 16, 16))
    o2 = F.structure_tensor_orientation(vol, (1, 1, 1), 1.0, 3.0, exclude=np.ones(vol.shape, bool))
    assert o2['status'] == 'NO_VALID_VOXELS' and o2['order_parameter'] is None
    f = F.extract_features(flat, (1, 1, 1))
    assert f['orientation']['status'] == 'NO_SIGNAL' and f['bright_signal_fraction']['status'] == 'NO_SIGNAL'
    assert f['compact_object_scale']['status'] == 'NO_SIGNAL'
    g = F.extract_features(vol, (1, 1, 1), exclude=np.ones(vol.shape, bool))
    assert all(g[k]['status'] == 'NO_VALID_VOXELS' for k in ('intensity_qc', 'bright_signal_fraction', 'orientation',
                                                           'compact_object_scale'))


def test_support_flags_follow_resolution():
    root = F.support_flags([62.05, 19.23, 19.14])
    assert not root['soma_scale_sampling_adequate']
    fine = F.support_flags([0.97, 0.30, 0.30])
    assert fine['soma_scale_sampling_adequate']
    for flags in (root, fine, F.support_flags([0.01, 0.01, 0.01])):         # never certified by sampling alone
        assert flags['single_neurite_resolution'].startswith('NOT_CERTIFIED')
        assert not any(v is True for k, v in flags.items() if 'neurite' in k)
    feats = F.extract_features(line_volume((1, 0, 0), (62, 19, 19), shape=(8, 24, 24)), (62, 19, 19))
    assert feats['compact_object_scale'] == {'status': 'NOT_SUPPORTED_AT_THIS_RESOLUTION'}


# ---------------- naming and pins

def test_forbidden_names_refused():
    with pytest.raises(F.Refusal):
        F.check_names({'ok': {'axon_length': 1}})
    with pytest.raises(F.Refusal):
        F.check_names({'gene_profile': 1})
    F.check_names(F.extract_features(line_volume((0, 1, 0), (1, 1, 1), shape=(16, 16, 16)), (1, 1, 1)))


def test_obj_vertices(tmp_path):
    p = tmp_path / 'm.obj'
    p.write_text('# x\nv 1 2 3\nv 4 5 6\nf 1 2 2\n')
    assert np.allclose(F.read_obj_vertices(p), [[1, 2, 3], [4, 5, 6]])


@pytest.mark.skipif(not (F.RAW / 'transform.txt').exists(), reason='raw pilot files not present')
def test_pins_detect_tampering(tmp_path):
    F.verify_pins(F.RAW)
    for n in F.PINS:
        shutil.copy(F.RAW / n, tmp_path / n)
    with (tmp_path / 'transform.txt').open('a') as f:
        f.write(' ')
    with pytest.raises(F.Refusal):
        F.verify_pins(tmp_path)


# ---------------- strict KTX1 decoder

def make_ktx(vol_zyxc, kv, spacing_xyz=(0.5, 0.5, 1.0), origin=(100.0, 200.0, 300.0), **override):
    import struct
    d, h, w, c = vol_zyxc.shape
    ext = np.array([w, h, d]) * np.asarray(spacing_xyz)
    M = np.eye(4); M[:3, :3] = np.diag(ext); M[:3, 3] = origin
    corners = [tuple(np.asarray(origin) + ext * np.array(t)) for t in
               [(i, j, k) for k in (0, 1) for j in (0, 1) for i in (0, 1)]]
    meta = {'distance_units': 'micrometers', 'xyz_from_texcoord_xform': str(M), 'corner_xyzs': str(corners),
            'number_of_channels': str(c)}
    meta.update(kv)
    blob = b''
    for k, v in meta.items():
        raw = k.encode() + b'\0' + v.encode() + b'\0'
        blob += struct.pack('<I', len(raw)) + raw + b'\0' * (-len(raw) % 4)
    hdr = dict(endianness=0x04030201, glType=0x1401, glTypeSize=1, glFormat=0x8227, glInternalFormat=0x822B,
               glBaseInternalFormat=0x8227, w=w, h=h, d=d, arr=0, faces=0, mips=1, kvbytes=len(blob))
    size_override = override.pop('size', None)
    hdr.update(override)
    data = vol_zyxc.astype(np.uint8).tobytes()
    size = size_override if size_override is not None else len(data)
    return F.KTX_IDENTIFIER + struct.pack('<13I', *hdr.values()) + blob + struct.pack('<I', size) + data


def test_ktx_roundtrip_and_geometry(tmp_path):
    rng = np.random.default_rng(0)
    vol = rng.integers(0, 255, (5, 6, 8, 2), dtype=np.uint8)
    p = tmp_path / 't.ktx'; p.write_bytes(make_ktx(vol, {'specimen_id': 'x'}))
    got, hd, kv = F.read_ktx1(p)
    assert np.array_equal(got, vol) and kv['specimen_id'] == 'x' and hd['pixelWidth'] == 8
    g = F.ktx_geometry(kv, got.shape[:3])
    assert np.allclose(g['spacing_um_xyz'], [0.5, 0.5, 1.0]) and np.allclose(g['origin_um_xyz'], [100, 200, 300])


@pytest.mark.parametrize('case', ['identifier', 'endianness', 'format', 'size', 'truncated', 'kv_length', 'faces'])
def test_ktx_strict_refusals(tmp_path, case):
    vol = np.zeros((2, 2, 4, 2), np.uint8)
    over = {'endianness': {'endianness': 0x01020304}, 'format': {'glInternalFormat': 0x8229},
            'size': {'size': 7}, 'kv_length': {'kvbytes': 3}, 'faces': {'faces': 6}}.get(case, {})
    b = make_ktx(vol, {}, **over)
    if case == 'identifier':
        b = b'XKTX' + b[4:]
    elif case == 'truncated':
        b = b[:-5]
    p = tmp_path / 'b.ktx'; p.write_bytes(b)
    with pytest.raises(F.Refusal):
        F.read_ktx1(p)


@pytest.mark.parametrize('case', ['units', 'corners'])
def test_ktx_geometry_refusals(tmp_path, case):
    vol = np.zeros((2, 2, 4, 2), np.uint8)
    kv = {'distance_units': 'millimeters'} if case == 'units' else {'corner_xyzs': str([(0.0, 0.0, 0.0)] * 8)}
    p = tmp_path / 'g.ktx'; p.write_bytes(make_ktx(vol, kv))
    _, _, meta = F.read_ktx1(p)
    with pytest.raises(F.Refusal):
        F.ktx_geometry(meta, (2, 2, 4))


@pytest.mark.skipif(not (F.KTX_ROUTE / F.KTX_PIN['file']).exists(), reason='ROOT KTX source not present')
def test_real_ktx_metadata_matches_receipt():
    vol, hd, kv = F.read_ktx1(F.KTX_ROUTE / F.KTX_PIN['file'])
    assert vol.shape == (150, 144, 224, 2) and vol.dtype == np.uint8
    g = F.ktx_geometry(kv, vol.shape[:3])
    assert np.allclose(g['spacing_um_xyz'], [0.6633740234375, 0.6529878472222221, 0.9988425], rtol=1e-8)
    assert kv['number_of_channels'] == '2' and 'rescaled intensity to 8 bits' in kv['relation_to_parent']


# ---------------- planar sheet versus tube (degenerate smallest eigenvalues)

def sheet_volume(normal_xyz, spacing_zyx, shape=(40, 64, 64), width_um=3.0):
    n = np.asarray(normal_xyz, float); n /= np.linalg.norm(n)
    z, y, x = np.indices(shape).astype(float)
    pz, py, px = z * spacing_zyx[0], y * spacing_zyx[1], x * spacing_zyx[2]
    c = np.array([px.mean(), py.mean(), pz.mean()])
    dist = np.abs((px - c[0]) * n[0] + (py - c[1]) * n[1] + (pz - c[2]) * n[2])
    return (1000 * np.exp(-(dist / width_um) ** 2) + 10).astype(np.uint16)


@pytest.mark.parametrize('normal', [(0, 0, 1), (1, 0, 0), (1, 1, 0)])
@pytest.mark.parametrize('spacing', [(1.0, 1.0, 1.0), (3.0, 1.0, 1.0)])
def test_planar_sheet_has_no_tube_direction(normal, spacing):
    o = F.structure_tensor_orientation(sheet_volume(normal, spacing), spacing, sigma_um=2.0, rho_um=4.0)
    assert o['status'] == 'NO_TUBE_DIRECTION' and o['dominant_direction_sample_xyz'] is None
    assert o['tube_gap_ratio'] < F.TUBE_GAP_MIN and o['planar_ratio'] > 0.8


@pytest.mark.parametrize('direction', [(0, 0, 1), (1, 0, 0), (0, 1, 1)])
def test_tube_has_large_gap_and_correct_direction(direction):
    o = F.structure_tensor_orientation(line_volume(direction, (1, 1, 1)), (1, 1, 1), sigma_um=2.0, rho_um=4.0)
    d = np.asarray(direction, float) / np.linalg.norm(direction)
    assert o['status'] == 'OK' and o['tube_gap_ratio'] > 0.5
    assert abs(np.dot(o['dominant_direction_sample_xyz'], d)) > 0.95


def test_tube_inside_weaker_sheet_still_points_along_tube():
    sp = (1.0, 1.0, 1.0)
    vol = (line_volume((1, 0, 0), sp).astype(float) + 0.3 * sheet_volume((0, 0, 1), sp, width_um=6.0)).astype(np.uint16)
    o = F.structure_tensor_orientation(vol, sp, sigma_um=2.0, rho_um=4.0)
    assert o['status'] == 'OK' and abs(o['dominant_direction_sample_xyz'][0]) > 0.9
