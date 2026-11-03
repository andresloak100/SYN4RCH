"""OBJ parsing and vertex-clustering simplification."""

import _path  # noqa: F401
import json
import unittest

from syn4rch.atlas import build, mesh


def grid_sphere(n=40, r=500.0):
    """UV sphere as a quad grid, written as OBJ text with a//n face tokens."""
    import math
    lines = []
    for i in range(n + 1):
        th = math.pi * i / n
        for j in range(n):
            ph = 2 * math.pi * j / n
            lines.append("v %f %f %f" % (r * math.sin(th) * math.cos(ph) + 5000,
                                         r * math.cos(th) + 4000, r * math.sin(th) * math.sin(ph) + 6000))
    for i in range(n):
        for j in range(n):
            a = i * n + j + 1
            b = i * n + (j + 1) % n + 1
            c = (i + 1) * n + (j + 1) % n + 1
            d = (i + 1) * n + j + 1
            lines.append("f %d//%d %d//%d %d//%d %d//%d" % (a, a, b, b, c, c, d, d))
    return lines


class ObjTest(unittest.TestCase):
    def test_parse_triangles_and_quads(self):
        v, t = mesh.parse_obj(["# c", "v 0 0 0", "v 1 0 0", "v 1 1 0", "v 0 1 0",
                               "vn 0 0 1", "f 1//1 2//2 3//3", "f 1 2 3 4"])
        self.assertEqual(len(v), 4)
        self.assertEqual(t, [(0, 1, 2), (0, 1, 2), (0, 2, 3)])

    def test_negative_indices(self):
        v, t = mesh.parse_obj(["v 0 0 0", "v 1 0 0", "v 0 1 0", "f -3 -2 -1"])
        self.assertEqual(t, [(0, 1, 2)])


class ClusterTest(unittest.TestCase):
    def setUp(self):
        self.v, self.t = mesh.parse_obj(grid_sphere())

    def test_reduces_and_stays_valid(self):
        v2, t2 = mesh.cluster(self.v, self.t, 100.0)
        self.assertLess(len(v2), len(self.v))
        for a, b, c in t2:
            self.assertTrue(len({a, b, c}) == 3)
            self.assertTrue(max(a, b, c) < len(v2))
        self.assertEqual(len(set(tuple(sorted(x)) for x in t2)), len(t2))
        used = {i for tri in t2 for i in tri}
        self.assertEqual(used, set(range(len(v2))))

    def test_geometry_stays_within_original_bbox(self):
        lo, hi = mesh.bbox(self.v)
        v2, _ = mesh.cluster(self.v, self.t, 150.0)
        lo2, hi2 = mesh.bbox(v2)
        for k in range(3):
            self.assertGreaterEqual(lo2[k], lo[k] - 1e-6)
            self.assertLessEqual(hi2[k], hi[k] + 1e-6)
            self.assertLess(hi[k] - hi2[k], 150.0 * 1.8)

    def test_simplify_to_budget(self):
        v2, t2, cell = mesh.simplify_to(self.v, self.t, 300, min_cell=5.0)
        self.assertLessEqual(len(v2), 300)
        self.assertGreater(len(v2), 100)
        self.assertGreater(cell, 5.0)

    def test_small_mesh_untouched_at_min_cell(self):
        v2, t2, cell = mesh.simplify_to(self.v, self.t, 10 ** 6, min_cell=1.0)
        self.assertEqual(cell, 1.0)
        self.assertEqual(len(v2), len(self.v) - 2 * 39)  # pole duplicates merge

    def test_rejects_bad_cell(self):
        with self.assertRaises(ValueError):
            mesh.cluster(self.v, self.t, 0)

    def test_encode_is_compact_integer_json(self):
        v2, t2 = mesh.cluster(self.v, self.t, 120.0)
        enc = mesh.encode(v2, t2)
        self.assertEqual(len(enc["p"]), 3 * enc["vertex_count"])
        self.assertEqual(len(enc["i"]), 3 * enc["triangle_count"])
        self.assertTrue(all(isinstance(x, int) for x in enc["p"]))
        json.dumps(enc)

    def test_inline_json_cannot_close_script(self):
        s = build.inline_json({"a": "</script><!--"})
        self.assertNotIn("</", s)
        self.assertNotIn("<!--", s)
        self.assertEqual(json.loads(s)["a"], "</script><!--")


if __name__ == "__main__":
    unittest.main()
