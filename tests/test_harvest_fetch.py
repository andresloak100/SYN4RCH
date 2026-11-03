"""Regression: a dropped connection is retried, not fatal to the harvest run."""
import http.client
import unittest
from unittest import mock

from _path import *  # noqa: F401,F403  (puts src/ on sys.path)

from syn4rch.harvest import cache as cmod


class _Resp:
    status = 200

    def __init__(self, body):
        self._body = body

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class FetchRetryTest(unittest.TestCase):
    def setUp(self):
        self.sleep = mock.patch.object(cmod.time, "sleep", lambda s: None).start()
        mock.patch.object(cmod, "_throttle", lambda url: None).start()
        self.addCleanup(mock.patch.stopall)

    def test_remote_disconnect_is_retried_then_succeeds(self):
        calls = {"n": 0}

        def flaky(req, timeout=None):
            calls["n"] += 1
            if calls["n"] == 1:
                raise http.client.RemoteDisconnected("closed without response")
            return _Resp(b"<ok/>")

        with mock.patch.object(cmod.urllib.request, "urlopen", flaky):
            body, cached = cmod.fetch("https://example.org/x", source="t", query="q")
        self.assertEqual(body, "<ok/>")
        self.assertFalse(cached)
        self.assertEqual(calls["n"], 2)

    def test_persistent_disconnect_becomes_fetcherror(self):
        def dead(req, timeout=None):
            raise http.client.RemoteDisconnected("closed without response")

        with mock.patch.object(cmod.urllib.request, "urlopen", dead):
            with self.assertRaises(cmod.FetchError):
                cmod.fetch("https://example.org/x", source="t", query="q", max_retries=2)

    def test_connection_reset_becomes_fetcherror(self):
        def reset(req, timeout=None):
            raise ConnectionResetError("reset by peer")

        with mock.patch.object(cmod.urllib.request, "urlopen", reset):
            with self.assertRaises(cmod.FetchError):
                cmod.fetch("https://example.org/x", source="t", query="q", max_retries=1)


if __name__ == "__main__":
    unittest.main()
