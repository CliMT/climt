"""Serve the climt pure wheel over a CORS-enabled localhost server.

The in-browser demo pages fetch the wheel from a different origin (port) than
the Quarto preview server, so the wheel host must send
`Access-Control-Allow-Origin`. Python's plain `http.server` does not, so use
this instead:

    python docs/radiative-transfer/_live/serve_wheel.py <dir-with-the-wheel> [port]

e.g., after building the pure wheel into /tmp/climt_wh:

    CLIMT_PURE_PYTHON=1 python -m pip wheel . --no-deps -w /tmp/climt_wh
    python docs/radiative-transfer/_live/serve_wheel.py /tmp/climt_wh 8912

Leave it running, then `quarto preview` the docs. The page front matter's
`pyodide: packages:` URL points at http://127.0.0.1:8912/climt-...whl.
"""
import functools
import http.server
import sys

DIRECTORY = sys.argv[1] if len(sys.argv) > 1 else "."
PORT = int(sys.argv[2]) if len(sys.argv) > 2 else 8912


class CORSHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


handler = functools.partial(CORSHandler, directory=DIRECTORY)
with http.server.ThreadingHTTPServer(("127.0.0.1", PORT), handler) as httpd:
    print(f"serving {DIRECTORY} on http://127.0.0.1:{PORT} (CORS *) — Ctrl-C to stop")
    httpd.serve_forever()
