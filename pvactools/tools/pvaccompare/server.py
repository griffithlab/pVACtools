from http.server import SimpleHTTPRequestHandler, HTTPServer
from pathlib import Path
import sys
import os


class CORSRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header(
            "Access-Control-Allow-Headers", "x-requested-with, content-type"
        )
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.end_headers()


def main():
    file_path = Path(__file__).resolve()
    html_dir = file_path.parent / "html_report"

    # Make the server serve files from the html_report directory
    os.chdir(html_dir)

    PORT = 8080
    server = HTTPServer(("localhost", PORT), CORSRequestHandler)
    print("Starting local server...")
    print(
        f"View reports at http://localhost:{PORT}/main.html"
    )
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server.")
        server.server_close()
        sys.exit(0)


if __name__ == "__main__":
    main()
