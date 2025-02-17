from http.server import SimpleHTTPRequestHandler, HTTPServer
import sys


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


def main():
    PORT = 8080
    server = HTTPServer(("localhost", PORT), CORSRequestHandler)
    print("Starting local server...")
    print(f"View reports at http://localhost:{PORT}/pvactools/tools/pvaccompare/html_report/main.html")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server.")
        server.server_close()
        sys.exit(0)

if __name__ == "__main__":
    main()
