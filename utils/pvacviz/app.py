from http.server import SimpleHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse
import os
import time

HOSTNAME = 'localhost'
PORT = 4200
INDEXFILE = 'index.html'
CLIENTDIR = os.path.join(os.path.dirname(__file__), 'client')


class MyHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        # Parse query data to find out what was requested
        parsedParams = urlparse(self.path)

        # See if the file requested exists
        if os.access(CLIENTDIR + os.sep + parsedParams.path, os.R_OK):
            # File exists, serve it up
            SimpleHTTPRequestHandler.do_GET(self)
        else:
            # send index.html, but don't redirect
            self.send_response(200)
            self.send_header('Content-Type', 'text/html')
            self.end_headers()
            with open(INDEXFILE, 'rb') as fin:
                self.copyfile(fin, self.wfile)


def main():
    os.chdir(CLIENTDIR)

    Handler = MyHandler
    httpd = HTTPServer((HOSTNAME, PORT), Handler)

    print(time.asctime(), "Server Starts - %s:%s" % (HOSTNAME, PORT))

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
