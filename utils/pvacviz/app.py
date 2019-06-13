from http.server import SimpleHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse
import threading
import os
import webbrowser
import time
import argparse
from utils.pvacapi.controllers.utils import getIpAddress

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
            try:
                SimpleHTTPRequestHandler.do_GET(self)
            except ConnectionError as err:
                print("%s occured attempting to open the path %s" % (err, parsedParams.path))
        else:
            # send index.html, but don't redirect
            self.send_response(200)
            self.send_header('Content-Type', 'text/html')
            self.end_headers()
            with open(INDEXFILE, 'rb') as fin:
                self.copyfile(fin, self.wfile)


def main():
    os.chdir(CLIENTDIR)
    parser = argparse.ArgumentParser(description='pVACviz serves the files for the pVACviz browser client')
    parser.add_argument('--ip-address', help='IP address to which the HTTP server will bind. If not provided, the default socket address will be used.')
    args = parser.parse_args()

    IP_ADDRESS = None
    if args.ip_address is None:
        IP_ADDRESS = getIpAddress()
    else:
        IP_ADDRESS = args.ip_address

    Handler = MyHandler
    httpd = HTTPServer((IP_ADDRESS, PORT), Handler)
    thread = threading.Thread(target=httpd.serve_forever)
    url = "http://{}:{}".format(IP_ADDRESS, PORT)

    try:
        print(time.asctime(), "Starting pVACviz client webserver")
        thread.start()
    except (KeyboardInterrupt, SystemExit):
        print(time.asctime(), "Stopping pVACviz client webserver")
        httpd.shutdown()
        return
    except OSError as err:
        print("OS error while starting pVACviz client server: {0}".format(err))
        httpd.shutdown()
        return

    print(time.asctime(), "pVACviz server started at {}".format(url))

    try:
        print(time.asctime(), "Opening pVACviz client at {} in default browser.".format(url))
        webbrowser.get().open(url, new=1, autoraise=True)
    except webbrowser.Error as err:
        print("No default browser found. Open pVACviz by visiting {} in a browser of your choosing.".format(url))

if __name__ == "__main__":
    main()
