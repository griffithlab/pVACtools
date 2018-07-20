from http.server import SimpleHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse
import threading
import os
import webbrowser
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

    Handler = MyHandler
    httpd = HTTPServer((HOSTNAME, PORT), Handler)
    thread = threading.Thread(target=httpd.serve_forever)

    try:
        print(time.asctime(), "Starting pVACviz client webserver")
        thread.start()

    except (KeyboardInterrupt, SystemExit):
        print(time.asctime(), "Stopping pVACviz client webserver")
        httpd.shutdown()
        pass

    except OSError as err:
        print("OS error while starting pVACviz client server: {0}".format(err))
        httpd.shutdown()
        pass

    finally:
        print(time.asctime(), "pVACviz server started at http://%s:%s" % (HOSTNAME, PORT))

    print(time.asctime(), "Opening pVACviz client at http://%s:%s in default browser." % (HOSTNAME, PORT))
    webbrowser.get().open("http://%s:%s" % (HOSTNAME, PORT), new=1, autoraise=True)


if __name__ == "__main__":
    main()
