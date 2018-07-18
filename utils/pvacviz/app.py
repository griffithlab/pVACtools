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
    thread = threading.Thread(target=httpd.serve_forever)

    try:
        print(time.asctime(), "%s:%s - Starting pVACviz client webserver" % (HOSTNAME, PORT))
        thread.start()

    except (KeyboardInterrupt, SystemExit):
        print(time.asctime(), "%s:%s - Stopping pVACviz client webserver" % (HOSTNAME, PORT))
        httpd.shutdown()
        pass

    except OSError as err:
        print("OS error while starting pVACviz client server: {0}".format(err))
        httpd.shutdown()
        pass

    finally:
        print(time.asctime(), "%s:%s - pVACviz server started." % (HOSTNAME, PORT))

    # import rpdb; rpdb.set_trace()
    print("Opening pVACviz client in default browser.")
    webbrowser.open("%s:%s" % (HOSTNAME, PORT), new=1, autoraise=True)


if __name__ == "__main__":
    main()
