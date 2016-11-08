#!/usr/bin/env python3

import connexion
import os
import sys
import http.server
import socketserver
from threading import Thread
from webbrowser import open_new_tab

def check_is_directory(directory):
    fullpath = os.path.abspath(directory)
    if os.path.isfile(fullpath):
        raise argparse.ArgumentTypeError("Path \"%s\" must be a directory"%directory)
    if not os.path.isdir(fullpath):
        os.makedirs(fullpath)
        # raise argparse.ArgumentTypeError("Path \"%s\" must be a directory"%directory)
    return fullpath

def main():
    app = connexion.App("pVAC-Seq Visualization Server", specification_dir=os.path.join(
        os.path.dirname(__file__),
        'swagger'
    ))
    app.add_api('swagger.yaml', arguments={'title': 'API to support pVacSeq user interface for generating reports on pipeline results'})
    app.app.secret_key = os.urandom(1024)

    #Eventually, have this open a browser to whatever the main page is
    # Thread(target=lambda:open_new_tab("localhost:8080/static/testpage"), daemon=True).start()
    app.run(port=8080)

if __name__ == '__main__':
    main()
