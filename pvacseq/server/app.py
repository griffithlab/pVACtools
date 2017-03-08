#!/usr/bin/env python3

import connexion
import os
import sys
import http.server
import socketserver
from threading import Thread
from webbrowser import open_new_tab
from flask_cors import CORS
from .controllers.utils import initialize
#FIXME: sanitize sample name
def main():

    app = connexion.App(
        "pVAC-Seq Visualization Server",
        specification_dir=os.path.join(
            os.path.dirname(__file__),
            'config'
        ),
    )

    from werkzeug.routing import IntegerConverter as BaseIntConverter
    class IntConverter(BaseIntConverter):
        regex = r'-?\d+'

    app.app.url_map.converters['int'] = IntConverter
    initialize(app.app) #initialize the app configuration
    app.add_api('swagger.yaml', arguments={'title': 'API to support pVacSeq user interface for generating reports on pipeline results'})
    app.app.secret_key = os.urandom(1024)

    #setup CORS
    CORS(
        app.app,
        #should match localhost at with any port, path, or protocol
        origins=r'^(.+://)?localhost(:\d+)?(/.*)?$'
    )

    #Eventually, have this open a browser to whatever the main page is
    # Thread(target=lambda:open_new_tab("localhost:8080/static/testpage"), daemon=True).start()
    app.run(port=8080)

if __name__ == '__main__':
    main()
