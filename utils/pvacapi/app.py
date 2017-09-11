#!/usr/bin/env python3

import connexion
import os
import sys
import http.server
import socketserver
from webbrowser import open_new_tab
from flask_cors import CORS
from utils.pvacapi.controllers.utils import initialize
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
    initialize(app.app, set(sys.argv)) #initialize the app configuration
    app.add_api('swagger.yaml', arguments={'title': 'API to support pVacSeq user interface for generating reports on pipeline results'})
    app.app.secret_key = os.urandom(1024)

    #setup CORS
    CORS(
        app.app,
        #should match localhost at with any port, path, or protocol
        origins=r'^(.+://)?localhost(:\d+)?(/.*)?$'
    )

    app.run(port=8080, debug='--debug' in sys.argv, threaded=True)

if __name__ == '__main__':
    main()
