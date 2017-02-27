#!/usr/bin/env python3

import connexion
import os
import sys
import http.server
import socketserver
from threading import Thread
from webbrowser import open_new_tab
import atexit
from flask_cors import CORS
from postgresql.exceptions import UndefinedTableError

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
    app.app.config['initialized'] = False

    def cleanup():
        print("Cleaning up observers and database connections")
        if app.app.config['initialized']:
            app.app.config['storage']['watcher'].stop()
            app.app.config['storage']['watcher'].join()
            if 'db-clean' in app.app.config:
                for table in app.app.config['db-clean']:
                    try:
                        app.app.config['storage']['db'].execute("DROP TABLE %s"%table)
                    except UndefinedTableError:
                        pass
            app.app.config['storage']['db'].close()

    atexit.register(cleanup)
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
