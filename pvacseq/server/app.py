#!/usr/bin/env python3

import connexion
import os
import sys
import http.server
import socketserver
from threading import Thread
from webbrowser import open_new_tab
import postgresql as psql
from postgresql.exceptions import UndefinedTableError
import atexit
from .controllers import watchdir

def check_is_directory(directory):
    fullpath = os.path.abspath(directory)
    if os.path.isfile(fullpath):
        raise argparse.ArgumentTypeError("Path \"%s\" must be a directory"%directory)
    if not os.path.isdir(fullpath):
        os.makedirs(fullpath)
        # raise argparse.ArgumentTypeError("Path \"%s\" must be a directory"%directory)
    return fullpath

def main():

    app = connexion.App(
        "pVAC-Seq Visualization Server",
        specification_dir=os.path.join(
            os.path.dirname(__file__),
            'swagger'
        ),
    )
    from werkzeug.routing import IntegerConverter as BaseIntConverter

    class IntConverter(BaseIntConverter):
        regex = r'-?\d+'

    # before routes are registered
    app.app.url_map.converters['int'] = IntConverter
    app.app.config['dropbox_dir'] = os.path.expanduser("~/Desktop/PVACSEQ_TEST")
    os.makedirs(app.app.config['dropbox_dir'], exist_ok=True)

    watcher = watchdir.Observe(app.app.config['dropbox_dir'])
    watcher.subscribe(lambda x:print("FS Event:", x))
    tmp = psql.open("localhost/postgres")
    if not len(tmp.prepare("SELECT 1 FROM pg_database WHERE datname = 'pvacseq'")()):
        tmp.execute("CREATE DATABASE pvacseq")
    tmp.close()
    db = psql.open("localhost/pvacseq")
    def cleanup():
        print("Cleaning up observers and database connections")
        watcher.stop()
        watcher.join()
        if 'db-clean' in app.app.config:
            for table in app.app.config['db-clean']:
                try:
                    db.execute("DROP TABLE %s"%table)
                except UndefinedTableError:
                    pass
        db.close()

    atexit.register(cleanup)

    app.app.config['watcher'] = watcher
    app.app.config['db'] = db
    app.app.config['initialized'] = False

    app.add_api('swagger.yaml', arguments={'title': 'API to support pVacSeq user interface for generating reports on pipeline results'})
    app.app.secret_key = os.urandom(1024)

    #Eventually, have this open a browser to whatever the main page is
    # Thread(target=lambda:open_new_tab("localhost:8080/static/testpage"), daemon=True).start()
    app.run(port=8080)

if __name__ == '__main__':
    main()
