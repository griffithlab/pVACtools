#!/usr/bin/env python3

import connexion
import os
import sys
import time
import argparse

from flask_cors import CORS
from utils.pvacapi.controllers.utils import initialize
from utils.pvacapi.controllers.utils import getIpAddress

#FIXME: sanitize sample name
def main():
    parser = argparse.ArgumentParser(description='parse pvacapi arguments')
    parser.add_argument('--ip_address', help='IP address for the HTTP server to bind')
    args = parser.parse_args()

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

    # determine IP address and setup CORS
    IPAddr = None
    if args.ip_address is None:
        IPAddr = getIpAddress()
    else:
        IPAddr = args.ip_address

    app.app.url_map.converters['int'] = IntConverter
    initialize(app.app, IPAddr) #initialize the app configuration
    app.add_api('swagger.yaml', arguments={'title': 'API to support pVacSeq user interface for generating reports on pipeline results'})
    app.app.secret_key = os.urandom(1024)

    CORS(
        app.app,
        # should match IP address at with any port, path, or protocol
        origins=r'^(.+://)?' + IPAddr + r'(:\d+)?(/.*)?$'
    )

    print(time.asctime(), "Starting pVACapi server at http://" + IPAddr + ":8080")
    app.run(port=8080, debug='--debug' in sys.argv, threaded=True)

if __name__ == '__main__':
    main()
