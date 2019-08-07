#!/usr/bin/env python3

import connexion
import os
import sys
import time
import argparse

from flask_cors import CORS
from flask import g
from utils.pvacapi.controllers.utils import initialize
from utils.pvacapi.controllers.utils import getIpAddress

def app_parser():
    parser = argparse.ArgumentParser(description='pVACapi provides a REST API to pVACtools')
    parser.add_argument('--ip-address', help='IP address for the HTTP server to bind. If not provided, the default socket address will be used.')
    parser.add_argument('--proxy-ip-address', help='IP address of proxy server or public IP address. If provided, server will send X-Forward headers required for Bokeh to properly work through a proxy server or with AWS private/public IP addresses.')
    parser.add_argument('--debug', default=False, action='store_true', help='Start sever in debug mode.')
    return parser

#FIXME: sanitize sample name
def main(args = None):
    if not args:
        parser = app_parser()
        args = parser.parse_args(sys.argv[1:])

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
    IP_ADDRESS = None
    if args.ip_address is None:
        IP_ADDRESS = getIpAddress()
    else:
        IP_ADDRESS = args.ip_address

    app.app.IP_ADDRESS = IP_ADDRESS

    # add forwarding headers if proxy_ip_address specified
    PROXY_IP_ADDRESS = None
    if args.proxy_ip_address is not None:
        PROXY_IP_ADDRESS = args.proxy_ip_address
        @app.app.after_request
        def apply_forwarding_headers(response):
            response.headers["X-Forwarded-Proto"] = "http"
            response.headers["X-Forwarded-Host"] = PROXY_IP_ADDRESS
            response.headers["X-Forwarded-For"] = IP_ADDRESS
            response.headers["X-Real-IP"] = IP_ADDRESS
            return response

    app.app.PROXY_IP_ADDRESS = PROXY_IP_ADDRESS

    app.app.url_map.converters['int'] = IntConverter
    initialize(app.app) #initialize the app configuration
    app.add_api('swagger.yaml', arguments={'title': 'API to support pVacSeq user interface for generating reports on pipeline results'})
    app.app.secret_key = os.urandom(1024)

    # remove all CORS restrictions
    CORS(app.app)

    # should match IP address at with any port, path, or protocol
    # CORS(
    #     app.app,
    #     origins=r'^(.+://)?' + IP_ADDRESS + r'(:\d+)?(/.*)?$'
    # )

    print(time.asctime(), "Starting pVACapi server at http://" + IP_ADDRESS + ":8080")
    app.run(port=8080, debug=args.debug, threaded=True)

if __name__ == '__main__':
    main()
