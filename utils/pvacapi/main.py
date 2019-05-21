import argparse
import sys
from utils.pvacapi import *

def main():
    parser = app.app_parser()
    subparsers = parser.add_subparsers()

    clear_cache_parser = subparsers.add_parser(
        "clear_cache",
        help="Clear pVACtool's api process cache",
        add_help=False
    )
    clear_cache_parser.set_defaults(func=clear_cache)

    args = parser.parse_args()
    try:
        args.func.main()
    except AttributeError as e:
        app.main(args)


if __name__ == '__main__':
    main()
