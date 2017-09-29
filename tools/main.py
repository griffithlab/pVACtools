import argparse
import sys
import pkg_resources

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    #add subcommands
    parser.add_argument(
        "-v", "--version",
        action="store_true",
        help="Display the currently installed pvactools version",
    )

    args = parser.parse_known_args()
    if args[0].version is True:
        print(pkg_resources.get_distribution("pvactools").version)
    else:
        try:
            args[0].func.main(args[1])
        except AttributeError as e:
            parser.print_help()
            print("Error: No command specified")
            sys.exit(-1)


if __name__ == '__main__':
    main()
