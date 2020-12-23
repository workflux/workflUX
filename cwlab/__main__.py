from __future__ import absolute_import
import argparse

def main():
    parser = argparse.ArgumentParser(
        prog="cwlab",
        description='CWLab: A platform-agnostic, cloud-ready framework for simplified' + \
            ' deployment of the Common Workflow Language using a graphical web interface '
    )

    parser.add_argument(
        '--version','-v',
        help="Return version of cwlab.",
        action='store_true'
    )

    subparser = parser.add_subparsers(
        help="CWLab sub-commands",
        dest='subcommand'
    )

    parser_up = subparser.add_parser(
        "up",
        help="Start the webserver."
    )
    parser_up.add_argument(
        '-c', '--config',
        help="Specify the path to a costum config file."
    )

    parser_up = subparser.add_parser(
        "print_config",
        help="Get an example config. Typical usage: cwlab print_config > example_config.yaml"
    )

    args = parser.parse_args()

    if args.version:
        from . import __version__
        print(f"cwlab {__version__}")
    if args.subcommand == "up":
        from . import create_app
        create_app(config_file=args.config, webapp=True)
    elif args.subcommand == "print_config":
        from .utils import output_example_config
        output_example_config()

if __name__ == "__main__":
    main()