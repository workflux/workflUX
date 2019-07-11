#!/usr/bin/env python
import argparse

def run():
    parser = argparse.ArgumentParser(
        prog="cwlab",
        description='CWLab: A platform-agnostic, cloud-ready framework for simplified' + \
            ' deployment of the Common Workflow Language using a graphical web interface '
    )

    subparser = parser.add_subparsers(
        help="CWLab sub-commands",
        dest='subcommand'
    )

    parser_up = subparser.add_parser(
        "up",
        help="Start the webserver."
    )
    parser_up.add_argument('-c', '--config',
        help="Specify the path to a costum config file."
    )

    args = parser.parse_args()

    if args.subcommand == "up":
        from . import up
        up(config_file=args.config)

if __name__ == "__main__":
    run()