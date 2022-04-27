#! /usr/bin/env python

import logging


try:
    import CIAlign.argP as argP
    import CIAlign.runCIAlign as runCIAlign
except ImportError:
    import argP
    import runCIAlign


def main():
    # Get and parse the argument parser
    parser = argP.getParser()
    args = parser.parse_args()

    # Set up logger
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    logfile = "%s_log.txt" % args.outfile_stem
    handler = logging.FileHandler(logfile)
    handler.setLevel(logging.INFO)

    # Create a logging format
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    log.addHandler(handler)

    # add the handlers to the logger
    log.addHandler(handler)

    log.info("\nInitial parameters:\n%s" % str(parser.format_values()))
    runCIAlign.run(args, log)


if __name__ == "__main__":
    main()
