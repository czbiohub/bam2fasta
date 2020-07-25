"""
bam2fasta command line.
"""
from __future__ import print_function
import sys
import argparse
import logging

from bam2fasta.cli import (
    percell, info, count_umis_percell, make_fastqs_percell)

usage = '''
bam2fasta   <command> [<args>]

** Commands include:
percell
info
count_umis_percell
make_fastqs_percell
'''


def main():
    """Group the bam2fasta commands
    currently percell, info are the commands under bam2fasta cli tool"""
    logging.basicConfig(
        format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s', # noqa
        datefmt='%Y-%m-%d:%H:%M:%S',
        stream=sys.stdout,
        level=logging.INFO)

    commands = {
        # keeping back compatibility
        # bam2fasta 1.0.1 conversion of bam2fsta was called convert
        'convert': percell,
        'percell': percell,
        'info': info,
        'count_umis_percell': count_umis_percell,
        'make_fastqs_percell': make_fastqs_percell}
    parser = argparse.ArgumentParser(
        description='conversion of 10x .bam file to several .fasta files')
    parser.add_argument('command', nargs='?')
    args = parser.parse_args(sys.argv[1:2])

    if not args.command:
        print(usage)
        sys.exit(1)

    if args.command not in commands:
        AssertionError('Unrecognized command')
        print(usage)
        sys.exit(1)

    cmd = commands.get(args.command)
    cmd(sys.argv[2:])


if __name__ == '__main__':
    main()
