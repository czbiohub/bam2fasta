"""
sourmash command line.
"""
from __future__ import print_function
import sys
import argparse
import logging

from bam2fasta.cli import convert

usage = '''
bam2fasta   <command> [<args>]

** Commands include:
convert
'''


def main():
    logging.basicConfig(
        format='%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s',
        datefmt='%Y-%m-%d:%H:%M:%S',
        stream=sys.stdout,
        level=logging.INFO)

    commands = {'convert': convert}
    parser = argparse.ArgumentParser(
        description='work with conversion of 10x .bam file to several .fasta files')
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
