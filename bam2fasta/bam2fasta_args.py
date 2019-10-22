import argparse


class Bam2FastaArgumentParser(argparse.ArgumentParser):
    """Specialize ArgumentParser for bam2Fasta.

    In particular, set up sharing of variables
    """
    def __init__(self, no_citation=False, **kwargs):
        super(Bam2FastaArgumentParser, self).__init__(add_help=False, **kwargs)

    def parse_args(self, args=None, namespace=None):
        args = super(Bam2FastaArgumentParser, self).parse_args(
            args=args,
            namespace=namespace)
        return args
