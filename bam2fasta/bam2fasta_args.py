import argparse

DEFAULT_LINE_COUNT = 1500
DEFAULT_DELIMITER = "X"
DEFAULT_PROCESSES = 2


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


def create_parser():
    parser = Bam2FastaArgumentParser()
    parser.add_argument('--filename', type=str, help="10x bam file")

    parser.add_argument('--count-valid-reads', default=0, type=int,
                        help="For 10x input only (i.e input-is-10x flag is True), "
                        "A barcode is only considered a valid barcode read "
                        "and its signature is written if number of umis are greater "
                        "than count-valid-reads. It is used to weed out cell barcodes "
                        "with few umis that might have been due to false rna enzyme reactions")
    parser.add_argument('--write-barcode-meta-csv', type=str,
                        help="For 10x input only (i.e input-is-10x flag is True), for each of the unique barcodes, "
                        "Write to a given path, number of reads and number of umis per barcode.")
    parser.add_argument('-p', '--processes', default=DEFAULT_PROCESSES, type=int,
                        help='For 10x input only (i.e input-is-10x flag is True, '
                        'Number of processes to use for reading 10x bam file')
    parser.add_argument('--delimiter', default=DEFAULT_DELIMITER, type=str,
                        help='delimiter between sequences')
    parser.add_argument('--save-fastas', default="", type=str,
                        help='For 10x input only (i.e input-is-10x flag is True), '
                        'save merged fastas for all the unique barcodes to {CELL_BARCODE}.fasta '
                        'in the absolute path given by this flag, By default, fastas are not saved')
    parser.add_argument('--line-count', type=int,
                        help='For 10x input only (i.e input-is-10x flag is True), line count for each bam shard',
                        default=DEFAULT_LINE_COUNT)
    parser.add_argument('--rename-10x-barcodes', type=str,
                        help="Tab-separated file mapping 10x barcode name "
                        "to new name, e.g. with channel or cell "
                        "annotation label", required=False)
    parser.add_argument('--barcodes-file', type=str,
                        help="Barcodes file if the input is unfiltered 10x bam file", required=False)
    return parser
