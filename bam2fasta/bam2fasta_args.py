import argparse

DEFAULT_PROCESSES = 2
DEFUALT_MIN_UMI_PER_BARCODE = 100
CELL_BARCODE_PATTERN = r'(CB|XC):Z:([ACGT]+)\-1'
MOLECULAR_BARCODE_PATTERN = '(UB|XB):Z:([ACGT]+)'


class Bam2FastaArgumentParser(argparse.ArgumentParser):
    """Specialize ArgumentParser for bam2Fasta."""
    def __init__(self, no_citation=False, **kwargs):
        super(Bam2FastaArgumentParser, self).__init__(add_help=False, **kwargs)

    def parse_args(self, args=None, namespace=None):
        args = super(Bam2FastaArgumentParser, self).parse_args(
            args=args,
            namespace=namespace)
        return args


def create_parser():
    """Returns after adding all arguments to Bam2FastaArgumentParser."""
    parser = Bam2FastaArgumentParser()
    parser.add_argument('--filename', type=str, help="10x bam file")

    parser.add_argument(
        '--min-umi-per-barcode', default=DEFUALT_MIN_UMI_PER_BARCODE, type=int,
        help="A barcode is only considered a valid barcode read "
        "and its fasta is written if number of umis are greater "
        "than min-umi-per-barcode. It is used to weed out cell barcodes "
        "with few umis that might have been due to false rna enzyme reactions")
    parser.add_argument(
        '--write-barcode-meta-csv', type=str,
        help="For each of the unique barcodes, "
        "Write to a given path, number of reads"
        "and number of umis per barcode.")
    parser.add_argument(
        '-p', '--processes', default=DEFAULT_PROCESSES, type=int,
        help='Number of processes to use for reading 10x bam file')
    parser.add_argument(
        '--save-fastas', default="", type=str,
        help='save merged fastas for all the unique'
        'barcodes to {CELL_BARCODE}.fasta '
        'in the absolute path given by this flag'
        'By default, fastas are saved in current directory')
    parser.add_argument(
        '--save-intermediate-files', default="", type=str,
        help='save intermediate fasta gzips')
    parser.add_argument(
        "--cell-barcode-pattern", type=str,
        help="Regular expressions for cell barcodes. Default is"
        " 10x Genomics 'CB:Z' tag",
        default=CELL_BARCODE_PATTERN)
    parser.add_argument(
        "--molecular-barcode-pattern", type=str,
        help="Regular expressions for molecular barcodes. "
             "Default is 10x Genomics 'UB:Z' tag",
        default=MOLECULAR_BARCODE_PATTERN)
    parser.add_argument(
        '--rename-10x-barcodes', type=str,
        help="Tab-separated file mapping 10x barcode name to new name"
        "e.g. with channel or cell "
        "annotation label", required=False)
    parser.add_argument(
        '--barcodes-file', type=str,
        help="Barcodes file if the input is unfiltered 10x bam file",
        required=False)
    return parser
