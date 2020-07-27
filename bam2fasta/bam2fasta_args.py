import argparse

DEFAULT_LINE_COUNT = 1500
DEFAULT_DELIMITER = "X"
DEFAULT_PROCESSES = 2
DEFAULT_MIN_UMI_PER_BARCODE = 0
CELL_BARCODE_PATTERN = r'(CB|XC):Z:(.*)[\t]'
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
    parser.add_argument(
        '--filename', type=str, help="10x bam file or fastq.gz file")

    parser.add_argument(
        '--min-umi-per-barcode', default=DEFAULT_MIN_UMI_PER_BARCODE, type=int,
        help="A barcode is only considered a valid barcode read "
        "and its fasta is written if number of umis are greater "
        "than min-umi-per-barcode. It is used to weed out cell barcodes "
        "with few umis that might have been due to false rna enzyme reactions",
        required=False)
    parser.add_argument(
        '--write-barcode-meta-csv', type=str,
        help="For each of the unique barcodes, "
        "Write to a given path, barcode and"
        "and number of umis per barcode.")
    parser.add_argument(
        '-p', '--processes', default=DEFAULT_PROCESSES, type=int,
        help='Number of processes to use for reading 10x bam file',
        required=False)
    parser.add_argument(
        '--delimiter', default=DEFAULT_DELIMITER, type=str,
        help='delimiter between sequences of the same barcode',
        required=False)
    parser.add_argument(
        '--save-fastas', default="", type=str,
        help='save merged fastas for all the unique'
        'barcodes to {CELL_BARCODE}.fasta '
        'in the absolute path given by this flag'
        'By default, fastas are saved in current directory',
        required=False)
    parser.add_argument(
        '--save-intermediate-files', default="/tmp/", type=str,
        help='save temporary fastas and chunks of bam files'
        'in the absolute path given by this flag'
        'By default, they are saved in temp directory.'
        'An important note is This might cause'
        'not enough space on the device left depending on the size of your'
        'bam file and harddisk space allocated for tmp folder on your machine,'
        'so its better to specify a directory.'
        'These files are deleted automatically at the end of the program',
        required=False)
    parser.add_argument(
        '--shard-size', '--line-count', type=int,
        help='Line/Alignment count for each bam shard, higher the count lesser'
        'the number of temporary bam files and shards created, increase the'
        'default number according to the size of the bam file',
        default=DEFAULT_LINE_COUNT, required=False)
    parser.add_argument(
        "--cell-barcode-pattern", type=str,
        help="Regular expressions for cell barcodes",
        default=CELL_BARCODE_PATTERN, required=False)
    parser.add_argument(
        "--molecular-barcode-pattern", type=str,
        help="Regular expressions for molecular barcodes.",
        default=MOLECULAR_BARCODE_PATTERN, required=False)
    parser.add_argument(
        '--rename-10x-barcodes', type=str,
        help="Tab-separated file mapping 10x barcode name to new name"
        "e.g. with channel or cell "
        "annotation label", required=False)
    parser.add_argument(
        '--barcodes-file', type=str,
        help="Barcodes file if the input is unfiltered 10x bam file",
        required=False)
    parser.add_argument(
        '--barcodes-significant-umis-file', type=str,
        help="Barcodes file with significant umi count",
        required=False)
    parser.add_argument(
        "--channel-id", type=str,
        help="Output prefix for fastqs", default="", required=False)
    parser.add_argument(
        '--output-format', type=str,
        help="Output prefix for fastqs, can be either fastq or fastq.gz",
        required=False, default="fastq")
    return parser
