from . import bam2fasta_tst_utils as utils
from bam2fasta import cli
from bam2fasta import bam2fasta_args
import os


def test_iter_split():
    expected = ['1', '2', '3']
    input_string = '1,2,3,'
    obtained = list(cli.iter_split(input_string, ","))
    assert expected == obtained

    obtained = list(cli.iter_split(input_string, None))
    assert [input_string] == obtained

    input_string = '/path/path2/1.fasta /path/path2/2.fasta /path/path2/3.fasta'
    obtained = list(cli.iter_split(input_string, None))
    expected = [
        '/path/path2/1.fasta',
        '/path/path2/2.fasta',
        '/path/path2/3.fasta']
    assert expected == obtained


def test_calculate_chunksize():
    tota_jobs_todo = 100
    processes = 1
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_bam2fasta_valid_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')
        csv_path = os.path.join(location, "all_barcodes_meta.csv")
        barcodes_path = utils.get_test_data('10x-example/barcodes.tsv')
        renamer_path = utils.get_test_data('10x-example/barcodes_renamer.tsv')
        fastas_dir = os.path.join(location, "fastas")
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)
        parser = bam2fasta_args.create_parser()
        args = [
            '--filename', testdata1,
            '--count-valid-reads', '10',
            '--write-barcode-meta-csv', csv_path,
            '--barcodes', barcodes_path,
            '--rename-10x-barcodes', renamer_path,
            '--save-fastas', fastas_dir,
        ]
        parser.parse_args(args)


def test_run_bam2fasta():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')
        csv_path = os.path.join(location, "all_barcodes_meta.csv")
        barcodes_path = utils.get_test_data('10x-example/barcodes.tsv')
        renamer_path = utils.get_test_data('10x-example/barcodes_renamer.tsv')
        fastas_dir = os.path.join(location, "fastas")
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)

        status, out, err = utils.run_shell_cmd(
            'bam2fasta convert --count-valid-reads 10 --filename ' +
            testdata1 + ' --write-barcode-meta-csv ' + csv_path +
            ' --barcodes ' + barcodes_path + ' --rename-10x-barcodes ' +
            renamer_path + ' --save-fastas ' + fastas_dir,
            in_directory=location)


def test_collect_reduce_temp_fastas():
    tota_jobs_todo = 100
    processes = 1
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_unfiltered_umi_to_fasta():
    tota_jobs_todo = 100
    processes = 1
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_filtered_umi_to_fasta():
    tota_jobs_todo = 100
    processes = 1
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_bam_to_fasta():
    tota_jobs_todo = 100
    processes = 1
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = cli.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained
