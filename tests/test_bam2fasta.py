from . import bam2fasta_tst_utils as utils
from bam2fasta import cli
import os
import screed


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

        assert status == 0
        with open(csv_path, 'rb') as f:
            data = [line.split() for line in f]
        assert len(data) == 9
        fasta_files = os.listdir(fastas_dir)
        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 1
        assert len(fasta_files) == 1
        assert barcodes[0] == 'lung_epithelial_cell|AAATGCCCAAACTGCT-1'
        count = 0
        fasta_file_name = os.path.join(fastas_dir, fasta_files[0])
        for record in screed.open(fasta_file_name):
            name = record.name
            sequence = record.sequence
            count += 1
            assert name.startswith('lung_epithelial_cell|AAATGCCCAAACTGCT-1')
            assert sequence.count(">") == 0
            assert sequence.count("X") == 0


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
