import os
import shutil

import screed

from bam2fasta.tests import bam2fasta_tst_utils as utils
from bam2fasta import cli
from bam2fasta import VERSION


def test_iter_split():
    expected = ['1', '2', '3']
    input_string = '1,2,3,'
    obtained = list(cli.iter_split(input_string, ","))
    assert expected == obtained

    obtained = list(cli.iter_split(input_string, None))
    assert [input_string] == obtained

    input_string = \
        '/path/path2/1.fasta /path/path2/2.fasta /path/path2/3.fasta'
    obtained = list(cli.iter_split(input_string, None))
    expected = [
        '/path/path2/1.fasta',
        '/path/path2/2.fasta',
        '/path/path2/3.fasta']
    assert expected == obtained


def test_bam2fasta_info():
    status, out, err = utils.run_shell_cmd('bam2fasta info')

    # no output to stderr
    assert not err
    print(status, out, err)
    assert "bam2fasta version" in out
    assert "loaded from path" in out
    assert VERSION in out


def test_bam2fasta_info_verbose():
    status, out, err = utils.run_shell_cmd('bam2fasta info -v')

    # no output to stderr
    assert not err
    assert "pathos version" in out
    assert "screed version" in out
    assert "pysam version" in out
    assert "loaded from path" in out


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


def test_run_bam2fasta_supply_all_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')
        csv_path = os.path.join(location, "all_barcodes_meta.csv")
        barcodes_path = utils.get_test_data('10x-example/barcodes.tsv')
        renamer_path = utils.get_test_data('10x-example/barcodes_renamer.tsv')
        fastas_dir = os.path.join(location, "fastas")
        temp_fastas_dir = os.path.join(
            os.path.dirname(testdata1), "temp_fastas/")
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)
        if not os.path.exists(temp_fastas_dir):
            os.makedirs(temp_fastas_dir)

        status, out, err = utils.run_shell_cmd(
            'bam2fasta convert --filename ' + testdata1 +
            ' --min-umi-per-barcode 10' +
            ' --write-barcode-meta-csv ' + csv_path +
            ' --save-intermediate-files ' + temp_fastas_dir +
            ' --barcodes ' + barcodes_path + ' --rename-10x-barcodes ' +
            renamer_path + ' --save-fastas ' + fastas_dir + " --processes 1",
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
        shutil.rmtree(temp_fastas_dir)


def test_run_bam2fasta_default_args():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')
        fastas_dir = os.path.join(location)
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)

        status, out, err = utils.run_shell_cmd(
            'bam2fasta convert --filename ' + testdata1,
            in_directory=location)

        assert status == 0
        fasta_files = os.listdir(fastas_dir)
        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 8


def test_run_convert():
    with utils.TempDirectory() as location:
        testdata1 = utils.get_test_data('10x-example/possorted_genome_bam.bam')
        fastas_dir = os.path.join(location)
        if not os.path.exists(fastas_dir):
            os.makedirs(fastas_dir)

        fasta_files = cli.convert(
            ['--filename', testdata1, '--save-fastas', location])

        barcodes = [filename.replace(".fasta", "") for filename in fasta_files]
        assert len(barcodes) == 8
