from . import bam2fasta_test_utils as utils
from bam2fasta import bam2fasta
import unittest


def test_iter_split():
    expected = ['1', '2', '3']
    input_string = '1,2,3,'
    obtained = list(bam2fasta.iter_split(input_string, ","))
    assert expected == obtained

    obtained = list(bam2fasta.iter_split(input_string, None))
    assert [input_string] == obtained

    input_string = '/path/path2/1.fasta /path/path2/2.fasta /path/path2/3.fasta'
    obtained = list(bam2fasta.iter_split(input_string, None))
    expected = [
        '/path/path2/1.fasta',
        '/path/path2/2.fasta',
        '/path/path2/3.fasta']
    assert expected == obtained


def test_calculate_chunksize():
    tota_jobs_todo = 100
    processes = 1
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_collect_reduce_temp_fastas():
    tota_jobs_todo = 100
    processes = 1
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_unfiltered_umi_to_fasta():
    tota_jobs_todo = 100
    processes = 1
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_filtered_umi_to_fasta():
    tota_jobs_todo = 100
    processes = 1
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained


def test_bam_to_fasta():
    tota_jobs_todo = 100
    processes = 1
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert tota_jobs_todo == obtained
    tota_jobs_todo = 51
    processes = 5
    expected = 11
    obtained = bam2fasta.calculate_chunksize(tota_jobs_todo, processes)
    assert expected == obtained
