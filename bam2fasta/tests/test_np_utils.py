import numpy as np
from bam2fasta import np_utils


def test_memmap():

    siglist = [1, 2, 2, 3, 34]
    memmapped, filename = np_utils.to_memmap(np.array(siglist))
    # Assert that the data didn't change as a result of memory-mapping
    np.testing.assert_array_equal(memmapped, siglist)
    assert filename.endswith(".mmap")
