import numpy as np
import tempfile


def to_memmap(array):
    """Return a memory-map to an array stored in a binary file on disk.
    Memory-mapped files are used for accessing small segments of
    large files on disk, without reading the entire file into memory.

    Parameters
    ----------
    array : np.array
        numpy array to memory map
    Returns
    -------
        memmap_array: np.array
            memmap memory mapped array
        filename: str
            name of the file that memory mapped array is written to
    """
    filename = tempfile.NamedTemporaryFile(
        prefix="array", suffix=".mmap", delete=False).name
    shape = array.shape
    f = np.memmap(filename, mode='w+', shape=shape, dtype=array.dtype)
    f[:] = array[:]
    del f
    memmap_array = np.memmap(filename, dtype=array.dtype, shape=shape)
    return memmap_array, filename
