from pkg_resources import Requirement, resource_filename, ResolutionError
import os
import tempfile
import shutil


class TempDirectory(object):
    def __init__(self):
        self.tempdir = tempfile.mkdtemp(prefix='sourmashtest_')

    def __enter__(self):
        return self.tempdir

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            shutil.rmtree(self.tempdir, ignore_errors=True)
        except OSError:
            pass

        if exc_type:
            return False


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(
            Requirement.parse(
                "sourmash"), "sourmash/sourmash/test-data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'test-data',
                                filename)
    return filepath
