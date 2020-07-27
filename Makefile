PYTHON ?= python

clean:
	$(PYTHON) setup.py clean --all

install:
	$(PYTHON) setup.py install

dist: FORCE
	$(PYTHON) setup.py sdist

test:
	$(PYTHON) setup.py install --user
	$(PYTHON) -m pytest

coverage:
	$(PYTHON) setup.py clean --all
	$(PYTHON) setup.py install --user
	$(PYTHON) -m pytest --cov=. --cov-report term-missing

FORCE: