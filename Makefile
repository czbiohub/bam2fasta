PYTHON ?= python

all:
	$(PYTHON) setup.py build_ext -i

clean:
	$(PYTHON) setup.py clean --all
	cd doc && make clean

install: all
	$(PYTHON) setup.py install

dist: FORCE
	$(PYTHON) setup.py sdist

wheel:
	$(PYTHON) setup.py bdist_wheel

test: all
	pip install -e '.[test]'
	$(PYTHON) -m pytest

coverage: all
	$(PYTHON) setup.py clean --all
	$(PYTHON) setup.py build_ext -i
	$(PYTHON) -m pytest --cov=. --cov-report term-missing

FORCE: