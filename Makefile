# shamelessly adapt https://github.com/qiime2/q2-emperor/blob/master/Makefile
.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	flake8 mixamp

test: all
	pytest

test-install: all
	# ensure the package is installed and the app is buildable. this test
	# is a passive verification that non-py essential files are part of the
	# installed entity.
	cd /  # go somewhere to avoid relative imports
	python -c "import mixamp"

test-cov: all
	pytest --cov=mixamp

install: all
	$(PYTHON) setup.py install

dev: all
	pip install -e .

clean: distclean

distclean: ;