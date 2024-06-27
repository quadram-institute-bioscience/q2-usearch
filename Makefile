.PHONY: all lint test install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test -vvv /data/qiime-dev/q2-usearch/q2_usearch/tests/test_fastx_truncate.py::TestFastxTruncate::test_fastx_truncate_custom_params

install: all
	pip install .

dev: all
	pip install -e .

clean: distclean

distclean: ;
