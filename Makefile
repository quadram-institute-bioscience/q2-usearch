.PHONY: all lint test install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8 --ignore=E501,W503

test: all
	py.test -vvv

install: all
	pip install .

dev: all
	pip install -e .

clean: distclean

distclean: ;
