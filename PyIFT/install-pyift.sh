#!/usr/bin/env bash

python3 parser/parser.py && python3 setup.py develop --user && mv _pyift.*.so pyift/
