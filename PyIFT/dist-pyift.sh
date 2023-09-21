#!/usr/bin/env bash

python3 setup.py build_ext && python3 setup.py build_py && python3 setup.py bdist_wheel
