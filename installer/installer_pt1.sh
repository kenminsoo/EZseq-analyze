#!/bin/bash

#this is the first part of the installer
#it will install a python package required to download larger files
#this requires python3
#run this first

python3 -m pip install --user pipx
python3 -m pipx ensurepath
