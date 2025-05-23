## Single-locus Species Delimitation using ASAP
## Date: 08 May 2025
## Author: C.E. Supsup

## Assemble Species by Automatic Partitioning (ASAP; Puillandre et al. 2021)
## Input - aligned sequence in .fas file or pairwise genetic distances in .csv format
## Analysis implemented in python or via web server (see https://bioinfo.mnhn.fr/abi/public/asap/)

## Required python modules
pip3 install PyQtWebEngine
pip3 install PyQt5

## Download ASAPy from https://github.com/iTaxoTools/ASAPy (Vences et al. 2021)
## Unzip file and navigate to ASAPy folder to install
cd /Users/christiansupsup/ASAPy-main

pip3 install -r requirements.txt
python3 setup.py build_ext --inplace
python3 setup.py build_qt

## Launch
python3 launcher.py

## Use Mega (https://www.megasoftware.net/) to generate the uncorrected pairwise genetic distances
## use genetic distances in .csv format as an input
## in ASAPy software, adjust the following settings
##	check "MEGA CSV" option
##	change the "Sequence Length" to 876 (bp) - the sequence length of Bronchocela mtDNA data
##	change the "Methond" to "Simple Distance"
