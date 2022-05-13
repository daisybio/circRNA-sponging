#!/usr/bin/bash
mkdir ./ext
cd ./ext
git clone https://github.com/Christina-hshi/psirc.git
cd psirc/psirc-quant
# make release
mkdir release
cd release
cmake ..
make psirc-quant
# the psirc-quant program can be found at "src/psirc-quant"
make install
