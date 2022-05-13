#!/usr/bin/bash
mkdir ./ext
cd ./ext
git clone https://github.com/Christina-hshi/psirc.git
cd psirc
cd psirc-quant
# you may need to compile htslib under "ext/htslib" by following the README there ("make install" is optional and only possible with admin permissions)
cd ext/htslib/
autoheader
autoconf
./configure
make
make install
cd ..
cd ..
# make release
mkdir release
cd release
cmake ..
make psirc-quant
# the psirc-quant program can be found at "src/psirc-quant"
make install
