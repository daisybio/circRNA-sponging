#!/usr/bin/bash
mkdir ./ext/PITA
cd ./ext/PITA

wget --no-check-certificate "https://genie.weizmann.ac.il/pubs/mir07/64bit_exe_pita_prediction.tar.gz"
tar xvfz *_pita_prediction.tar.gz
make install
# make script compatible with newer perl versions
sed -i -E "s/(=~\s\S+)\{HOME\}(.\S+)/\1\\\{HOME\\\}\2/" ./lib/libfile.pl
sed -i -E "s/(if\()defined\((@\S+)\)(.*)/\1\2\3/" ./lib/join.pl
