#!/bin/bash

base_path="${CMSSW_BASE}/src/myMonoJetCode"       
cfg_path="${base_path}/config"
src_path="${base_path}/src"

cd $src_path
options="-nw"

./tmp/main ${cfg_path}/synch_80X_signal.txt $options
./tmp/main ${cfg_path}/synch_80X_zmumu.txt $options
./tmp/main ${cfg_path}/synch_80X_zee.txt $options
./tmp/main ${cfg_path}/synch_80X_wmunu.txt $options
./tmp/main ${cfg_path}/synch_80X_wenu.txt $options
./tmp/main ${cfg_path}/synch_80X_gamma.txt $options
