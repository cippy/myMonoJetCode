#!/bin/bash
# this was supposed to be used in Tier2, for lxplus it is better not to use bsub since the queue 8nh is always very busy
#
# in Tier2 one can do 
# bsub -q cmsan -J <jobname> -o <output-logfile> "source <script_name.sh>"
#
#
base_path="${CMSSW_BASE}/src/myMonoJetCode"       
cfg_path="${base_path}/config"
src_path="${base_path}/src"

cd $src_path

# add common options
commonOptions=""

# add specific options
options_sig="$commonOptions "
options_zmumu="$commonOptions "
options_zee="$commonOptions "
options_wmnu="$commonOptions "
options_wenu="$commonOptions "
options_gam="$commonOptions "

echo
echo
./tmp/main ${cfg_path}/monojet_signalRegion_config.txt $options_sig
echo
echo
./tmp/main ${cfg_path}/zmumujets_ControlRegion_config.txt $options_zmumu
echo
echo
./tmp/main ${cfg_path}/wmunujets_ControlRegion_config.txt $options_wmunu
echo
echo
./tmp/main ${cfg_path}/zeejets_ControlRegion_config.txt $options_zee
echo
echo
./tmp/main ${cfg_path}/wenujets_ControlRegion_config.txt $options_wenu
echo
echo
./tmp/main ${cfg_path}/gammajets_ControlRegion_config.txt $options_gam
echo
echo
echo "producing plots for monojet (exclusive)"
./../macro/distribution 0
echo
echo
echo "producing plots for monoV"
./../macro/distribution 1
echo
echo
echo "========================"
echo "===  End of Analysis ==="
echo "========================"
#echo "Checking for running jobs"
#bjobs
echo
