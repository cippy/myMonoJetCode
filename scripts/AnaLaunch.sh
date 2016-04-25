#!/bin/bash
# this was supposed to be used in Tier2, for lxplus it is better not to use bsub since the queue 8nh is always very busy
#
# in Tier2 one can do 
# bsub -q cmsan -J <jobname> -o <output-logfile> "source <script_name.sh>"
#
#
MY_PATH=$CMSSW_BASE/src/myMonoJetCode
MY_CFG_PATH=$MY_PATH/config
MY_MACRO_PATH=$MY_PATH/macro
cd $MY_PATH/src
echo
echo
./tmp/main $MY_CFG_PATH/monojet_signalRegion_config.txt
echo
echo
./tmp/main $MY_CFG_PATH/zmumujets_ControlRegion_config.txt
echo
echo
./tmp/main $MY_CFG_PATH/wmunujets_ControlRegion_config.txt
echo
echo
./tmp/main $MY_CFG_PATH/zeejets_ControlRegion_config.txt
echo
echo
./tmp/main $MY_CFG_PATH/wenujets_ControlRegion_config.txt
echo
echo
#echo "producing plots for monojet (exclusive)"
./../macro/distribution 0
echo
echo
#echo "producing plots for monoV"
./../macro/distribution 1
echo
echo
echo "========================"
echo "===  End of Analysis ==="
echo "========================"
#echo "Checking for running jobs"
#bjobs
echo
