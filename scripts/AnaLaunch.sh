#!/bin/bash
# this was supposed to be used in Tier2, for lxplus it is better not to use bsub since the queue 8nh is always very busy
#
# in Tier2 one can do 
# bsub -q cmsan -J <jobname> -o <output-logfile> "source <script_name.sh>"
#
#
echo
#./tmp/main config/monojet_signalRegion_config_spring15_25ns.txt
pwd
ls
echo
echo
#./tmp/main config/zmumujets_ControlRegion_config_spring15_25ns.txt
echo
echo
#./tmp/main config/wmunujets_ControlRegion_config_spring15_25ns.txt
echo
echo
#./tmp/main config/zeejets_ControlRegion_config_spring15_25ns.txt
#echo
#echo
#./tmp/main config/zeejets_ControlRegion_config_spring15_25ns.txt -calibEle
#echo
#echo
#./tmp/main config/wenujets_ControlRegion_config_spring15_25ns.txt
#echo
#echo
#./tmp/main config/wenujets_ControlRegion_config_spring15_25ns.txt -calibEle
#echo
#echo
#echo "producing plots for monojet (exclusive)"
#./macro/distribution 0
echo
echo
#echo "producing plots for monoV"
#./macro/distribution 1
echo
echo
echo "========================"
echo "===  End of Analysis ==="
echo "========================"
echo "Checking for running jobs"
bjobs
echo
