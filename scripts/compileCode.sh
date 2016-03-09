#!/bin/bash
# ok, it's working
#launch as --> source <script_name>.sh
#
#be sure that environment is set
cd $CMSSW_BASE/src
cmsenv
echo
cd myMonoJetCode/src
#now compile main
echo "Compiling code with make command"
make
echo
# prepare to compile distribution.C macro
cd ..
MACRO_FOLDER="macro"
MACRO_NAME="distribution.C"
EXE_NAME="distribution"
echo "Compiling $MACRO_FOLDER/$MACRO_NAME"
g++ -Wall -pedantic -lm -o $MACRO_FOLDER/$EXE_NAME $MACRO_FOLDER/$MACRO_NAME `rootlib`
echo
echo "Now you are in folder:"
echo $PWD
echo
