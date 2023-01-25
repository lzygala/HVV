#!/bin/sh -e

export X509_USER_PROXY=/afs/cern.ch/user/l/lzygala/.x509up_u115293
voms-proxy-info

source /cvmfs/sft.cern.ch/lcg/views/LCG_97/x86_64-centos7-gcc8-opt/setup.sh

OUTPUTFILE=$2;
INPUTFILE=$3;
YEAR=$4;
NAME=$1;
XS=$5
SW=$6

xrdcp --nopbar $INPUTFILE ./input_file.root
root -l -b -q  "nanoAODAnalyzer_SelectionTTreeMaker.C++(\"input_file.root\",\"$YEAR\",\"$NAME\",$XS,$SW)"

xrdcp -f --nopbar  TTree_tmp.root $OUTPUTFILE 

rm ./*
