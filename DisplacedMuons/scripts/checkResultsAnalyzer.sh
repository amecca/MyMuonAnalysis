#!/bin/sh

rootmacro=$CMSSW_BASE/src/MyMuonAnalysis/DisplacedMuons/macro/checkResults_MuonGeneralAnalyzer.C
if [ $# -ge 1 ] ; then
    root -l -q $rootmacro\(\"$1\"\)
else
    root -l -q $rootmacro
fi
