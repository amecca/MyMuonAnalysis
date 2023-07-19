#!/bin/sh

set -o pipefail

cmsRun recoMuonGeneralAnalyzer_cfg.py -n 1 2>&1 | gzip --fast > cmsrun.log.gz
exitStatus=$?
bell-ring
exit $exitStatus
