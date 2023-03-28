import FWCore.ParameterSet.Config as cms

_baseDir = '/afs/cern.ch/work/a/amecca/Muons/CMSSW_12_0_2_patch1/src/RecoMuon/MuonSeedGenerator/test/production/ZToMuMu_M-50To120/AAAOK/'
_fileName = 'seedRebuild.root'

from subprocess import Popen, PIPE
from os import path
_p1 = Popen(['ls', _baseDir], stdout=PIPE)
_p2 = Popen(['grep', 'Chunk_'], stdin=_p1.stdout, stdout=PIPE)
_p1.stdout.close()
_files = [ 'file:'+path.join(_baseDir, chunk, _fileName) for chunk in  _p2.communicate()[0].decode('ASCII').rstrip('\n').split('\n') ]

# _files = _files[:10]

source = cms.Source("PoolSource",
                    fileNames=cms.untracked.vstring(_files)
)
