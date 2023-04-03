import FWCore.ParameterSet.Config as cms

from os import environ

process = cms.Process('MUONANALYSIS')

process.source = cms.Source('PoolSource', 
                            # fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/a/amecca/public/for_Daniele/seedRebuild_RECO.root'),
                            # fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/amecca/Muons/CMSSW_12_0_2_patch1/src/RecoMuon/MuonSeedGenerator/test/ZToMuMu_M-50To120_2files.root'),
                            fileNames = cms.untracked.vstring('file:{}/src/MyMuonAnalysis/DisplacedMuons/test/ZToMuMu.root'.format(environ['CMSSW_BASE'])),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5))


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') 
from Configuration.AlCa.GlobalTag import GlobalTag 
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_mc_FULL', '')

process.muonAnalysisOLD = cms.EDAnalyzer('MuonGeneralAnalyzer', 
                                         verbose          = cms.untracked.bool(True),
                                         printOutput      = cms.untracked.bool(False),
                                         puInfoTag        = cms.InputTag('addPileupInfo'),
                                         gpCollectionTag  = cms.InputTag('genParticles'),
                                         muonTag          = cms.InputTag('muons'),
                                         saTracksTag      = cms.InputTag('standAloneMuons::RECO'),
                                         saUpdTracksTag   = cms.InputTag('standAloneMuons:UpdatedAtVtx:RECO'),
                                         globalTracks     = cms.InputTag('globalMuons::RECO'),
                                         genParticleMatch = cms.InputTag('muonSimClassifier:toPrimaries:RECO')
                                     )

process.muonAnalysisNEW = cms.EDAnalyzer('MuonGeneralAnalyzer', 
                                         verbose          = cms.untracked.bool(True),
                                         printOutput      = cms.untracked.bool(False),
                                         puInfoTag        = cms.InputTag('addPileupInfo'),
                                         gpCollectionTag  = cms.InputTag('genParticles'),
                                         muonTag          = cms.InputTag('muons'),
                                         saTracksTag      = cms.InputTag('standAloneMuons::RESTA'),
                                         saUpdTracksTag   = cms.InputTag('standAloneMuons:UpdatedAtVtx:RESTA'),
                                         globalTracks     = cms.InputTag('globalMuons::RESTA')
                                     )

process.p = cms.Path(
    process.muonAnalysisOLD
    + process.muonAnalysisNEW
)

# setattr(process,     'muonAnalysisOLD', muonAnalysisOLD)
# setattr(process,     'muonAnalysisNEW', muonAnalysisNEW)
# setattr(process, 'run_muonAnalysis', cms.Path(muonAnalysisOLD + muonAnalysisNEW))

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('muonAnalysis.root'),
                                   closeFileFast = cms.untracked.bool(False))

# process.load('MyMuonAnalysis.DisplacedMuons.files_ZToMuMu_cff')
print('INFO: # of files:', len(process.source.fileNames))
