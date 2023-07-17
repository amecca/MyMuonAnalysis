import FWCore.ParameterSet.Config as cms

from os import environ

process = cms.Process('MUONANALYSIS')

process.source = cms.Source('PoolSource', 
                            # fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/a/amecca/public/for_Daniele/seedRebuild_RECO.root'),
                            # fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/a/amecca/Muons/CMSSW_12_0_2_patch1/src/RecoMuon/MuonSeedGenerator/test/ZToMuMu_M-50To120_2files.root'),
                            fileNames = cms.untracked.vstring('file:{}/src/MyMuonAnalysis/DisplacedMuons/test/ZToMuMu.root'.format(environ['CMSSW_BASE'])),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') 
from Configuration.AlCa.GlobalTag import GlobalTag 
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_mc_FULL', '')

process.muonAnalysis = cms.EDAnalyzer('MuonGeneralAnalyzer', 
                                         verbose          = cms.untracked.bool(True),
                                         printOutput      = cms.untracked.bool(False),
                                         puInfoTag        = cms.InputTag('addPileupInfo'),
                                         gpCollectionTag  = cms.InputTag('genParticles'),
                                         muonTag          = cms.InputTag('muons'),
                                         oldSaTracksTag   = cms.InputTag('standAloneMuons::RECO'),
                                         oldSaUpdTracksTag= cms.InputTag('standAloneMuons:UpdatedAtVtx:RECO'),
                                         oldGlobalTracks  = cms.InputTag('globalMuons::RECO'),
                                         newSaTracksTag   = cms.InputTag('standAloneMuons::RESTA'),
                                         newSaUpdTracksTag= cms.InputTag('standAloneMuons:UpdatedAtVtx:RESTA'),
                                         newGlobalTracks  = cms.InputTag('globalMuons::RESTA'),
                                         genParticleMatch = cms.InputTag('muonSimClassifier:toPrimaries:RECO'),
                                         hitCounterParams = cms.untracked.PSet(
                                             debug = cms.untracked.bool(True),
                                             matchingFractionCut = cms.double(0.8)
                                         )
                                     )


process.p = cms.Path(
    process.muonAnalysis
)

# setattr(process,     'muonAnalysisOLD', muonAnalysisOLD)
# setattr(process,     'muonAnalysisNEW', muonAnalysisNEW)
# setattr(process, 'run_muonAnalysis', cms.Path(muonAnalysisOLD + muonAnalysisNEW))

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('muonAnalysis.root'),
                                   closeFileFast = cms.untracked.bool(False))

# process.load('MyMuonAnalysis.DisplacedMuons.files_ZToMuMu_cff')
print('INFO: # of files:', len(process.source.fileNames))
