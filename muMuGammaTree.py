import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2022F/ParkingDoubleMuonLowMass0/MINIAOD/PromptReco-v1/000/360/390/00000/05607e2d-f8e5-41a7-8392-bb5d44c16a6d.root'
 )
)

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

#process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
#process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("mmgTree.root")
)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '')

#L1Info = ["L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6","L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6(5)","L1_DoubleMu0er1p5(4)_SQ_OS_dR_Max1p4"]
L1Info = ["L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4","L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6","L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6","L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu4_SQ_OS_dR_Max1p2"]
process.tree = cms.EDAnalyzer('MuMuGammaTreeMaker',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Info),
                                      muons            = cms.InputTag("slimmedMuons"),
                                      electrons        = cms.InputTag("slimmedElectrons"),
                                      photons          = cms.InputTag("slimmedPhotons"),
                                      pfcands          = cms.InputTag("packedPFCandidates"),
                                      primaryVertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                      displacedVertices  = cms.InputTag("slimmedSecondaryVertices"),
                                  )

#process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)
process.p = cms.Path(process.tree)
