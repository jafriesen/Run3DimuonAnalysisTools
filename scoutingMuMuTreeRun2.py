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
        '/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/325/101/00000/1EDC0B34-D167-DB4E-81DA-EDDF7C3F0CF6.root'
 )
)

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scoutRun2.root")
)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun2_v2', '')

L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p5_SQ_dR_Max1p4","L1_DoubleMu0er2p0_SQ_dR_Max1p4","L1_DoubleMu4p5_SQ_OS","L1_TripleMu_5_3_3","L1_TripleMu_5_5_3","L1_QuadMu0","L1_SingleMu22"]

process.scoutingTree = cms.EDAnalyzer('ScoutingMuMuTreeMakerRun2',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Info),
                                      muons            = cms.InputTag("hltScoutingMuonPackerCalo"),
                                      pfjets           = cms.InputTag("hltScoutingCaloPacker"),
                                      tracks           = cms.InputTag("hltScoutingTrackPacker"),
                                      primaryVertices  = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
                                      displacedVertices  = cms.InputTag("hltScoutingMuonPackerCalo","displacedVtx"),
                                      rho         = cms.InputTag("hltScoutingCaloPacker","rho"),
                                  )

process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)