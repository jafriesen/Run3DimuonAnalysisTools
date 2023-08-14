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
        #'/store/data/Run2018D/ScoutingCaloMuon/RAW/v1/000/320/570/00000/F2AC3A45-9494-E811-9CD1-FA163E390D83.root',
        #'/store/data/Run2022F/ScoutingPFRun3/RAW/v1/000/362/167/00000/bfafc181-f546-4f07-8a99-d56c3890e4a5.root',
        #'/store/data/Run2022F/ScoutingPFRun3/RAW/v1/000/361/303/00000/37f1ab1d-94f9-4177-91e5-db46490bc69a.root',
        '/store/data/Run2022F/ScoutingPFRun3/RAW/v1/000/362/091/00000/839020c4-51a8-4c21-a6d0-5b73bf0cbe14.root'
        #'/store/data/Run2022F/ScoutingPFRun3/RAW/v1/000/362/091/00000/82a8207c-b475-4f93-9290-6ef1e84b42c4.root'
 )
)

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("scout.root")
)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '')

L1Info = ["L1_DoubleMu_12_5","L1_DoubleMu_15_7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7","L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18","L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p5_SQ_dR_Max1p4","L1_DoubleMu0er2p0_SQ_dR_Max1p4","L1_DoubleMu4p5_SQ_OS","L1_TripleMu_5_3_3","L1_TripleMu_5_5_3","L1_QuadMu0","L1_SingleMu22"]

process.scoutingTree = cms.EDAnalyzer('ScoutingMuMuGammaTreeMakerRun3',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Info),
                                      muons            = cms.InputTag("hltScoutingMuonPacker"),
                                      electrons        = cms.InputTag("hltScoutingEgammaPacker"),
                                      photons          = cms.InputTag("hltScoutingEgammaPacker"),
                                      pfcands          = cms.InputTag("hltScoutingPFPacker"),
                                      pfjets           = cms.InputTag("hltScoutingPFPacker"),
                                      tracks           = cms.InputTag("hltScoutingTrackPacker"),
                                      primaryVertices  = cms.InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"),
                                      displacedVertices  = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
                                      pfMet            = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
                                      pfMetPhi         = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
                                      rho         = cms.InputTag("hltScoutingPFPacker","rho"),
                                  )

process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)