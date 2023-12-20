import CRABClient
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'muMuGammaTree_Run3Summer22MiniAODv3_InclusiveDileptonMinBias_MINIAODSIM_19Dec2023'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muMuGammaTree.py'

config.Data.inputDataset = "/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/Run3Summer22MiniAODv3-Pilot_124X_mcRun3_2022_realistic_v12-v5/MINIAODSIM"
config.Data.inputDBS = 'global'
#config.Data.outputDatasetTag = config.General.requestName
#config.Data.outputDatasetTag = 'skimEtaMuMuGamma_InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8_MINIAODSIM'
#config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-v1/RAW'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'

# Where the output files will be transmitted to
config.Site.storageSite = 'T3_CH_CERNBOX'
