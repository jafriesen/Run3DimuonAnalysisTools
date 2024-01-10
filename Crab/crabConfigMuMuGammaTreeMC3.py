import CRABClient
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'muMuGammaTree_Run3Summer22MiniAODv3_InclusiveDileptonMinBias_MINIAODSIM_19Dec2023'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muMuGammaTreeMC.py'

dir = "root://cmsxrootd.fnal.gov/"
file_list = open("fileList.txt", "r")
input_files = []
for f in file_list :
	input_files.append(dir+f)
config.Data.userInputFiles = input_files
config.Data.inputDBS = 'global'
#config.Data.outputDatasetTag = config.General.requestName
config.Data.outputDatasetTag = 'skimMuMuGammaTree_Run3Summer22MiniAODv3_InclusiveDileptonMinBias_MINIAODSIM'
config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-v1/RAW'

config.Debug.extraJDL = ['+REQUIRED_OS="rhel9"']

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'

# Where the output files will be transmitted to
config.Site.storageSite = 'T3_CH_CERNBOX'
