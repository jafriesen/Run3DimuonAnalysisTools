import CRABClient
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'muMuGammaTree_2022Test_MINIAOD_pt5to25_21Sep2023'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muMuGammaTree.py'

#config.Data.userInputFiles = ['root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(700)]+['root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_2/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(700)]+['root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_3/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(3500)]
config.Data.userInputFiles = ['root://cmsxrootd.fnal.gov//store/user/bgreenbe/EtaToMuMuGamma/Run3_2022_MINIAOD_3/EtaToMuMuGamma_2022Test_MINIAOD_%d.root'%i for i in range(3500)]
config.Data.inputDBS = 'global'
#config.Data.outputDatasetTag = config.General.requestName
config.Data.outputDatasetTag = 'skimEtaMuMuGamma_pythiaGun_MINIAOD_pt5to25'
config.Data.outputPrimaryDataset = 'CRAB_UserFiles'
config.Data.inputDBS = 'global'
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-v1/RAW'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_US_MIT'
