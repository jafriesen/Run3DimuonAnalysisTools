import CRABClient
from CRABClient.UserUtilities import config
config = config()

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'muMuGammaTree.py'

#config.Data.inputDataset = '/ParkingDoubleMuonLowMass1/Run2022F-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-v1/RAW'

# These values only make sense for processing data
#    Select input data based on a lumi mask
config.Data.lumiMask = 'Cert_Collisions2022_355100_362760_Golden.json'

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_US_MIT'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    datasets = [
        '/ParkingDoubleMuonLowMass0/Run2022G-PromptReco-v1/MINIAOD',
        '/ParkingDoubleMuonLowMass1/Run2022G-PromptReco-v1/MINIAOD',
        '/ParkingDoubleMuonLowMass2/Run2022G-PromptReco-v1/MINIAOD',
        '/ParkingDoubleMuonLowMass3/Run2022G-PromptReco-v1/MINIAOD',
        '/ParkingDoubleMuonLowMass4/Run2022G-PromptReco-v1/MINIAOD',
        '/ParkingDoubleMuonLowMass5/Run2022G-PromptReco-v1/MINIAOD',
        '/ParkingDoubleMuonLowMass6/Run2022G-PromptReco-v1/MINIAOD',
        '/ParkingDoubleMuonLowMass7/Run2022G-PromptReco-v1/MINIAOD',
    ]
    
    for dataset in datasets:
        config.Data.inputDataset = dataset
        config.General.requestName = "muMuGamma_3Oct2023_" + dataset.split('/')[1] + "_" + dataset.split('/')[2]
        crabCommand('submit', config = config)
