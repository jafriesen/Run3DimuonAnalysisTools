#### Run on all ParkingDoubleMuonLowMass datasets
Start parallel CRAB jobs on Run2022F/ParkingDoubleMuonLowMass0-7 datasets:
(Modify datasets and job name in crabConfigMuMuGammaTree.py as necessary)
```
cmsenv
source /cvmfs/cms.cern.ch/common/crab-setup.sh
voms-proxy-init --voms cms
python3 crabConfigMuMuGammaTree.py
```
