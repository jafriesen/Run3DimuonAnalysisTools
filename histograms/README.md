How to make histograms from ntuples using HTCondor:

#### Setup and creating lists of ntuple paths
Make sure to correctly set up your proxy and run `cmsenv` before running Condor jobs:
```
voms-proxy-init --voms cms --valid 168:00 -out $HOME/private/.proxy
export X509_USER_PROXY=$HOME/private/.proxy
cmsenv
```
The plotting scripts here use lists of paths to ntuples. `makeList.py` is a crude script which uses gfal to make a list of files split into directories at a site like the T2 MIT, so modify the host and directories as necessary if you want to use it.

Running makeList.py:
```
eval `scram unsetenv -sh`
python makeList.py
cmsenv
```
`scram unsetenv -sh` may not be necessary but can help when accessing sites with gfal.

#### Run Condor jobs
`fillMMGHistograms.py` makes histograms of variables using combinations of input photon collections, photon selection criteria, and cuts on events. Make sure your list of ntuples is correct:
```
listDir = "/afs/cern.ch/user/.../CMSSW_12_4_2/src/Run3DimuonAnalysisTools/Plotting/FillHistogram"	
files = [ "root://cmsxrootd.fnal.gov/"+line for line in open(listDir+"/muMuGammaTree_ntuples.txt")]
```
These commands submit 2000 jobs running `fillMMGHistograms.py` on the ntuples listed in `muMuGammaTree_ntuples.txt`:
```
mkdir /eos/user/.../histos_etaToMuMuGamma_mass_minDr
python3 submitFillHistogram.py -o h_misc -e /eos/user/.../histos_etaToMuMuGamma_mass_minDr -t histos_etaToMuMuGamma_mass_minDr -n 2000 -s fillMMGHistograms.py -l muMuGammaTree_ntuples.txt
```
`/eos/user/.../histos_etaToMuMuGamma_mass_minDr` is the output directory.

Add `-d` to set up a dry run or check that the script works by running a single job locally (out of 5000 here to limit the number of input files):
```
python3 /afs/cern.ch/user/.../CMSSW_12_4_2/src/Run3DimuonAnalysisTools/Plotting/FillHistogram/fillMMGHistograms.py -i None -o h_mass -n 5000 -j 5
```
Run `hadd` to combine the outputs into one ROOT file:
```
hadd /eos/user/.../histos_etaToMuMuGamma_mass_minDr/output.root /eos/user/.../histos_etaToMuMuGamma_mass_minDr/h_mass*
```
