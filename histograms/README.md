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
`fillDimuHistograms-binned.py` makes dimuon histograms in bins of pT, eta, and soft muon MVA. These commands submit 2000 jobs running `fillDimuHistograms-binned.py` with the 2022-24 ntuples listed in `muMuGammaTree_ntuples_22_23_24.txt`:

```
python3 submitFillHistogram.py -o h -e /eos/user/.../histos -t histos -n 2000 -s fillDimuHistograms-binned.py -l muMuGammaTree_ntuples_22_23_24.txt
```
`/eos/user/.../histos` is the output directory.

Add `-d` to set up a dry run or check that the script works by running a single job locally (out of 5000 here to limit the number of input files):
```
python3 fillDimuHistograms-16Dec2024.py -o h -n 5000 -j 1 -l muMuGammaTree_ntuples_22_23_24.txt
```
Run `hadd` to combine the outputs into one ROOT file:
```
hadd /eos/user/.../histos/output.root /eos/user/.../histos/h*
```
