import gfal2

input_dirs = { # label for organization : path to dataset
    ## Run2022C
    "ParkingDoubleMuonLowMass0-Run2022C-10Dec2022-v2-231004_024726-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass0/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass0/231004_024726/0000",
    "ParkingDoubleMuonLowMass1-Run2022C-10Dec2022-v3-231004_024729-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass1/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass1/231004_024729/0000",
    "ParkingDoubleMuonLowMass2-Run2022C-10Dec2022-v2-231004_024732-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass2/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass2/231004_024732/0000",
    "ParkingDoubleMuonLowMass3-Run2022C-10Dec2022-v2-231004_024735-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass3/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass3/231004_024735/0000",
    "ParkingDoubleMuonLowMass4-Run2022C-10Dec2022-v2-231004_024738-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass4/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass4/231004_024738/0000",
    "ParkingDoubleMuonLowMass5-Run2022C-10Dec2022-v2-231004_024741-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass5/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass5/231004_024741/0000",
    "ParkingDoubleMuonLowMass6-Run2022C-10Dec2022-v2-231004_024744-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass6/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass6/231004_024744/0000",
    "ParkingDoubleMuonLowMass7-Run2022C-10Dec2022-v2-231004_024747-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass7/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass7/231004_024747/0000",
    ## Run2022D
    "ParkingDoubleMuonLowMass0-Run2022D-10Dec2022-v2-231004_025429-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass0/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass0_Run2022D-10Dec2022-v2/231004_025429/0000",
    "ParkingDoubleMuonLowMass1-Run2022D-10Dec2022-v3-231004_025431-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass1/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass1_Run2022D-10Dec2022-v3/231004_025431/0000",
    "ParkingDoubleMuonLowMass2-Run2022D-10Dec2022-v2-231004_025434-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass2/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass2_Run2022D-10Dec2022-v2/231004_025434/0000",
    "ParkingDoubleMuonLowMass3-Run2022D-10Dec2022-v2-231004_025437-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass3/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass3_Run2022D-10Dec2022-v2/231004_025437/0000",
    "ParkingDoubleMuonLowMass4-Run2022D-10Dec2022-v2-231004_025440-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass4/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass4_Run2022D-10Dec2022-v2/231004_025440/0000",
    "ParkingDoubleMuonLowMass5-Run2022D-10Dec2022-v2-231004_025444-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass5/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass5_Run2022D-10Dec2022-v2/231004_025444/0000",
    "ParkingDoubleMuonLowMass6-Run2022D-10Dec2022-v2-231004_025447-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass6/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass6_Run2022D-10Dec2022-v2/231004_025447/0000",
    "ParkingDoubleMuonLowMass7-Run2022D-10Dec2022-v2-231004_025450-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass7/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass7_Run2022D-10Dec2022-v2/231004_025450/0000",
    ## Run2022D
    "ParkingDoubleMuonLowMass0-Run2022E-10Dec2022-v2-231004_030206-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass0/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass0_Run2022E-10Dec2022-v2/231004_030206/0000",
    "ParkingDoubleMuonLowMass1-Run2022E-10Dec2022-v2-231004_030209-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass1/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass1_Run2022E-10Dec2022-v2/231004_030209/0000",
    "ParkingDoubleMuonLowMass2-Run2022E-10Dec2022-v2-231004_030212-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass2/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass2_Run2022E-10Dec2022-v2/231004_030212/0000",
    "ParkingDoubleMuonLowMass3-Run2022E-10Dec2022-v2-231004_030215-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass3/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass3_Run2022E-10Dec2022-v2/231004_030215/0000",
    "ParkingDoubleMuonLowMass4-Run2022E-10Dec2022-v2-231004_030218-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass4/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass4_Run2022E-10Dec2022-v2/231004_030218/0000",
    "ParkingDoubleMuonLowMass5-Run2022E-10Dec2022-v2-231004_030221-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass5/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass5_Run2022E-10Dec2022-v2/231004_030221/0000",
    "ParkingDoubleMuonLowMass6-Run2022E-10Dec2022-v2-231004_030225-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass6/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass6_Run2022E-10Dec2022-v2/231004_030225/0000",
    "ParkingDoubleMuonLowMass7-Run2022E-10Dec2022-v2-231004_030227-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass7/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass7_Run2022E-10Dec2022-v2/231004_030227/0000",
    ## Run2022F
    "ParkingDoubleMuonLowMass0-Run2022F-PromptReco-v1-230815_203407-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass0/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass0/230815_203407/0000",
    "ParkingDoubleMuonLowMass1-Run2022F-PromptReco-v1-230815_203421-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass1/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass1/230815_203421/0000",
    "ParkingDoubleMuonLowMass2-Run2022F-PromptReco-v1-230815_203435-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass2/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass2/230815_203435/0000",
    "ParkingDoubleMuonLowMass3-Run2022F-PromptReco-v1-230820_224715-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass3/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass3/230820_224715/0000",
    "ParkingDoubleMuonLowMass4-Run2022F-PromptReco-v1-230815_203502-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass4/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass4/230815_203502/0000",
    "ParkingDoubleMuonLowMass5-Run2022F-PromptReco-v1-230815_203517-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass5/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass5/230815_203517/0000",
    "ParkingDoubleMuonLowMass5-Run2022F-PromptReco-v1-230815_203517-0001" : "/store/user/jfriesen/ParkingDoubleMuonLowMass5/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass5/230815_203517/0001",
    "ParkingDoubleMuonLowMass6-Run2022F-PromptReco-v1-230815_203531-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass6/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass6/230815_203531/0000",
    "ParkingDoubleMuonLowMass6-Run2022F-PromptReco-v1-230815_203531-0001" : "/store/user/jfriesen/ParkingDoubleMuonLowMass6/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass6/230815_203531/0001",
    "ParkingDoubleMuonLowMass7-Run2022F-PromptReco-v1-230815_203544-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass7/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass7/230815_203544/0000",
    "ParkingDoubleMuonLowMass7-Run2022F-PromptReco-v1-230815_203544-0001" : "/store/user/jfriesen/ParkingDoubleMuonLowMass7/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass7/230815_203544/0001",
    ## Run2022G
    "ParkingDoubleMuonLowMass0-Run2022G-PromptReco-v1-231004_031058-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass0/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass0_Run2022G-PromptReco-v1/231004_031058/0000",
    "ParkingDoubleMuonLowMass1-Run2022G-PromptReco-v1-231004_031101-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass1/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass1_Run2022G-PromptReco-v1/231004_031101/0000",
    "ParkingDoubleMuonLowMass2-Run2022G-PromptReco-v1-231004_031104-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass2/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass2_Run2022G-PromptReco-v1/231004_031104/0000",
    "ParkingDoubleMuonLowMass3-Run2022G-PromptReco-v1-231004_031106-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass3/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass3_Run2022G-PromptReco-v1/231004_031106/0000",
    "ParkingDoubleMuonLowMass4-Run2022G-PromptReco-v1-231004_031109-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass4/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass4_Run2022G-PromptReco-v1/231004_031109/0000",
    "ParkingDoubleMuonLowMass5-Run2022G-PromptReco-v1-231004_031112-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass5/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass5_Run2022G-PromptReco-v1/231004_031112/0000",
    "ParkingDoubleMuonLowMass6-Run2022G-PromptReco-v1-231004_031115-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass6/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass6_Run2022G-PromptReco-v1/231004_031115/0000",
    "ParkingDoubleMuonLowMass7-Run2022G-PromptReco-v1-231004_031117-0000" : "/store/user/jfriesen/ParkingDoubleMuonLowMass7/crab_muMuGamma_3Oct2023_ParkingDoubleMuonLowMass7_Run2022G-PromptReco-v1/231004_031117/0000",
    ## Run2023B
    "ParkingDoubleMuonLowMass0-Run2023B-PromptReco-v1-231001_173013-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass0/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass0/231001_173013/0000",
    "ParkingDoubleMuonLowMass1-Run2023B-PromptReco-v1-231001_173017-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass1/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass1/231001_173017/0000",
    "ParkingDoubleMuonLowMass2-Run2023B-PromptReco-v1-231001_173020-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass2/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass2/231001_173020/0000",
    "ParkingDoubleMuonLowMass3-Run2023B-PromptReco-v1-231001_173028-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass3/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass3/231001_173028/0000",
    "ParkingDoubleMuonLowMass4-Run2023B-PromptReco-v1-231001_173031-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass4/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass4/231001_173031/0000",
    "ParkingDoubleMuonLowMass5-Run2023B-PromptReco-v1-231001_173034-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass5/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass5/231001_173034/0000",
    "ParkingDoubleMuonLowMass6-Run2023B-PromptReco-v1-231001_173038-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass6/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass6/231001_173038/0000",
    "ParkingDoubleMuonLowMass7-Run2023B-PromptReco-v1-231001_173041-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass7/crab_muMuGamma_01Oct2023_ParkingDoubleMuonLowMass7/231001_173041/0000",
    ## Run2023C
    # v1
    "ParkingDoubleMuonLowMass0-Run2023C-PromptReco-v1-231001_174125-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass0/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass0/231001_174125/0000",
    "ParkingDoubleMuonLowMass1-Run2023C-PromptReco-v1-231001_174128-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass1/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass1/231001_174128/0000",
    "ParkingDoubleMuonLowMass2-Run2023C-PromptReco-v1-231001_174131-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass2/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass2/231001_174131/0000",
    "ParkingDoubleMuonLowMass3-Run2023C-PromptReco-v1-231001_174134-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass3/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass3/231001_174134/0000",
    "ParkingDoubleMuonLowMass4-Run2023C-PromptReco-v1-231001_174137-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass4/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass4/231001_174137/0000",
    "ParkingDoubleMuonLowMass5-Run2023C-PromptReco-v1-231001_174140-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass5/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass5/231001_174140/0000",
    "ParkingDoubleMuonLowMass6-Run2023C-PromptReco-v1-231001_174144-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass6/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass6/231001_174144/0000",
    "ParkingDoubleMuonLowMass7-Run2023C-PromptReco-v1-231001_174147-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass7/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass7/231001_174147/0000",
    # v2
    "ParkingDoubleMuonLowMass0-Run2023C-PromptReco-v2-231001_174337-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass0/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass0Run2023C-PromptReco-v2/231001_174337/0000",
    "ParkingDoubleMuonLowMass1-Run2023C-PromptReco-v2-231001_174340-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass1/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass1Run2023C-PromptReco-v2/231001_174340/0000",
    "ParkingDoubleMuonLowMass2-Run2023C-PromptReco-v2-231001_174343-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass2/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass2Run2023C-PromptReco-v2/231001_174343/0000",
    "ParkingDoubleMuonLowMass3-Run2023C-PromptReco-v2-231001_174346-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass3/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass3Run2023C-PromptReco-v2/231001_174346/0000",
    "ParkingDoubleMuonLowMass4-Run2023C-PromptReco-v2-231001_174349-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass4/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass4Run2023C-PromptReco-v2/231001_174349/0000",
    "ParkingDoubleMuonLowMass5-Run2023C-PromptReco-v2-231001_174352-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass5/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass5Run2023C-PromptReco-v2/231001_174352/0000",
    "ParkingDoubleMuonLowMass6-Run2023C-PromptReco-v2-231001_174355-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass6/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass6Run2023C-PromptReco-v2/231001_174355/0000",
    "ParkingDoubleMuonLowMass7-Run2023C-PromptReco-v2-231001_174358-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass7/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass7Run2023C-PromptReco-v2/231001_174358/0000",
    # v3
    "ParkingDoubleMuonLowMass0-Run2023C-PromptReco-v3-231001_174401-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass0/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass0Run2023C-PromptReco-v3/231001_174401/0000",
    "ParkingDoubleMuonLowMass1-Run2023C-PromptReco-v3-231001_174404-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass1/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass1Run2023C-PromptReco-v3/231001_174404/0000",
    "ParkingDoubleMuonLowMass2-Run2023C-PromptReco-v3-231001_174407-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass2/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass2Run2023C-PromptReco-v3/231001_174407/0000",
    "ParkingDoubleMuonLowMass3-Run2023C-PromptReco-v3-231001_174410-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass3/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass3Run2023C-PromptReco-v3/231001_174410/0000",
    "ParkingDoubleMuonLowMass4-Run2023C-PromptReco-v3-231001_174413-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass4/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass4Run2023C-PromptReco-v3/231001_174413/0000",
    "ParkingDoubleMuonLowMass5-Run2023C-PromptReco-v3-231001_174416-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass5/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass5Run2023C-PromptReco-v3/231001_174416/0000",
    "ParkingDoubleMuonLowMass6-Run2023C-PromptReco-v3-231001_174419-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass6/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass6Run2023C-PromptReco-v3/231001_174419/0000",
    "ParkingDoubleMuonLowMass7-Run2023C-PromptReco-v3-231001_174422-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass7/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass7Run2023C-PromptReco-v3/231001_174422/0000",
    # v4
    "ParkingDoubleMuonLowMass0-Run2023C-PromptReco-v4-231001_174425-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass0/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass0Run2023C-PromptReco-v4/231001_174425/0000",
    "ParkingDoubleMuonLowMass1-Run2023C-PromptReco-v4-231001_174428-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass1/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass1Run2023C-PromptReco-v4/231001_174428/0000",
    "ParkingDoubleMuonLowMass2-Run2023C-PromptReco-v4-231001_174431-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass2/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass2Run2023C-PromptReco-v4/231001_174431/0000",
    "ParkingDoubleMuonLowMass3-Run2023C-PromptReco-v4-231001_174434-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass3/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass3Run2023C-PromptReco-v4/231001_174434/0000",
    "ParkingDoubleMuonLowMass4-Run2023C-PromptReco-v4-231001_174437-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass4/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass4Run2023C-PromptReco-v4/231001_174437/0000",
    "ParkingDoubleMuonLowMass5-Run2023C-PromptReco-v4-231001_174440-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass5/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass5Run2023C-PromptReco-v4/231001_174440/0000",
    "ParkingDoubleMuonLowMass6-Run2023C-PromptReco-v4-231001_174443-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass6/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass6Run2023C-PromptReco-v4/231001_174443/0000",
    "ParkingDoubleMuonLowMass7-Run2023C-PromptReco-v4-231001_174446-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass7/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass7Run2023C-PromptReco-v4/231001_174446/0000",
    ## Run2023D
    # v1
    "ParkingDoubleMuonLowMass0-Run2023D-PromptReco-v1-231001_174639-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass0/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass0Run2023D-PromptReco-v1/231001_174639/0000",
    "ParkingDoubleMuonLowMass1-Run2023D-PromptReco-v1-231001_174642-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass1/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass1Run2023D-PromptReco-v1/231001_174642/0000",
    "ParkingDoubleMuonLowMass2-Run2023D-PromptReco-v1-231001_174645-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass2/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass2Run2023D-PromptReco-v1/231001_174645/0000",
    "ParkingDoubleMuonLowMass3-Run2023D-PromptReco-v1-231001_174648-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass3/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass3Run2023D-PromptReco-v1/231001_174648/0000",
    "ParkingDoubleMuonLowMass4-Run2023D-PromptReco-v1-231001_174651-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass4/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass4Run2023D-PromptReco-v1/231001_174651/0000",
    "ParkingDoubleMuonLowMass5-Run2023D-PromptReco-v1-231001_174655-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass5/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass5Run2023D-PromptReco-v1/231001_174655/0000",
    "ParkingDoubleMuonLowMass6-Run2023D-PromptReco-v1-231001_174658-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass6/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass6Run2023D-PromptReco-v1/231001_174658/0000",
    "ParkingDoubleMuonLowMass7-Run2023D-PromptReco-v1-231001_174701-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass7/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass7Run2023D-PromptReco-v1/231001_174701/0000",
    # v2
    "ParkingDoubleMuonLowMass0-Run2023D-PromptReco-v2-231001_174704-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass0/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass0Run2023D-PromptReco-v2/231001_174704/0000",
    "ParkingDoubleMuonLowMass1-Run2023D-PromptReco-v2-231001_174707-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass1/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass1Run2023D-PromptReco-v2/231001_174707/0000",
    "ParkingDoubleMuonLowMass2-Run2023D-PromptReco-v2-231001_174710-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass2/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass2Run2023D-PromptReco-v2/231001_174710/0000",
    "ParkingDoubleMuonLowMass3-Run2023D-PromptReco-v2-231001_174712-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass3/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass3Run2023D-PromptReco-v2/231001_174712/0000",
    "ParkingDoubleMuonLowMass4-Run2023D-PromptReco-v2-231001_174715-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass4/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass4Run2023D-PromptReco-v2/231001_174715/0000",
    "ParkingDoubleMuonLowMass5-Run2023D-PromptReco-v2-231001_174718-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass5/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass5Run2023D-PromptReco-v2/231001_174718/0000",
    "ParkingDoubleMuonLowMass6-Run2023D-PromptReco-v2-231001_174722-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass6/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass6Run2023D-PromptReco-v2/231001_174722/0000",
    "ParkingDoubleMuonLowMass7-Run2023D-PromptReco-v2-231001_174725-0000" : "/store/user/secholak/ParkingDoubleMuonLowMass7/crab_muMuGamma_C_01Oct2023_ParkingDoubleMuonLowMass7Run2023D-PromptReco-v2/231001_174725/0000",
}

ntuple_list = {}
nfiles = 0

for d in input_dirs :
    gfal2_ctx = gfal2.creat_context()
    dir_path = "davs://xrootd.cmsaf.mit.edu:1094/" + input_dirs[d]
    dir_ntuples = [ (input_dirs[d] + file) for file in gfal2_ctx.listdir(dir_path) ]
    nfiles += len(dir_ntuples)
    print(str(len(dir_ntuples)) + " files in " + dir_path)
    ntuple_list[d] = dir_ntuples

print("making list with " + str(nfiles) + " files")
with open("muMuGammaTree_ntuples_fullRun3.json", "w") as write_file:
    json.dump(ntuple_list, write_file)
