import gfal2

dirs = ["ParkingDoubleMuonLowMass0/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass0/","ParkingDoubleMuonLowMass1/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass1/","ParkingDoubleMuonLowMass2/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass2/","ParkingDoubleMuonLowMass3/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass3/","ParkingDoubleMuonLowMass4/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass4/","ParkingDoubleMuonLowMass5/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass5/","ParkingDoubleMuonLowMass6/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass6/","ParkingDoubleMuonLowMass7/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass7/"]
ntuples = []

for d in dirs :
    gfal2_ctx = gfal2.creat_context()
    mainDir = "davs://xrootd.cmsaf.mit.edu:1094/store/user/jfriesen/" + d
    subDir1 = gfal2_ctx.listdir(mainDir)[0] + "/"
    mainDir = mainDir + subDir1
    subDir2 = gfal2_ctx.listdir(mainDir)[0] + "/"
    mainDir = mainDir + subDir2
    inputDir = "/store/user/jfriesen/" + d + subDir1 + subDir2
    print(inputDir)
    print(mainDir)
    ntuples.extend([ (inputDir + "/" + file) for file in gfal2_ctx.listdir(mainDir) ])

print(ntuples)
file = open('muMuGammaTree_ntuples.txt','w')
file.writelines("\n".join(ntuples))
file.close()
