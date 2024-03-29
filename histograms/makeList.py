import gfal2

dirs = ["ParkingDoubleMuonLowMass0/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass0/","ParkingDoubleMuonLowMass1/crab_muMuGamma_15Aug2023_ParkingD\
oubleMuonLowMass1/","ParkingDoubleMuonLowMass2/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass2/","ParkingDoubleMuonLowMass3/crab_muMuGamma_15Aug2\
023_ParkingDoubleMuonLowMass3/","ParkingDoubleMuonLowMass4/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass4/","ParkingDoubleMuonLowMass5/crab_muMu\
Gamma_15Aug2023_ParkingDoubleMuonLowMass5/","ParkingDoubleMuonLowMass6/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass6/","ParkingDoubleMuonLowMas\
s7/crab_muMuGamma_15Aug2023_ParkingDoubleMuonLowMass7/"]
ntuples = []

for d in dirs :
    gfal2_ctx = gfal2.creat_context()
    dir_path = "davs://xrootd.cmsaf.mit.edu:1094/store/user/jfriesen/" + d
    print(dir_path)
    for subdir in gfal2_ctx.listdir(dir_path) :
        subdir_path = dir_path + "/" + subdir
        print(" "+subdir_path)
        for subsubdir in gfal2_ctx.listdir(subdir_path) :
            subsubdir_path = subdir_path + "/" + subsubdir
            print("     "+subsubdir_path)
            output = "/store/user/jfriesen/" + d + subdir + "/" +  subsubdir
            ntuples.extend([ (output + "/" + file) for file in gfal2_ctx.listdir(subsubdir_path) ])

print(ntuples)
file = open('muMuGammaTree_ntuples.txt','w')
file.writelines("\n".join(ntuples))
file.close()