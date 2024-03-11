import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000000)
)

dir = "root://xrootd-cms.infn.it//"

n = [ 24, 26, 27, 30, 31, 34, 35, 40 ]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        ["file:/eos/user/a/asterenb/EtaToGammaApr_AprToMuMu_2022/EtaToGammaApr_AprToMuMu_2022_MINIAOD_" + str(int(i)) + ".root" for i in range(100)]
    )
)

'''
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/f365faac-9a86-45cc-b1b4-87f036b1a2fa.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/fad35f5f-e549-40d0-a48f-d2b5e2086041.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/c95208fb-fa54-40a4-a000-5a65877fb3c2.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/7019e6c1-2adb-4fb0-b04e-f97449a0dfa4.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/9d7dd0b0-a52b-41a7-8197-fe0f7d051af2.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/678a8d19-b4d8-48be-89d6-c98bca514b69.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/839d6b3e-e85b-4c44-9c64-8632971e300f.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/1acebde6-6bee-464f-8775-2b161f1380ec.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/5351effc-3ca0-44f0-9acc-e1e355860c43.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/23aa49a6-a411-4e54-9212-f5ec90c0cb39.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/163eb519-49a0-4dc0-9aad-55a53add0a5e.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/209cd86b-b3a3-4d8f-88ff-1916399e27a7.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/f2432d3e-b8c2-48cd-a206-d0bd38507d4b.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/4c114f30-acb4-44d6-9390-34991c4d1f5a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/a3285029-4cde-42f2-9f87-71be739778a8.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/48c83db8-b0f1-4343-b2e4-e85edd5c6fba.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/0502f452-0df8-4d85-890d-443e015e31eb.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/c7590f38-fbbf-473f-ac12-6a39ae090b79.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/e8c46d5e-ba1f-411a-a0c4-ad8792f34a10.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/033d3edc-1ccd-4817-889b-42b051cbcd93.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/20db678c-8d7c-4178-9d4f-f1be0a1cea5a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/f60926d1-ac49-4e2d-82af-413131c51eec.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/c5c28304-6cdf-4826-9044-390367934329.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/52767c90-1518-4f68-8086-871f192b3f90.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/5e63e12c-e91e-47cc-8be5-7c405f59fccf.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/6a766825-17c2-42bc-9043-0aaa3efc28f8.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/10fde694-e715-40d1-a62a-b739e9ee8ffd.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/46c4b7ee-9fd7-4b21-90a3-68040787e422.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/6ef49640-a2a7-4dc6-b436-6e150c4f8b7f.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/6d0f448e-bf7a-42e8-940b-7f4b49133c3a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/6379e704-dc95-4cef-8f31-c67025fc2283.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/88edbf06-93c0-4b19-92e5-023daf632497.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/53bb35c7-b08a-4db7-bc3d-a3d3932179fd.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/d76f4623-eef5-4464-82ed-80c55f9e57a5.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/ac3ba788-8455-46f2-bf25-9e5b881c5b80.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/c461521f-6921-447e-be15-7453316167f1.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/d80dbb75-ddbd-4fcf-96e7-b64906f7aab1.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/e41cff72-5b92-4389-8f9b-2e15c0f1a5a6.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/b6a40e02-e656-4e73-9f50-4900b4ba451e.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/5669bc8e-e7ed-49fa-9fbd-995c363c05b3.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/325a76ee-d933-4f7e-be91-69f21397254c.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/9bf181f9-87f9-4202-90e5-ff7c9ff9097a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/5849399c-009e-45e1-9f08-6e85f5c97f33.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/33cf0a7f-9619-42f8-9de2-77dd9013f64e.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/f57c96d3-dec6-4f10-8a02-25e8729abe73.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/2db97211-1e3e-4e45-bbb7-a7e28f18e8bb.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/bae8d381-189a-42cc-a0fe-a26f2eb5584c.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/3e5bd985-dfd3-4401-982b-59c58942ae54.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/60a9c65a-3256-4a6a-97be-f46ee0fbe46b.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/6dbacede-e601-4a2c-8e8c-870b57449a70.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/233b52e0-fc2f-4f8b-b578-d5f46ce50b25.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/27de452a-9cff-45df-98ef-5c7ac814a8e4.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/a8b185c0-c490-4e87-bccf-31a294e8566d.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/dae5608b-81c2-427e-96ca-47b9c28b8b13.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/04fd9007-b793-4f02-aba9-a7a077a7830f.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/acadb9c5-65a9-4844-9837-cb080550f22a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/ff476702-f226-4f6a-90b6-17b6e6ed4c1a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/67679309-1416-4d2f-a5e8-09134e0cc5c7.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/6961ed1a-6077-42b3-a553-3e17ba849f2b.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/95ecc754-ecf3-44c4-91ec-804f8537f9c2.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/70000/44d77e2b-45cc-4a37-b14c-2d77e033b627.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/6228acfb-be8b-4c99-87e3-3529acc18b2e.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/e23a60c3-816f-461a-b268-92a5b9c7f0e0.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/68a45379-aa93-4d4f-ad5f-0aedeb0d1faf.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/c19aee31-1b45-45cd-abe6-0162cee20f29.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/5181c86d-aa5d-46d1-8163-89dbf27675cb.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/13cff5cc-e7e1-4ace-b5d5-05332c2400a4.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/5ce51d56-7cb7-49c3-866c-2c85f1cf7146.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/c08e9e9f-740b-496e-8ba4-0cd18fc0bc9a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/3999175d-90c9-4d9d-b7fa-1eeaea80b376.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/b11d36b0-defd-41b5-b432-b7d94bc83431.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/2b5b0696-5cfa-4911-bc19-dbee1272ac51.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/ef85ed89-9d92-47e9-99a2-936d47f57b85.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/14fd788e-2f1a-4a1b-a416-5673e6fdf76a.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/a7eb947d-6046-43f7-a692-31b2b19607e3.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/0bd97657-7ab9-466b-a2b4-51710347f798.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/9c223eff-f409-4b1a-9522-96b3e8206618.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/6f4fca18-2094-4c25-ac5a-78af2c5aaaea.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/777fa73f-e977-45b8-89ff-f0df4d963989.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/febf84d4-a26d-4eb2-9b7c-43116ba80005.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/e3c033fa-adc8-4eab-b21d-3326a63ff8d0.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/08460517-c8a6-460e-9fb8-dd82c3a310fe.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/8026708e-d782-4b53-b9e4-9741ac325734.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/2fc5d207-6f26-4e6e-af87-308cf4e05924.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/01f15ab6-678b-4765-a4f7-d909af09a809.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2830000/e3888142-dae1-44c5-8863-07e93749310e.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/383afc71-df31-4d93-bc8d-4c03bcdccc65.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/cf204d04-cd5b-40c0-8793-9ce235c3fdd7.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/60a990c2-3b65-4ded-a33e-6bf672c452d2.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/64dc669b-8255-4066-8f1d-830ae1ebca6f.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/c23ca7f3-1d46-4b52-9890-c3a3a03feb19.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/b52fa251-b850-420d-b000-fe1b0532359c.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/ef070ba2-d9ca-47e5-865c-d9453d1c6065.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/ea9a05c4-056c-479a-b4fd-83d6f93c9b06.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/826bb022-dda9-4f05-92c1-44077ccda396.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/3f190968-8eb9-45f1-903e-d67946236b52.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/eff23f70-00b8-4b1e-bdd7-ed77fb32fc78.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/e627de21-f30b-4876-8141-13d713f2f72c.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/4da60109-d559-4986-abf1-fca297dbccb6.root",
dir+"/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/2c0ec819-fba2-4216-a000-cd36f8a87e6a.root"
 )
)
'''

#process.load("Run3ScoutingAnalysisTools.ScoutingFilter.ScoutingFilter_cff")

#process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
#process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("mmgTree.root")
)

#process.ScoutingFilterPath = cms.Path(process.scoutingFilter)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4', '')

#L1Info = ["L1_DoubleMu4_SQ_OS_dR_Max1p2","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6","L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6(5)","L1_DoubleMu0er1p5(4)_SQ_OS_dR_Max1p4"]
L1Info = ["L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4","L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6","L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6","L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5","L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4","L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4","L1_DoubleMu4p5_SQ_OS_dR_Max1p2","L1_DoubleMu4_SQ_OS_dR_Max1p2"]
process.tree = cms.EDAnalyzer('MuMuGammaTreeMaker',
                                      triggerresults   = cms.InputTag("TriggerResults", "", "HLT"),
                                      ReadPrescalesFromFile = cms.bool( False ),
                                      AlgInputTag       = cms.InputTag("gtStage2Digis"),
                                      l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
                                      doL1 = cms.bool( True ),
                                      l1Seeds           = cms.vstring(L1Info),
                                      muons            = cms.InputTag("slimmedMuons"),
                                      electrons        = cms.InputTag("slimmedElectrons"),
                                      photons          = cms.InputTag("slimmedPhotons"),
                                      pfcands          = cms.InputTag("packedPFCandidates"),
                                      prunedGenParticles  = cms.InputTag("prunedGenParticles"),
                                      packedGenParticles  = cms.InputTag("packedGenParticles"),
                                      doGEN = cms.bool( True ),
                                      doFullGEN = cms.bool( True ),
                                      primaryVertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                      displacedVertices  = cms.InputTag("slimmedSecondaryVertices"),
                                  )

#process.p = cms.Path(process.gtStage2Digis+process.scoutingTree)
process.p = cms.Path(process.tree)
