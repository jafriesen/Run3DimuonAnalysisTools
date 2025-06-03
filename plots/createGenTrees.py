#!/usr/bin/python3                                                                                                                                            
#-----------------------------------------------   

import ROOT
from math import *
import math
from array import array

Z_MASS = 91.1876
ETA_MASS = 0.547862
MU_MASS = 0.105658

ROOT.gROOT.SetBatch()

itree = ROOT.TChain("tree/tree")
mc_path = "/eos/user/j/jfriesen/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/crab_muMuGammaTree_Run3Summer22MiniAODv3_InclusiveDileptonMinBias_MINIAODSIM_19Dec2023/240110_204550"
itree.Add(f"{mc_path}/mmgTree0000.root")
itree.Add(f"{mc_path}/mmgTree0001.root")

l_level = (
	"gen",
	"gen_matched",
	"gen_partiallymatched",
	"reco_matched",
	"reco_partiallymatched",
	"delta",
	"delta_over_truth"
)
l_var = (
	"mass_mumu",
	"pt_mumu",
	"pt_leading",
	"pt_trailing",
	"dr",
	"eta1",
	"eta2"
)


save = False


i_ev = 0
for ev in itree :

	if i_ev % 5000000 == 0 :
		if save : 
			if i_ev > 0 :
				print("saving")
				ofile.Write()
				ofile.Close()
			ofile = ROOT.TFile("MCTree_" + str(int(i_ev / 5000000)) + ".root", "RECREATE")
		otree = ROOT.TTree("tree_" + str(int(i_ev / 5000000)), "tree" + str(int(i_ev / 5000000)))
		arrs = {}
		arrs["decay"] = array("i",[0])
		otree.Branch("decay", arrs["decay"], "decay/I")
		for l in l_level :
			arrs[l] = {}
			for v in l_var :
				name = l + "_" + v
				arrs[l][v] = array("d",[0])
				otree.Branch(name, arrs[l][v], name+"/D")
	if i_ev % 100000 == 0 : print(i_ev)
	i_ev += 1

	#print(i_ev, ev.motherID1, ev.motherID2, ev.mumu_mass, ev.gen_motherID, ev.gen_mumu_mass)
	for i_dimu in range(len(ev.gen_motherID)) :

		nPhotons = len(ev.gen_matchedPhotonPt[i_dimu])
		'''
		decay = "other"
		if ev.gen_isEta2MuMu[i_dimu] : decay = "isEta2MuMu"
		if ev.gen_isEta2MuMuGamma[i_dimu] : decay = "isEta2MuMuGamma"
		if ev.gen_isEtaPrime2MuMu[i_dimu] : decay = "isEtaPrime2MuMu"
		if ev.gen_isEtaPrime2MuMuGamma[i_dimu] : decay = "isEtaPrime2MuMuGamma"
		if ev.gen_isOmega2MuMu[i_dimu] : decay = "isOmega2MuMu"
		if abs(ev.gen_motherID[i_dimu]) == 223 and (nPhotons == 1) : decay = "isOmega2MuMuGamma"
		if ev.gen_isOmega2Pi0MuMu[i_dimu] : decay = "isOmega2Pi0MuMu"
		if ev.gen_isRho2MuMu[i_dimu] : decay = "isRho2MuMu"
		if abs(ev.gen_motherID[i_dimu]) == 113 and (nPhotons == 1) : decay = "isRho2MuMuGamma"
		if ev.gen_isPhi2MuMu[i_dimu] : decay = "isPhi2MuMu"
		if abs(ev.gen_motherID[i_dimu]) == 333 and (nPhotons == 1) : decay = "isPhi2MuMuGamma"
		if ev.gen_isPhi2KK[i_dimu] : decay = "isPhi2KK"
		'''
		decay_flags = (
			1 * ev.gen_isEta2MuMu[i_dimu],
			2 * ev.gen_isEta2MuMuGamma[i_dimu],
			3 * ev.gen_isEtaPrime2MuMu[i_dimu],
			4 * ev.gen_isEtaPrime2MuMuGamma[i_dimu],
			5 * ev.gen_isOmega2MuMu[i_dimu],
			6 * ( abs(ev.gen_motherID[i_dimu]) == 223 and (nPhotons == 1) ),
			7 * ev.gen_isOmega2Pi0MuMu[i_dimu],
			8 * ev.gen_isRho2MuMu[i_dimu],
			9 * ( abs(ev.gen_motherID[i_dimu]) == 113 and (nPhotons == 1) ),
			10 * ev.gen_isPhi2MuMu[i_dimu],
			11 * ( abs(ev.gen_motherID[i_dimu]) == 333 and (nPhotons == 1) ),
			12 * ev.gen_isPhi2KK[i_dimu]
		)
		decay = 0
		for flag in decay_flags :
			if flag > 0 : decay = flag
		arrs["decay"][0] = decay
		#print(decay)

		v_mu1_gen = ROOT.Math.PtEtaPhiMVector(ev.gen_mu1_pt[i_dimu], ev.gen_mu1_eta[i_dimu], ev.gen_mu1_phi[i_dimu], MU_MASS)
		v_mu2_gen = ROOT.Math.PtEtaPhiMVector(ev.gen_mu2_pt[i_dimu], ev.gen_mu2_eta[i_dimu], ev.gen_mu2_phi[i_dimu], MU_MASS)
		v_dimu_gen = v_mu1_gen+v_mu2_gen
		dr_gen = abs(ROOT.Math.VectorUtil.DeltaR(v_mu1_gen, v_mu2_gen))
		variables = {
			"mass_mumu" : {
				"gen" :		ev.gen_mumu_mass[i_dimu],
				"reco" : 	ev.mumu_mass
			},
			"pt_mumu" : {
				"gen" :		v_dimu_gen.Pt(),
				"reco" :	ev.mumu_pt
			},
			"dr" : {
				"gen" : 	dr_gen,
				"reco" :	ev.mumu_deltaR
			}
		}
		pt_gen = ( ev.gen_mu1_pt[i_dimu], ev.gen_mu2_pt[i_dimu] )
		pt_reco = ( ev.mu1_pt, ev.mu2_pt )
		i_leading = 0 if ev.gen_mu1_pt[i_dimu] > ev.gen_mu2_pt[i_dimu] else 1
		variables["pt_leading"] = {
			"gen" : 	pt_gen[i_leading],
			"reco" : 	pt_reco[i_leading],
		}
		variables["pt_trailing"] = {
			"gen" : 	pt_gen[1-i_leading],
			"reco" : 	pt_reco[1-i_leading],
		}
		eta_gen = ( ev.gen_mu1_eta[i_dimu], ev.gen_mu2_eta[i_dimu] )
		eta_reco = ( ev.mu1_eta, ev.mu2_eta )
		variables["eta1"] = {
			"gen" : 	eta_gen[i_leading],
			"reco" : 	eta_reco[i_leading],
		}
		variables["eta2"] = {
			"gen" : 	eta_gen[1-i_leading],
			"reco" : 	eta_reco[1-i_leading],
		}

		if ev.mumu_mass > 0 and ev.mumu_mass < 0.6 :
			print()
			v_mu1 = ROOT.Math.PtEtaPhiMVector(ev.mu1_pt, ev.mu1_eta, ev.mu1_phi, MU_MASS)
			v_mu2 = ROOT.Math.PtEtaPhiMVector(ev.mu2_pt, ev.mu2_eta, ev.mu2_phi, MU_MASS)
			dr1 = ROOT.Math.VectorUtil.DeltaR(v_mu1_gen, v_mu1)
			dr2 = ROOT.Math.VectorUtil.DeltaR(v_mu2_gen, v_mu2)
			print("decay", decay, "mass", ev.mumu_mass)
			print("match	", ev.gen_mu2_recoMatch[i_dimu], ev.gen_mu1_recoMatch[i_dimu])
			print("dr	", dr1, dr2)
			print("reco 1	", ev.mu1_pt, ev.mu1_eta, ev.mu1_phi)
			print("gen 1	", ev.gen_mu1_pt[i_dimu], ev.gen_mu1_eta[i_dimu], ev.gen_mu1_phi[i_dimu])
			print("reco 2	", ev.mu2_pt, ev.mu2_eta, ev.mu2_phi)
			print("gen 2	", ev.gen_mu2_pt[i_dimu], ev.gen_mu2_eta[i_dimu], ev.gen_mu2_phi[i_dimu])


		for v in l_var :
			if ev.gen_mu1_recoMatch[i_dimu] or ev.gen_mu2_recoMatch[i_dimu] :
				arrs["gen_partiallymatched"][v][0] = variables[v]["gen"]
				arrs["reco_partiallymatched"][v][0] = variables[v]["reco"]
			if ev.gen_mu1_recoMatch[i_dimu] and ev.gen_mu2_recoMatch[i_dimu] :
				arrs["gen_matched"][v][0] = variables[v]["gen"]
				arrs["reco_matched"][v][0] = variables[v]["reco"]
				arrs["delta"][v][0] = variables[v]["gen"] - variables[v]["reco"]
				arrs["delta_over_truth"][v][0] = ( variables[v]["gen"] - variables[v]["reco"] ) / variables[v]["gen"]
			else :
				for l in l_level :
					arrs[l][v][0] = -9999
			arrs["gen"][v][0] = variables[v]["gen"]

		otree.Fill()


if save :
	print("saving")
	ofile.Write()
	ofile.Close()