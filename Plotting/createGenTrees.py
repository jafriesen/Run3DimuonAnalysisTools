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

itree = TChain("tree/tree")
itree.Add("")
ofile = ROOT.TFile("MCTree.root", "RECREATE")
otree = ROOT.TTree("tree", "tree")

l_decay = (
	"other",
	"isEta2MuMu",
	"isEta2MuMuGamma",
	"isEtaPrime2MuMu",
	"isEtaPrime2MuMuGamma",
	"isOmega2MuMu",
	"isOmega2Pi0MuMu",
	"isRho2MuMu",
	"isPhi2MuMu",
	"isPhi2KK",
)
l_var = (
	"gen_mumu_mass",
	"mass_mumu_match1",
	"mass_mumu_match2",
	"delta_mass_mumu"
	"delta_mass_mumu_over_truth",
)

arrs = {}
for d in l_decay :
	arrs[d] = {}
	for v in l_var :
		name = d + "_" + v
		arrs[d][v] = array("d",[0])
		t.Branch(name, arrs[d][v], name+"/D")

for ev in itree :
	for i_dimu in range(len(ev.gen_motherID)) :

		for d in l_decay :
			for v in l_var :
				arrs[d][v][0] = -1
		decay = "other"
		if ev.gen_isEta2MuMu[i_dimu] : decay = "isEta2MuMu"
		if ev.gen_isEta2MuMuGamma[i_dimu] : decay = "isEta2MuMuGamma"
		if ev.gen_isEtaPrime2MuMu[i_dimu] : decay = "isEtaPrime2MuMu"
		if ev.gen_isEtaPrime2MuMuGamma[i_dimu] : decay = "isEtaPrime2MuMuGamma"
		if ev.gen_isOmega2MuMu[i_dimu] : decay = "isOmega2MuMu"
		if ev.gen_isOmega2Pi0MuMu[i_dimu] : decay = "isOmega2Pi0MuMu"
		if ev.gen_isRho2MuMu[i_dimu] : decay = "isRho2MuMu"
		if ev.gen_isPhi2MuMu[i_dimu] : decay = "isPhi2MuMu"
		if ev.gen_isPhi2KK[i_dimu] : decay = "isPhi2KK"

		arrs[decay]["gen_mumu_mass"] = ev.gen_mumu_mass[i_dimu]
		if ev.gen_mu1_recoMatch[i_dimu] or ev.gen_mu2_recoMatch[i_dimu] : arrs[decay]["mass_mumu_match1"] = ev.gen_mumu_mass[i_dimu]
		if ev.gen_mu1_recoMatch[i_dimu] and ev.gen_mu2_recoMatch[i_dimu] :
			arrs[decay]["mass_mumu_match2"] = ev.gen_mumu_mass[i_dimu]
			arrs[decay]["delta_mass_mumu"] = ev.gen_mumu_mass[i_dimu] - ev.mumu_mass[i_dimu]
			arrs[decay]["delta_mass_mumu_over_truth"] = ( ev.gen_mumu_mass[i_dimu] - ev.mumu_mass[i_dimu] ) / ev.gen_mumu_mass[i_dimu]

		t.Fill()

ofile.Write()
ofile.Close()

