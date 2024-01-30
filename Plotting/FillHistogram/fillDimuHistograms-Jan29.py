#!/usr/bin/python3                                                                                                                                            
#-----------------------------------------------   

import ROOT
from ROOT import *
from math import *
#import gfal2                                                                                                      
import sys, os, pwd, subprocess, glob, fnmatch
import optparse, shlex, re
import time
from time import gmtime, strftime
import math
from array import array

Z_MASS = 91.1876
ETA_MASS = 0.547862
MU_MASS = 0.105658

#define function for parsing options
def parseOptions():

	usage = ('usage: %prog [options]\n'
			+ '%prog -h for help')
	parser = optparse.OptionParser(usage)

	# input options
	parser.add_option('-i', '--input', dest='INPUT', type='string', help='input file')
	parser.add_option('-o', '--output', dest='OUTPUT', type='string', help='output file')
	parser.add_option('-n', '--njobs', dest='NJOBS', type=int, help='njobs')
	parser.add_option('-j', '--job', dest='JOB', type=int, help='job')
	parser.add_option('-l', dest='LIST', help='ntuple list', default="/afs/cern.ch/user/j/jfriesen/CMSSW_13_0_10/src/Run3DimuonAnalysisTools/Plotting/FillHistogram/muMuGammaTree_ntuples_fullRun3.txt")

	# store options and arguments as global variables
	global opt, args
	(opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quite = 0):
	status, output = subprocess.getstatusoutput(cmd)
	if (status !=0 and not quite):
		print('Error in processing command:\n 	['+cmd+']')
		print('Output:\n   ['+output+'] \n')
	return output

def fillHistogram():

	global opt, args
	parseOptions()

	ROOT.gROOT.SetBatch()
	print(opt.INPUT)

	files = [ "root://cmsxrootd.fnal.gov//"+line for line in open(str(opt.LIST))]

	N = len(files)

	first = int(float(N)/float(opt.NJOBS)*float(opt.JOB-1))
	last = int(float(N)/float(opt.NJOBS)*float(opt.JOB))

	print(first, last)

	t_scoutMuon = TChain("tree/tree")
	for i in range(len(files)):
		if (i<first or i>=last): continue
		print(files[i])
		t_scoutMuon.Add(files[i])
		print(t_scoutMuon.GetEntries())

	bin_width = 0.001
	mmg_low = round(0./bin_width)*bin_width
	mmg_high = round(1./bin_width)*bin_width
	mmg_bins = round((mmg_high - mmg_low)/bin_width)

	verbose = True

	muon_selections = (
		"isHighPtMuon",
		"isLooseMuon",
		"isMediumMuon",
		"isSoftMuon",
		"isTightMuon"
	)

	eta_cut = {
		"barrel" : ( 0, 1.48 ),
		"endcap" : ( 1.48, 3 )
	}

	pt_bins = [ 7, 9, 11, 13, 15, 19, 23, 27 ]

	config = {}
	for i_bin in range(len(pt_bins)-1) :
		config[i_bin] = {}
		for cut in eta_cut :
			config[i_bin][cut] = {}
			for selection in muon_selections :
				name = "massDimu_pt" + str(pt_bins[i_bin]) + "to" + str(pt_bins[i_bin+1]) + "_" + cut + "_" + selection
				config[i_bin][cut][selection] = ROOT.TH1F(name,name,mmg_bins,mmg_low,mmg_high)

	i_event = 0
	for ev in t_scoutMuon :
		if (ev.pt < pt_bins[0] or ev.pt > pt_bins[-1]) : continue
		if(verbose or i_event%10000==0): print("passed event:",i_event)
		i_event+=1
		for i_bin in range(len(pt_bins)-1) :
			if ev.pt < pt_bins[i_bin] or ev.pt > pt_bins[i_bin+1] : continue
			for cut in eta_cut :
				if abs(ev.eta) < eta_cut[cut][0] or abs(ev.eta) > eta_cut[cut][1] : continue
				for i_selection in range(len(muon_selections)) :
					if not (ev.muonID1[i_selection] and ev.muonID2[i_selection]) : continue
					if(verbose) : print(ev.pt, ev.eta, ev.muonID1[i_selection], ev.muonID2[i_selection], "massDimu_pt" + str(pt_bins[i]) + "to" + str(pt_bins[i+1]) + "_" + cut + "_" + muon_selections[i_selection])
					config[i_bin][cut][muon_selections[i_selection]].Fill(ev.mass)

	print("saving as "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	outfile = ROOT.TFile(str(opt.OUTPUT)+str(opt.JOB)+".root", "recreate")
	for i_bin in range(len(pt_bins)-1) :
		for cut in eta_cut :
			for selection in muon_selections :
				name = "massDimu_pt" + str(pt_bins[i_bin]) + "to" + str(pt_bins[i_bin+1]) + "_" + cut + "_" + selection
				print("saving as " + name + " with",config[i_bin][cut][selection].GetEntries(),"entries")
				outfile.WriteObject(config[i_bin][cut][selection], name)

	outfile.Close()

if __name__ == "__main__":
	fillHistogram()

