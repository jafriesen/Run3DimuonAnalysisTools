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

	verbose = False
	vtxCut = False

	muon_selections = [
		[0, "isHighPtMuon"],
		[1, "isLooseMuon"],
		[2, "isMediumMuon"],
		[3, "isSoftMuon"],
		[4, "isTightMuon"]
	]
	vtxProb_cuts = [
		[0.005, "0.005"],
		[0.01, "0p01"],
		[0.03, "0p03"],
		[0.06, "0p06"],
		[0.1, "0p1"],
		[0.15, "0p15"],
		[0.2, "0p2"],
		[0.25, "0p25"],
		[0.3, "0p3"]
	]

	config = {}
	for cut in vtxProb_cuts :
		config[cut[1]] = {}
		name = "massDimu_anyMuon_minProbVtx" + cut[1]
		config[cut[1]]["any"] = ROOT.TH1F(name,name,mmg_bins,mmg_low,mmg_high)
		for selection in muon_selections :
			name = "massDimu_" + selection[1] + "_minProbVtx" + cut[1]
			config[cut[1]][selection[1]] = ROOT.TH1F(name,name,mmg_bins,mmg_low,mmg_high)

	i_event = 0
	for ev in t_scoutMuon :
		if (ev.mass > 1.0) : continue
		if(verbose or i_event%10000==0): print("m2 < 1 event:",i_event)
		i_event+=1
		for cut in vtxProb_cuts :
			if ev.probVtx < cut[0] : break
			if(verbose) : print("	probVtx",ev.probVtx,"cut",cut[0])
			config[cut[1]]["any"].Fill(ev.mass)
			for selection in muon_selections :
				if(verbose) : print("		selection",selection[0],ev.muonID1[selection[0]],ev.muonID2[selection[0]])
				if(ev.muonID1[selection[0]] and ev.muonID2[selection[0]]) : 
					if(verbose) : print("			fill")
					config[cut[1]][selection[1]].Fill(ev.mass)

	print("saving as "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	outfile = ROOT.TFile(str(opt.OUTPUT)+str(opt.JOB)+".root", "recreate")
	for cut in vtxProb_cuts :
		name = "massDimu_anyMuon_minProbVtx" + cut[1]
		print("saving as " + name + " with",config[cut[1]]["any"].GetEntries(),"entries")
		outfile.WriteObject(config[cut[1]]["any"], name)
		for selection in muon_selections :
			name = "massDimu_" + selection[1] + "_minProbVtx" + cut[1]
			print("saving as " + name + " with",config[cut[1]][selection[1]].GetEntries(),"entries")
			outfile.WriteObject(config[cut[1]][selection[1]], name)

	outfile.Close()

if __name__ == "__main__":
	fillHistogram()

