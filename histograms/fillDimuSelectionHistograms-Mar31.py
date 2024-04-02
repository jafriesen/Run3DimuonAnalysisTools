#!/usr/bin/python3
#-----------------------------------------------

import subprocess, optparse
from pathlib import Path
import ROOT, math
from array import array

# potentially useful constants
ETA_MASS = 0.547862
MU_MASS = 0.105658

# define function for parsing options
def parseOptions() :
	usage = ('usage: %prog [options]\n'
			+ '%prog -h for help')
	parser = optparse.OptionParser(usage)
	# input options
	parser.add_option('-i', '--input', dest='INPUT', type='string', help='input file', default="")
	parser.add_option('-o', '--output', dest='OUTPUT', type='string', help='output file')
	parser.add_option('-n', '--njobs', dest='NJOBS', type=int, help='total njobs')
	parser.add_option('-j', '--job', dest='JOB', type=int, help='job index')
	parser.add_option('-l', dest='LIST', help='input file list', default="/afs/cern.ch/user/j/jfriesen/CMSSW_13_0_10/src/Run3DimuonAnalysisTools/Plotting/FillHistogram/muMuGammaTree_ntuples_fullRun3.txt")
	# store options and arguments as global variables
	global opt, args
	(opt, args) = parser.parse_args()

# define function for processing the external os commands
def processCmd(cmd, quite = 0) :
	status, output = subprocess.getstatusoutput(cmd)
	if (status !=0 and not quite):
		print('Error in processing command:\n 	['+cmd+']')
		print('Output:\n   ['+output+'] \n')
	return output

# define function for job
def fillHistograms() :

	verbose = False
	ROOT.gROOT.SetBatch()

	global opt, args
	parseOptions()

	# get input files from list (parsed option)
	print( "Opening list of input files", opt.LIST )
	files = Path(opt.LIST).read_text().splitlines()
	print( "	Found", len(files), "items in list")



	njobs_actual = min( opt.NJOBS, len(files) )
	file_range = (
		int( ( (float(opt.JOB)-1) * len(files) ) / njobs_actual ),
		int( ( (float(opt.JOB)) * len(files) ) / njobs_actual )
	)

	N = len(files)
	file_range = ( int(float(N)/float(opt.NJOBS)*float(opt.JOB-1)), int(float(N)/float(opt.NJOBS)*float(opt.JOB)) )
	
	print( "Job", opt.JOB, "of", opt.NJOBS, "using file range", file_range )

	redirector = "root://cmsxrootd.fnal.gov//"
	print("Using redirector", redirector)

	itree_name = "tree/tree"
	itree = TChain(itree_name)
	for i in range(len(files)):
		if (i<first or i>=last): continue
		print("Getting", itree_name, "from", redirector+files[i])
		itree.Add(redirector+files[i])
		print(itree.GetEntries(), "total entries in TChain")

	print("Creating "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	outfile = ROOT.TFile(str(opt.OUTPUT)+str(opt.JOB)+".root", "recreate")

	# define histogram binning
	bin_width = 0.0001 # GeV
	mumu_mass_range = ( round(0./bin_width) * bin_width, round(1.4/bin_width) * bin_width )
	mumu_mass_bins = round( (mumu_mass_range[1] - mumu_mass_range[0])/bin_width )



	muon_IDs = {
		"anyID" : {
			"ID" : -1,
		},
		"isHighPtMuon" : {
			"ID" : 0,
		},
		"isLooseMuon" : {
			"ID" : 1,
		},
		"isMediumMuon" : {
			"ID" : 2,
		},
		"isSoftMuon" : {
			"ID" : 3,
		},
		"isTightMuon" : {
			"ID" : 4,
		},
	}

	vtxProb_bins = [ 0.005, 0.0075, 0.01, 0.015, 0.02, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]

	for vtxProb in vtxProb_bins :
		for ID in muon_IDs :
			if vtxProb == vtxProb_bins[-1] :
				name = "mumu_mass_" + ID




	config = {}
	for cut in vtxProb_cuts :


		config[cut[1]] = {}
		name = "massDimu_anyMuon_minProbVtx" + cut[1]
		config[cut[1]]["any"] = ROOT.TH1F(name,name,mumu_mass_bins,mumu_mass_low,mumu_mass_high)
		for selection in muon_selections :
			name = "massDimu_" + selection[1] + "_minProbVtx" + cut[1]
			config[cut[1]][selection[1]] = ROOT.TH1F(name,name,mumu_mass_bins,mumu_mass_low,mumu_mass_high)

	print("Processing events...")
	i_event = 0
	for ev in itree :
		if (ev.mass > mumu_mass_high) : continue
		if(verbose or i_event%10000==0): print("mumu_mass < 1 event:",i_event)
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

	for cut in vtxProb_cuts :
		name = "massDimu_anyMuon_minProbVtx" + cut[1]
		print("Saving histogram " + name + " with",config[cut[1]]["any"].GetEntries(),"entries")
		outfile.WriteObject(config[cut[1]]["any"], name)
		for selection in muon_selections :
			name = "massDimu_" + selection[1] + "_minProbVtx" + cut[1]
			print("Saving histogram " + name + " with",config[cut[1]][selection[1]].GetEntries(),"entries")
			outfile.WriteObject(config[cut[1]][selection[1]], name)

	outfile.Close()

if __name__ == "__main__":
	fillHistograms()

