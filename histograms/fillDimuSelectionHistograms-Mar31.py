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
	file_list = Path(opt.LIST).read_text().splitlines()
	print( "	Found", len(file_list), "items in list")

	njobs_actual = min( opt.NJOBS, len(file_list) )
	file_range = (
		int( ( (float(opt.JOB)-1) * len(file_list) ) / njobs_actual ),
		int( ( (float(opt.JOB)) * len(file_list) ) / njobs_actual )
	)
	
	print( "	Job", opt.JOB, "of", opt.NJOBS, "using file range", file_range )

	redirector = "root://cmsxrootd.fnal.gov//"
	print( "	Using redirector", redirector )

	itree_name = "tree/tree"
	itree = ROOT.TChain(itree_name)
	for i in range(file_range[0], file_range[1]) :
		print("Getting", itree_name, "from", redirector+file_list[i])
		itree.Add(redirector+file_list[i])
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

	for i_vtxProb in range(len(vtxProb_bins)) :
		for ID in muon_IDs :
			if vtxProb_bins[i_vtxProb] == vtxProb_bins[-1] :
				name = "mumu_mass_" + ID + "_vtxProbMin" + str(vtxProb_bins[i_vtxProb]).replace(".","p")
			else :
				name = "mumu_mass_" + ID + "_vtxProb" + str(vtxProb_bins[i_vtxProb]).replace(".","p") + "to" + str(vtxProb_bins[i_vtxProb+1]).replace(".","p")
			print("Creating histogram", name)
			muon_IDs[ID][i_vtxProb] = ROOT.TH1F(name,name,mumu_mass_bins,mumu_mass_low,mumu_mass_high)

	print("Processing events...")
	i_event = 0
	for ev in itree :
		if (ev.mass > mumu_mass_high) : continue
		if(verbose or i_event%10000==0): print("mumu_mass < mumu_mass_high event:",i_event)
		i_event+=1
		for i_vtxProb in range(len(vtxProb_bins)) :
			if vtxProb_bins[i_vtxProb] != vtxProb_bins[-1] and ev.probVtx > vtxProb_bins[i_vtxProb + 1] : continue
			if ev.probVtx < vtxProb_bins[i_vtxProb] : break
			for ID in muon_IDs : 
				if muon_IDs[ID]["ID"] == -1 or ( ev.muonID1[muon_IDs[ID]["ID"]] and ev.muonID2[muon_IDs[ID]["ID"]] ) :
					muon_IDs[ID][i_vtxProb].Fill(ev.mass)

	for i_vtxProb in range(len(vtxProb_bins)) :
		for ID in muon_IDs :
			name = muon_IDs[ID][i_vtxProb].GetName()
			print( "Saving", name, "to", str(opt.OUTPUT)+str(opt.JOB)+".root", "with", muon_IDs[ID][i_vtxProb].GetEntries(), "entries" )
			outfile.WriteObject( muon_IDs[ID][i_vtxProb], name )

	outfile.Close()

if __name__ == "__main__":
	fillHistograms()

