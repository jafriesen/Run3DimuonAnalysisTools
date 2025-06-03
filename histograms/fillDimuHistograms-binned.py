#!/usr/bin/python3
#-----------------------------------------------

import subprocess, optparse
from pathlib import Path
import ROOT, math
from array import array
import bisect

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

pt_bins_bb = (7.36, 8.36, 9.00, 9.65, 10.32, 11.10, 12.04, 13.33, 15.28, 18.75)
pt_bins_ee = (7.06, 8.30, 9.48, 10.98, 13.58)


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
	
	print( f"	Job {opt.JOB} of {njobs_actual} using files ({file_range[0]}, {file_range[1]}]")

	redirector = "root://cmsxrootd.fnal.gov//"
	redirector = "root://xrootd.cmsaf.mit.edu:1094//"
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
	mumu_mass_range = ( round(0./bin_width) * bin_width, round(0.6/bin_width) * bin_width )
	mumu_mass_bins = round( (mumu_mass_range[1] - mumu_mass_range[0])/bin_width )

	histos = {}

	print("Processing events...")
	i_event = 0
	#fill = 0
	for ev in itree :
		if (ev.mass > mumu_mass_range[1]) : continue
		#if ( ev.softMvaRun3Value1 < 0.3 or ev.softMvaRun3Value2 < 0.3) : continue
		if(verbose or i_event%1000==0): print( "mumu_mass <", mumu_mass_range[1], "event", i_event )
		i_event+=1

		#print(f"mass: {ev.mass}, pt: {ev.pt}, eta1: {ev.eta1}, eta2: {ev.eta2}, softMvaRun3Value1: {ev.softMvaRun3Value1}, {ev.custom_softMvaRun3Value1}, softMvaRun3Value2: {ev.softMvaRun3Value2}, {ev.custom_softMvaRun3Value2}")

		if abs(ev.eta1) < 1.48 and abs(ev.eta2) < 1.48 :
			pt_key = bisect.bisect(pt_bins_bb, ev.pt)
			histo_name = f"Cat_{pt_key}_bb"
		elif abs(ev.eta1) > 1.48 and abs(ev.eta2) > 1.48 and abs(ev.eta1) < 3 and abs(ev.eta2) < 3 :
			pt_key = bisect.bisect(pt_bins_ee, ev.pt)
			histo_name = f"Cat_{pt_key}_ee"
		else :
			continue

		if ev.custom_softMvaRun3Value1 > 0.6 and ev.custom_softMvaRun3Value2 > 0.6 :
			histo_name = histo_name + "_SR"
		else :
			histo_name = histo_name + "_CR"

		#print("Filling histogram", histo_name, "with pt", ev.pt, "and mass", ev.mass)

		if not histo_name in histos :
			histos[histo_name] = ROOT.TH1F(histo_name,histo_name,mumu_mass_bins,mumu_mass_range[0],mumu_mass_range[1])
		histos[histo_name].Fill(ev.mass)

	print("Saving histograms...")
	for pt_key in histos :
		histo = histos[pt_key]
		name = histo.GetName()
		if histo.GetEntries() > 0 :
			print( "	Saving", name, "to", str(opt.OUTPUT)+str(opt.JOB)+".root", "with", histo.GetEntries(), "entries" )
			outfile.WriteObject( histo, name )
		else : print( "	NOT saving", name )

	outfile.Close()

if __name__ == "__main__":
	fillHistograms()