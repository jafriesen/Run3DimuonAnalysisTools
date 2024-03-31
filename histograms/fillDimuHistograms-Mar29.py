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
from pathlib import Path

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

	print("Getting list of input files from ", opt.LIST)
	files = Path(opt.LIST).read_text().splitlines()

	N = len(files)
	first = int(float(N)/float(opt.NJOBS)*float(opt.JOB-1))
	last = int(float(N)/float(opt.NJOBS)*float(opt.JOB))
	print("Total number of input files:", N)
	print("Job", opt.JOB, "of", opt.NJOBS, "using files", first, "through", last-1)

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

	variables = {
		"event" : { "arr" : array("i",[0]) },
		"luminosityBlock" : { "arr" : array("i",[0]) },
		"run" : { "arr" : array("i",[0]) },
	}
	otree = ROOT.TTree("t_event_info","t_event_info")
	for v in variables : otree.Branch(v, variables[v]["arr"], v+"/I")

	verbose = False

	print("Processing events...")
	i_event = 0
	for ev in itree :
		if verbose : print( "event", ev.eventNum, "luminosityBlock", ev.lumiSec, "run", ev.runNum )
		variables["event"]["arr"][0] = ev.eventNum
		variables["luminosityBlock"]["arr"][0] = ev.lumiSec
		variables["run"]["arr"][0] = ev.runNum
		otree.Fill()

	print("Saving output to "+str(opt.OUTPUT)+str(opt.JOB)+".root")
	outfile.Write()
	outfile.Close()

if __name__ == "__main__":
	fillHistogram()

