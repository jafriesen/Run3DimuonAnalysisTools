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
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

ofile = ROOT.TFile( "data/histosMar29_v2.root", "RECREATE" )

itree = ROOT.TChain()
itree.Add("data/histosMar29.root?#t_event_info")

run_range = ( 355794, 372415 )
h_runs = ROOT.TH1F( "h_runs", "h_runs", run_range[1]-run_range[0], run_range[0], run_range[1] )

lumi_block_range = ( 0, 5000 )
h_lumi_blocks = ROOT.TH1F( "h_lumi_blocks", "h_lumi_blocks", lumi_block_range[1]-lumi_block_range[0], lumi_block_range[0], lumi_block_range[1] )

nentries = itree.GetEntries()

i_ev = 0
for ev in itree :
	if not i_ev % 100000 : print("event", i_ev, "of", nentries)
	h_runs.Fill(ev.run)
	h_lumi_blocks.Fill(ev.luminosityBlock)
	i_ev += 1

ofile.WriteObject(h_runs, "h_runs")
ofile.WriteObject(h_lumi_blocks, "h_lumi_blocks")
ofile.Close()
