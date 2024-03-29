#!/usr/bin/env python3

import ROOT, math, os
from array import array
import utils.CMSGraphics as CMSGraphics

ROOT.gStyle.SetLabelOffset(0.02, "X")
ROOT.gStyle.SetTitleOffset(2.5, "X")
ROOT.gStyle.SetTitleOffset(1.8, "Y")
ROOT.gStyle.SetTitleSize(0.025, "X")
ROOT.gStyle.SetTitleSize(0.04, "Y")
ROOT.gROOT.ForceStyle()

ETA_MASS = 0.547862

itree = ROOT.TChain()
for i in range(9)[1:]:
	itree.Add("data/MC Efficiency Studies/MCTree_" + str(int(i)) + ".root?#tree_" + str(int(i)))

config_var = {
	"mass_mumu" : {
		"nbins" : 120,
		"range" : (-0.02, 0.02),
		"max" : 0.04,
		"label" : "#frac{m_{#mu#mu}^{GEN}-m_{#mu#mu}^{RECO}}{m_{#mu#mu}^{GEN}}",
		"bins" : (	
			{ "max" : 0.4, "label" : "m_{#mu#mu}^{GEN}<0.4 GeV" },
			{ "max" : 0.6, "label" : "0.4<m_{#mu#mu}^{GEN}<0.6 GeV" },
			{ "max" : 0.8, "label" : "0.6<m_{#mu#mu}^{GEN}<0.8 GeV" },
			{ "max" : 1.1, "label" : "0.8<m_{#mu#mu}^{GEN}<1.1 GeV" },
			{ "max" : 1.4, "label" : "1.1<m_{#mu#mu}^{GEN}<1.4 GeV" },
		)
	},
	"pt_mumu" : {
		"nbins" : 60,
		"range" : (-0.05, 0.05),
		"max" : 0.08,
		"label" : "Dimuon #frac{p_{T}^{GEN}-p_{T}^{RECO}}{p_{T}^{GEN}}",
		"bins" : (	
			{ "max" : 8, "label" : "p_{T}^{GEN}<8 GeV" },
			{ "max" : 14, "label" : "8<p_{T}^{GEN}<14 GeV" },
			{ "max" : 22, "label" : "14<p_{T}^{GEN}<22 GeV" },
			{ "max" : 32, "label" : "22<p_{T}^{GEN}<32 GeV" },
			{ "max" : 100, "label" : "32 GeV<p_{T}^{GEN}" },
		)
	},
	"pt_leading" : {
		"nbins" : 60,
		"range" : (-0.05, 0.05),
		"max" : 0.08,
		"label" : "Leading muon #frac{p_{T}^{GEN}-p_{T}^{RECO}}{p_{T}^{GEN}}",
		"bins" : (	
			{ "max" : 4, "label" : "p_{T}^{GEN}<4 GeV" },
			{ "max" : 6, "label" : "4<p_{T}^{GEN}<6 GeV" },
			{ "max" : 8, "label" : "6<p_{T}^{GEN}<8 GeV" },
			{ "max" : 10, "label" : "8<p_{T}^{GEN}<10 GeV" },
			{ "max" : 100, "label" : "10 GeV<p_{T}^{GEN}" },
		)
	},
	"pt_trailing" : {
		"nbins" : 60,
		"range" : (-0.05, 0.05),
		"max" : 0.08,
		"label" : "Trailing muon #frac{p_{T}^{GEN}-p_{T}^{RECO}}{p_{T}^{GEN}}",
		"bins" : (	
			{ "max" : 4, "label" : "p_{T}^{GEN}<4 GeV" },
			{ "max" : 6, "label" : "4<p_{T}^{GEN}<6 GeV" },
			{ "max" : 8, "label" : "6<p_{T}^{GEN}<8 GeV" },
			{ "max" : 10, "label" : "8<p_{T}^{GEN}<10 GeV" },
			{ "max" : 100, "label" : "10 GeV<p_{T}^{GEN}" },
		)
	},
	"dr" : {
		"nbins" : 80,
		"range" : (-0.05, 0.05),
		"max" : 0.12,
		"label" : "Dimuon #frac{#DeltaR^{GEN}-#DeltaR^{RECO}}{#DeltaR^{GEN}}",
		"bins" : (	
			{ "max" : 0.05, "label" : "#DeltaR^{GEN}<0.05" },
			{ "max" : 0.1, "label" : "0.05<#DeltaR^{GEN}<0.1" },
			{ "max" : 0.15, "label" : "0.1<#DeltaR^{GEN}<0.15" },
			{ "max" : 0.2, "label" : "0.15<#DeltaR^{GEN}<0.2" },
			{ "max" : 100, "label" : "0.2<#DeltaR^{GEN}" },
		)
	},
}

colors = (
	ROOT.kOrange,
	ROOT.kAzure,
	ROOT.kRed,
	ROOT.kViolet,
	ROOT.kGray,
	ROOT.kPink,
	ROOT.kTeal,
	ROOT.kMagenta,
	ROOT.kSpring,
	ROOT.kBlue,
	ROOT.kYellow,
)

p = "delta_over_truth"
for v in config_var :
	for b in range(len(config_var[v]["bins"])) :
		name = p+"_"+v+"_"+str(b)
		if b == 0 :
			bin_range = ("0", str(config_var[v]["bins"][b]["max"]))
			condition = bin_range[0] + " < " + "gen_matched" + "_" + v + " && " + bin_range[1] + " > " + "gen_matched" + "_" + v
		else :
			bin_range = (str(config_var[v]["bins"][b-1]["max"]), str(config_var[v]["bins"][b]["max"]))
		if str(config_var[v]["bins"][b]["max"]) == 100 :
			condition = bin_range[0] + " < " + "gen_matched" + "_" + v
		else :
			condition = bin_range[0] + " < " + "gen_matched" + "_" + v + " && " + bin_range[1] + " > " + "gen_matched" + "_" + v
		config_var[v]["bins"][b]["histo"] = ROOT.TH1F( name, name, config_var[v]["nbins"], config_var[v]["range"][0], config_var[v]["range"][1] )
		itree.Draw( p + "_" + v + ">>" + name, condition )
		config_var[v]["bins"][b]["histo"].Sumw2()
		if config_var[v]["bins"][b]["histo"].GetEntries() > 0 : config_var[v]["bins"][b]["histo"].Scale(1/config_var[v]["bins"][b]["histo"].GetEntries())

for v in config_var :
	canvas = CMSGraphics.makeCMSCanvas(v,v,1400,1200)
	canvas.SetLeftMargin(0.15)
	canvas.SetBottomMargin(0.15)
	canvas.cd()
	legend = CMSGraphics.makeLegend(nentries=len(config_var[v]["bins"]), left=0.68, margin=0.2)
	it = range(len(config_var[v]["bins"])) #if v != "mass_mumu" else reversed(range(len(config_var[v]["bins"])))
	for b in it :
		print(config_var[v]["bins"][b]["histo"].GetEntries())
		if "max" in config_var[v] :
			config_var[v]["bins"][b]["histo"].SetMaximum(config_var[v]["max"])
		config_var[v]["bins"][b]["histo"].SetMinimum(0)
		config_var[v]["bins"][b]["histo"].SetLineWidth(3)
		config_var[v]["bins"][b]["histo"].SetLineColor(colors[b])
		if b == 0 :
			config_var[v]["bins"][b]["histo"].SetFillColorAlpha(colors[b],0.25)
			config_var[v]["bins"][b]["histo"].SetFillStyle(3011)
		if b == 4 :
			config_var[v]["bins"][b]["histo"].SetFillColorAlpha(colors[b],0.25)
		config_var[v]["bins"][b]["histo"].GetXaxis().SetTitle(config_var[v]["label"])
		config_var[v]["bins"][b]["histo"].GetYaxis().SetTitle("Normalized Events")
		config_var[v]["bins"][b]["histo"].Draw("H same")
		legend.AddEntry(config_var[v]["bins"][b]["histo"], config_var[v]["bins"][b]["label"], "f")
	CMSGraphics.printLumiLeft(canvas, extraText="Simulation Preliminary")
	legend.Draw("same")
	canvas.Update()
	canvas.SaveAs("plots/resolution_"+v+".png")

