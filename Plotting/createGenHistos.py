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

itree = ROOT.TChain()
itree.Add("MCTree_1.root?#tree_1")
itree.Add("MCTree_2.root?#tree_2")
#itree.Add("MCTree_2.root")
#itree.Add("MCTree_3.0.root")

histos_to_plot = {
	"gen" : {
		"ylabel" : "Events",
	},
	"gen_matched" : {
		"ylabel" : "Events",
	},
	"efficiency" : {
		"ylabel" : "Efficiency",
	},
	"reco_matched" : {
		"ylabel" : "Events",
	},
	"delta_over_truth" : {
		"ylabel" : "Events",
	},
}

vars_to_plot = {
	"mass_mumu" : {
		"nbins" : 75,
		"range" : (0.2, 1.5),
		"nbins2" : 50,
		"xlabel" : "m(#mu#mu) [GeV]",
		"xlabel_delta_over_truth" : "#Deltam(#mu#mu) / m_{truth}(#mu#mu)",
		"bins_delta_over_truth" : [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 100],
		"min" : 2,
	},
	"pt_mumu" : {
		"nbins" : 50,
		"range" : (0.0, 40),
		"nbins2" : 50,
		"xlabel" : "p_{T}(#mu#mu) [GeV]",
		"xlabel_delta_over_truth" : "#Deltap_{T}(#mu#mu) / p_{truth}(#mu#mu)",
		"bins_delta_over_truth" : [2, 4, 8, 14, 22, 32, 100],
		"min" : 1,
	},
	"pt_leading" : {
		"nbins" : 50,
		"range" : (2.0, 20),
		"nbins2" : 50,
		"xlabel" : "Leading p_{T}(#mu) [GeV]",
		"xlabel_delta_over_truth" : "Leading #Deltap_{T}(#mu) / p_{truth}(#mu)",
		"bins_delta_over_truth" : [2, 4, 6, 8, 10, 10, 100],
		"min" : 1,
	},
	"pt_trailing" : {
		"nbins" : 50,
		"range" : (2.0, 6),
		"nbins2" : 50,
		"xlabel" : "Trailing p_{T}(#mu) [GeV]",
		"xlabel_delta_over_truth" : "Trailng #Deltap_{T}(#mu) / p_{truth}(#mu)",
		"bins_delta_over_truth" : [2, 4, 6, 8, 10, 100],
		"min" : 1,
	},
	"dr" : {
		"nbins" : 25,
		"range" : (0.0, 0.3),
		"nbins2" : 50,
		"xlabel" : "#DeltaR(#mu#mu)",
		"xlabel_delta_over_truth" : "#Delta#DeltaR(#mu#mu) / #DeltaR_{truth}(#mu#mu)",
		"bins_delta_over_truth" : [0.05, 0.1, 0.15, 0.2, 0.3, 100],
		"min" : 1,
	}
}

colors = (
	ROOT.kOrange,
	ROOT.kAzure,
	ROOT.kRed,
	ROOT.kCyan,
	ROOT.kGreen,
	ROOT.kViolet,
	ROOT.kPink,
	ROOT.kTeal,
	ROOT.kMagenta,
	ROOT.kSpring,
	ROOT.kBlue,
	ROOT.kYellow,
)

decays = {
	0: {
		"color" : ROOT.kGray,
		"label" : "Other",
	},
	1: {
		"color" : ROOT.kOrange,
		"label" : "#eta#rightarrow#mu#mu",
	},
	2: {
		"color" : ROOT.kGreen+3,
		"label" : "#eta#rightarrow#mu#mu#gamma",
	},
	4: {
		"color" : ROOT.kOrange-3,
		"label" : "#eta'#rightarrow#mu#mu#gamma",
	},
	5: {
		"color" : ROOT.kViolet,
		"label" : "#omega#rightarrow#mu#mu",
	},
	6: {
		"color" : ROOT.kSpring,
		"label" : "#omega#rightarrow#mu#mu#gamma",
	},
	7: {
		"color" : ROOT.kCyan,
		"label" : "#omega#rightarrow#pi_{0}#mu#mu",
	},
	8: {
		"color" : ROOT.kAzure,
		"label" : "#rho#rightarrow#mu#mu",
	},
	9: {
		"color" : ROOT.kYellow+2,
		"label" : "#rho#rightarrow#mu#mu#gamma",
	},
	10: {
		"color" : ROOT.kMagenta,
		"label" : "#phi#rightarrow#mu#mu",
	},
	11: {
		"color" : ROOT.kTeal,
		"label" : "#phi#rightarrow#mu#mu#gamma",
	},
}

decays = {
	1: {
		"color" : ROOT.kRed,
		"label" : "#eta#rightarrow#mu#mu",
	},
	2: {
		"color" : ROOT.kGreen,
		"label" : "#eta#rightarrow#mu#mu#gamma",
	},
	5: {
		"color" : ROOT.kViolet,
		"label" : "#omega#rightarrow#mu#mu",
	},
	7: {
		"color" : ROOT.kCyan,
		"label" : "#omega#rightarrow#pi_{0}#mu#mu",
	},
}

histos = {}
for v in vars_to_plot :
	histos[v] = {}

	for h in histos_to_plot :
		name = h + "_" + v
		if h == "delta_over_truth" :
			h_range = ( -0.1, 0.1 )
			nbins = vars_to_plot[v]["nbins2"]
		else :
			h_range = ( vars_to_plot[v]["range"][0], vars_to_plot[v]["range"][1] )
			nbins = vars_to_plot[v]["nbins"]
		histos[v][h] = {
			"total" : ROOT.TH1F(
				name,
				name,
				nbins,
				h_range[0],
				h_range[1],
			)
		}
		for d in decays :
			histos[v][h]["decay_" + str(d)] = ROOT.TH1F(
				name + "_decay_" + str(d),
				name + "_decay_" + str(d),
				nbins,
				h_range[0],
				h_range[1],
			)
			histos[v][h]["decay_" + str(d)].SetLineColor(decays[d]["color"])
			histos[v][h]["decay_" + str(d)].SetLineWidth(2)
			if h == "efficiency" :
				itree.Draw("gen_matched" + "_" + v + ">>" + name + "_decay_" + str(d), "decay==" + str(d))
			else :
				itree.Draw(name + ">>" + name + "_decay_" + str(d), "decay==" + str(d))
			histos[v][h]["decay_" + str(d)].Sumw2()
		for b in range(len(vars_to_plot[v]["bins_delta_over_truth"])-1) :
			histos[v][h][b] = ROOT.TH1F(
				name + "_" + str(b),
				name + "_" + str(b),
				nbins,
				h_range[0],
				h_range[1],
			)
			histos[v][h][b].SetLineColor(colors[b])
			histos[v][h][b].SetLineWidth(2)
			histos[v][h][b].Sumw2()
			condition = str(vars_to_plot[v]["bins_delta_over_truth"][b]) + "<" + "gen_matched" + "_" + v + " && " + str(vars_to_plot[v]["bins_delta_over_truth"][b+1]) + ">" + "gen_matched" + "_" + v
			if h == "efficiency" :
				itree.Draw("gen_matched" + "_" + v + ">>" + name + "_" + str(b), condition)
			else : 
				#print(name + ">>" + name + "_" + str(b), condition)
				itree.Draw(name + ">>" + name + "_" + str(b), condition)
			histos[v][h][b].Sumw2()
		if h == "efficiency" :
			itree.Draw("gen_matched" + "_" + v + ">>" + name)
			histos[v][h]["total"].Sumw2()
			for p in histos[v][h] :
				histos[v][h][p].Divide(histos[v]["gen"][p])
		else :
			itree.Draw(name + ">>" + name)
			histos[v][h]["total"].Sumw2()
		histos[v][h]["total"].SetLineColor(ROOT.kBlack)
		histos[v][h]["total"].SetLineWidth(2)
		histos[v][h]["total"].GetXaxis().SetLabelOffset(0.02)
		histos[v][h]["total"].GetXaxis().SetTitleOffset(1.9)
		histos[v][h]["total"].GetYaxis().SetTitleOffset(1.5)
		histos[v][h]["total"].GetYaxis().SetTitle(histos_to_plot[h]["ylabel"])
		if h == "delta_over_truth" : 
			histos[v][h]["total"].GetXaxis().SetTitle(vars_to_plot[v]["xlabel_delta_over_truth"])
		else :
			histos[v][h]["total"].GetXaxis().SetTitle(vars_to_plot[v]["xlabel"])
		

		### Decay plots
		if h == "gen" or h == "gen_matched" or h == "reco_matched" or h == "efficiency" :
			c_name = name + "_decays"
			c = ROOT.TCanvas("c_" + c_name, "c_" + c_name, 1200, 1200)
			c.cd()
			if h == "efficiency" :
				histos[v][h]["total"].SetMinimum(0)
				histos[v][h]["total"].SetMaximum(1.0)
			else :
				c.SetLogy()
				histos[v][h]["total"].SetMinimum(vars_to_plot[v]["min"])
			c.SetTopMargin(0.1)
			c.SetLeftMargin(0.12)
			c.SetRightMargin(0.08)
			c.SetBottomMargin(0.17)
			legend = ROOT.TLegend (0.68, 0.73, .92, .89)
			legend.SetMargin(0.5)
			legend.SetNColumns(2)
			legend.SetTextSize (0.015)
			#legend.SetHeader("#bf{Decay}")
			histos[v][h]["total"].Draw("same")
			legend.AddEntry (histos[v][h]["total"], "Total", "f")
			for d in decays :
				histos[v][h]["decay_" + str(d)].Draw("same")
				legend.AddEntry (histos[v][h]["decay_" + str(d)], decays[d]["label"], "f")
			legend.SetLineWidth (0)
			legend.Draw("same")
			c.Update()
			c.SaveAs("Plots/"+c_name+".png")

		### Bin plot
		if h == "delta_over_truth" :
			c_name = name + "_bins"
			c = ROOT.TCanvas("c_" + c_name, "c_" + c_name, 1200, 1200)
			c.cd()
			#if v != "efficiency" : c.SetLogy()
			c.SetTopMargin(0.1)
			c.SetLeftMargin(0.12)
			c.SetRightMargin(0.08)
			c.SetBottomMargin(0.17)
			legend = ROOT.TLegend (0.65, 0.67, .92, .89)
			legend.SetTextSize (0.015)
			histos[v][h]["total"].Sumw2()
			if histos[v][h]["total"].Integral() > 0 : histos[v][h]["total"].Scale(1/histos[v][h]["total"].Integral())
			histos[v][h]["total"].SetMaximum(1.25*histos[v][h]["total"].GetMaximum())
			histos[v][h]["total"].Draw("same")
			legend.AddEntry (histos[v][h]["total"], "Total", "L")
			hh = []
			for b in range(len(vars_to_plot[v]["bins_delta_over_truth"])-1) :
				bin_width = (vars_to_plot[v]["range"][1]-vars_to_plot[v]["range"][0]) / vars_to_plot[v]["nbins2"]
				label = str(vars_to_plot[v]["bins_delta_over_truth"][b]) + "<" + vars_to_plot[v]["xlabel"] + "<" + str(vars_to_plot[v]["bins_delta_over_truth"][b+1])
				histos[v][h][b].Sumw2()
				if histos[v][h][b].Integral() > 0 : histos[v][h][b].Scale(1/histos[v][h][b].Integral())
				histos[v][h][b].Draw("same")
				legend.AddEntry (histos[v][h][b], label, "L")
			legend.SetLineWidth (0)
			legend.Draw("same")
			c.Update()
			c.SaveAs("Plots/"+c_name+".png")

