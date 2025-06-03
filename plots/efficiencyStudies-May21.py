#!/usr/bin/env python3

import ROOT, math, os
from array import array
import utils.CMSGraphics as CMSGraphics
import numpy as np

ETA_MASS = 0.547862

pt_leading_bins =	( 3.05, 4.20, 4.41, 4.62, 4.82, 5.03, 5.24, 5.47, 5.71,  5.99,  6.29,  6.62,  6.99,  7.43,  7.94,  8.57,  9.40,  10.53, 12.28, 15.72, 213.18 )
pt_mumu_bins =		( 0, 6.04, 7.81, 8.18, 8.50, 8.82, 9.14, 9.47, 9.81, 10.17, 10.56, 10.97, 11.43, 11.95, 12.54, 13.24, 14.10, 15.18, 16.67, 18.91, 23.29, 320.68 )

itree = ROOT.TChain()
for i in range(9)[1:]:
	itree.Add("data/MC Efficiency Studies/MCTree_" + str(int(i)) + ".root?#tree_" + str(int(i)))

label = "#font[12]{m}_{#it{#mu#mu}} [GeV]"

binning = [ 0.211 + 0.001 * i for i in range(20) ]
binning = np.arange(0.21,0.6,0.005)
plots = ( "gen", "gen_matched" )

colors = [
	ROOT.kRed,
	ROOT.kBlue,
	ROOT.kViolet,
	ROOT.kGray,
	ROOT.kOrange,
	ROOT.kPink,
	ROOT.kGray,
	ROOT.kAzure,
	ROOT.kOrange
]



config_decay = {
	"etaToMuMuGamma": {
		"decay" : 2,
		"color" : ROOT.kAzure-4,
		"label" : "#eta#rightarrow#mu#mu#gamma",
	},
	"etaToMuMu" : {
		"decay" : 1,
		"color" : ROOT.kRed,
		"label" : "#eta#rightarrow#mu#mu",
	},
	"omegaToPiMuMu": {
		"decay" : 6,
		"color" : ROOT.kGreen+2,
		"label" : "#omega#rightarrow#pi_{0}#mu#mu",
		"mass_mumu" : { "bins" : [ 0.21 + i * 0.04 for i in range(9) ] + [0.58, 0.63, 0.68] },
		"pt_mumu" : { "bins" : [ i * 0.5 for i in range(20) ] + [ 10 + i * 1 for i in range(5) ]  + [ 15 + i * 2 for i in range(5) ] }
	},
	"omegaToMuMu": {
		"decay" : 5,
		"color" : ROOT.kViolet,
		"label" : "#omega#rightarrow#mu#mu",
		"mass_mumu" : { "bins" : [ 0.75 + i * 0.01 for i in range(20) ] },
		"pt_mumu" : { "bins" : [ i * 0.5 for i in range(20) ] + [ 10 + i * 1 for i in range(5) ]  + [ 15 + i * 2 for i in range(5) ] }
	},
}




for decay in [3] :

	binning = np.arange(0.21,0.6,0.005)

	if decay == "etaToMuMu" : 
		binning = np.arange(0.5478,0.5479,0.000001)

	i_color = 0
	plot_pt_bins = [0,1,2,5,10]
	histos_pt_mumu = []

	histos_pt_mumu.append({
		"color" : ROOT.kBlack,
		"label" : f">{pt_mumu_bins[plot_pt_bins[-1]]}",
		"histo_gen" : ROOT.TH1F( "gen_all", "gen_all", len(binning)-1, array('d', binning) ),
		"histo_gen_matched" : ROOT.TH1F( "gen_matched_all", "gen_matched_all", len(binning)-1, array('d', binning) )
	})
	#d = config_decay[decay]["decay"]
	condition = f"gen_pt_mumu > {pt_mumu_bins[plot_pt_bins[-1]]}"# + f" && decay=={d}"
	print("gen")
	itree.Draw( "gen_mass_mumu >> gen_all", condition )
	histos_pt_mumu[-1]["histo_gen"].Sumw2()
	print("gen_matched")
	itree.Draw( "gen_matched_mass_mumu >> gen_matched_all", condition )
	histos_pt_mumu[-1]["histo_gen_matched"].Sumw2()
	histos_pt_mumu[-1]["histo_efficiency"] = histos_pt_mumu[-1]["histo_gen_matched"].Clone("efficiency_all")
	histos_pt_mumu[-1]["histo_efficiency"].Divide(histos_pt_mumu[-1]["histo_gen"])
	print(histos_pt_mumu[-1]["histo_gen_matched"].Integral(),histos_pt_mumu[-1]["histo_gen"].Integral())


	for ipt in plot_pt_bins :
		print()
		r_str = ( str(pt_mumu_bins[ipt]), str(pt_mumu_bins[ipt+1]) )

		#d = config_decay[decay]["decay"]
		condition = "gen_pt_mumu > " + r_str[0] + " && " + "gen_pt_mumu < " + r_str[1]# + f" && decay=={d}"
		print(condition)

		name = "pt_mumu_" + (r_str[0]).replace(".", "p") + "to" + (r_str[1]).replace(".", "p")
		histos_pt_mumu.append({
			"color" : colors[i_color],
			"label" : r_str[0] + "-" + r_str[1],
			"histo_gen" : ROOT.TH1F( "gen_" + name, "gen_" + name, len(binning)-1, array('d', binning) ),
			"histo_gen_matched" : ROOT.TH1F( "gen_matched_" + name, "gen_matched_" + name, len(binning)-1, array('d', binning) )
		})
		i_color += 1

		print("gen")
		itree.Draw( "gen_mass_mumu >> gen_" + name, condition )
		histos_pt_mumu[-1]["histo_gen"].Sumw2()

		print("gen_matched")
		itree.Draw( "gen_matched_mass_mumu >> gen_matched_" + name, condition )
		histos_pt_mumu[-1]["histo_gen_matched"].Sumw2()
		
		histos_pt_mumu[-1]["histo_efficiency"] = histos_pt_mumu[-1]["histo_gen_matched"].Clone("efficiency_"+name)
		histos_pt_mumu[-1]["histo_efficiency"].Divide(histos_pt_mumu[-1]["histo_gen"])

		print(histos_pt_mumu[-1]["histo_gen_matched"].Integral(),histos_pt_mumu[-1]["histo_gen"].Integral())

	### MASS DIMU PLOTS
	# Dimuon mass efficiency by pt_mumu
	canvas = CMSGraphics.makeCMSCanvas("c","c",1200,1200)
	canvas.SetLeftMargin(0.12)
	canvas.SetBottomMargin(0.13)
	canvas.cd()
	legend = CMSGraphics.makeLegend(nentries=len(plot_pt_bins), left=0.8, margin=0.3, scale=1.0)
	legend.SetTextSize(0.015)
	legend.SetNColumns(1)
	legend.SetHeader("p_{T} [GeV]")
	for h in reversed(histos_pt_mumu) :
		h["histo_efficiency"].SetMaximum(1.0)
		h["histo_efficiency"].SetMinimum(0.0)
		h["histo_efficiency"].SetLineWidth(3)
		h["histo_efficiency"].SetLineColor(h["color"])
		h["histo_efficiency"].SetMarkerStyle(ROOT.kDot)
		h["histo_efficiency"].GetXaxis().SetTitle(label)
		h["histo_efficiency"].GetYaxis().SetTitle("Efficiency")
		h["histo_efficiency"].Draw("E1 same")
		legend.AddEntry(h["histo_efficiency"], h["label"], "l")
	CMSGraphics.printLumiOut(canvas, extraText="Simulation Preliminary")
	legend.Draw("same")
	canvas.Update()
	canvas.SaveAs(f"plots/efficiency_mass_dimu_bins_pt_mumu_{decay}.png")


