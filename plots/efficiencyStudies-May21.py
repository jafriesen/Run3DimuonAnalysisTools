#!/usr/bin/env python3

import ROOT, math, os
from array import array
import utils.CMSGraphics as CMSGraphics

ETA_MASS = 0.547862

pt_leading_bins =	( 3.05, 4.20, 4.41, 4.62, 4.82, 5.03, 5.24, 5.47, 5.71,  5.99,  6.29,  6.62,  6.99,  7.43,  7.94,  8.57,  9.40,  10.53, 12.28, 15.72, 213.18 )
pt_mumu_bins =		( 6.04, 7.81, 8.18, 8.50, 8.82, 9.14, 9.47, 9.81, 10.17, 10.56, 10.97, 11.43, 11.95, 12.54, 13.24, 14.10, 15.18, 16.67, 18.91, 23.29, 320.68 )

itree = ROOT.TChain()
for i in range(9)[1:]:
	itree.Add("data/MC Efficiency Studies/MCTree_" + str(int(i)) + ".root?#tree_" + str(int(i)))

label = "#font[12]{m}_{#it{#mu#mu}} [GeV]"

bins = [ 0.2 + 0.02 * i for i in range(20) ]
plots = ( "gen", "gen_matched" )

histos_pt_mumu = []
for ipt in range(len(pt_mumu_bins)-1) :
	r_str = ( str(pt_mumu_bins[ipt]), str(pt_mumu_bins[ipt+1]) )
	condition = "gen_pt_mumu > " + r_str[0] + " && " + "gen_pt_mumu < " + r_str[1]
	print(condition)

	name = "pt_mumu_" + (r_str[0]).replace(".", "p") + "to" + (r_str[1]).replace(".", "p")
	histos_pt_mumu.append({
		"color" : ROOT.kAzure - 9 + ipt,
		"label" : r_str[0] + "-" + r_str[1],
		"histo_gen" : ROOT.TH1F( "gen_" + name, "gen_" + name, len(bins)-1, array('d', bins) ),
		"histo_gen_matched" : ROOT.TH1F( "gen_matched_" + name, "gen_matched_" + name, len(bins)-1, array('d', bins) )
	})

	print("gen")
	itree.Draw( "gen_mass_mumu >> gen_" + name, condition )
	histos_pt_mumu[-1]["histo_gen"].Sumw2()

	print("gen_matched")
	itree.Draw( "gen_matched_mass_mumu >> gen_matched_" + name, condition )
	histos_pt_mumu[-1]["histo_gen_matched"].Sumw2()
	
	histos_pt_mumu[-1]["histo_efficiency"] = histos_pt_mumu[-1]["histo_gen_matched"].Clone("efficiency_"+name)
	histos_pt_mumu[-1]["histo_efficiency"].Divide(histos_pt_mumu[-1]["histo_gen"])

	print(histos_pt_mumu[-1]["histo_gen_matched"].GetEntries(),histos_pt_mumu[-1]["histo_gen"].GetEntries())

### MASS DIMU PLOTS
# Dimuon mass efficiency by pt_mumu
canvas = CMSGraphics.makeCMSCanvas("c","c",1200,1200)
canvas.SetLeftMargin(0.12)
canvas.SetBottomMargin(0.13)
canvas.cd()
legend = CMSGraphics.makeLegend(nentries=5, left=0.2, margin=0.4, scale=0.3)
legend.SetTextSize(0.015)
legend.SetNColumns(4)
legend.SetHeader("p_{T} [GeV]")
for h in histos_pt_mumu :
	h["histo_efficiency"].SetMaximum(1.0)
	h["histo_efficiency"].SetMinimum(0.0)
	h["histo_efficiency"].SetLineWidth(3)
	h["histo_efficiency"].SetLineColor(h["color"])
	h["histo_efficiency"].GetXaxis().SetTitle(label)
	h["histo_efficiency"].GetYaxis().SetTitle("Efficiency")
	h["histo_efficiency"].Draw("H same")
	legend.AddEntry(h["histo_efficiency"], h["label"], "l")
CMSGraphics.printLumiOut(canvas, extraText="Simulation Preliminary")
legend.Draw("same")
canvas.Update()
canvas.SaveAs("plots/efficiency_mass_dimu_bins_pt_mumu.png")


