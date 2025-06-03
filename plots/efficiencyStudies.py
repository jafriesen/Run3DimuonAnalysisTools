#!/usr/bin/env python3

import ROOT, math, os
from array import array
import utils.CMSGraphics as CMSGraphics

ETA_MASS = 0.547862

itree = ROOT.TChain()
for i in range(9)[1:]:
	itree.Add("data/MC Efficiency Studies/MCTree_" + str(int(i)) + ".root?#tree_" + str(int(i)))

config_decay = {
	"any" : {
		"decay" : -1,
		"color" : ROOT.kGray+1,
		"label" : "Any",
	},
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

config_var = {
	"mass_mumu" : {
		"nbins" : int((0.3-0.21)/0.005),
		"range" : (0.21, 0.3),
		"label" : "#font[12]{m}_{#it{#mu#mu}} [GeV]",
	},
	"pt_mumu" : {
		"nbins" : 40,
		"range" : (5.0, 25.0),
		"bins" : [ 5 + i * 0.5 for i in range(10) ] + [ 10 + i * 1 for i in range(6) ]  + [ 16 + i * 2 for i in range(6) ],
		"label" : "Dimuon #font[12]{p}_{T} [GeV]",
	},
	"pt_leading" : {
		"nbins" : 50,
		"range" : (2.0, 20),
		"bins" : [ 0 + i * 0.5 for i in range(20) ] + [ 10 + i * 1 for i in range(6) ]  + [ 16 + i * 2 for i in range(6) ] + [30, 34],
		"label" : "Leading #font[12]{p}_{T} [GeV]",
	},
	"pt_trailing" : {
		"nbins" : 40,
		"range" : (2.0, 6),
		"label" : "Trailing #font[12]{p}_{T} [GeV]",
	},
	"dr" : {
		"nbins" : 60,
		"range" : (0.0, 0.4),
		"label" : "Dimuon #DeltaR",
	},
}

config_var = {
	"mass_mumu" : {
		"nbins" : int((0.611-0.211)/0.01),
		"range" : (0.211, 0.611),
		"label" : "#font[12]{m}_{#it{#mu#mu}} [GeV]",
	}
}

width = 0.4/500
bins = [ 0.2 ]
bins += [ bins[-1] + (i + 1) * (0.4 - bins[-1]) / 20 for i in range(20) ]
bins += [ bins[-1] + (i + 1) * (0.5 - bins[-1]) / 5 for i in range(5) ]
#bins += [ bins[-1] + (i + 1) * (0.5 - bins[-1]) / 2 for i in range(2) ]
bins += [ 0.55 ]
print(bins)
if False : 
	config_decay["etaToMuMuGamma"]["mass_mumu"] = {"bins" : bins}



gen_match_level = "gen_matched"
bins_pt = [ 6.04, 7.81, 8.18, 8.50, 8.82, 9.14, 9.47, 9.81, 10.17, 10.56, 10.97, 11.43, 11.95, 12.54, 13.24, 14.10, 15.18, 16.67, 18.91, 23.29, 320 ]
for b_low,b_up in zip(bins_pt[:-1],bins_pt[1:]) : 
	for v in config_var :
		plots = ( "gen", gen_match_level )
		for c in config_decay :
			if not v in config_decay[c] :
				config_decay[c][v] = {}
			pt_condition = f"gen_pt_mumu > {b_low} && gen_pt_mumu < {b_up} && gen_pt_trailing > 3" # gen_pt_mumu > 6.04 && gen_pt_mumu < 7.81 
			condition = pt_condition+" && decay=="+str(config_decay[c]["decay"]) if config_decay[c]["decay"] > -1 else pt_condition
			for p in plots :
				name = v+"_"+p+"_"+c
				if "bins" in config_var[v] :
					config_decay[c][v][p] = ROOT.TH1F( name, name, len(config_var[v]["bins"])-1, array('d', config_var[v]["bins"]) )
				elif "bins" in config_decay[c][v] :
					config_decay[c][v][p] = ROOT.TH1F( name, name, len(config_decay[c][v]["bins"])-1, array('d', config_decay[c][v]["bins"]) )
				else :
					config_decay[c][v][p] = ROOT.TH1F( name, name, config_var[v]["nbins"], config_var[v]["range"][0], config_var[v]["range"][1] )
				itree.Draw( p + "_" + v + ">>" + name, condition )
				config_decay[c][v][p].Sumw2()
			config_decay[c][v]["efficiency"] = config_decay[c][v][gen_match_level].Clone("efficiency_"+c)
			config_decay[c][v]["efficiency"].Divide(config_decay[c][v]["gen"])
			print(gen_match_level,c,v,config_decay[c][v][gen_match_level].Integral())
			print("gen",c,v,config_decay[c][v]["gen"].Integral())
	### MASS DIMU PLOTS
	v = "mass_mumu"
	# Dimuon mass efficiency
	#canvas = CMSGraphics.makeCMSCanvas(v,v,1000,2000,L=0.08,R=0.01,T=0.05,B=0.2)
	#canvas.Divide(1,2,0,0)
	canvas = CMSGraphics.makeCMSCanvas(v,v,1400,1200)
	canvas.cd(1)
	legend = CMSGraphics.makeLegend(nentries=4, left=0.8,right=1.0, margin=0.2)
	legend.SetHeader("GEN Decay")
	to_plot = ( "any", "etaToMuMuGamma", "etaToMuMu", "omegaToMuMu", "omegaToPiMuMu" )
	for c in to_plot :
		config_decay[c][v]["efficiency"].SetMaximum(1.3)
		config_decay[c][v]["efficiency"].SetMinimum(0.0)
		config_decay[c][v]["efficiency"].SetLineWidth(3)
		config_decay[c][v]["efficiency"].SetLineColor(config_decay[c]["color"])
		config_decay[c][v]["efficiency"].GetXaxis().SetTitle(config_var[v]["label"])
		config_decay[c][v]["efficiency"].GetYaxis().SetTitle("Efficiency")
		if c=="any" or c=="etaToMuMuGamma" or c=="etaToMuMu" : 
			config_decay[c][v]["efficiency"].SetFillColorAlpha(config_decay[c]["color"],0.5)
			config_decay[c][v]["efficiency"].SetFillStyle(1001)
			config_decay[c][v]["efficiency"].Draw("H same")
		else :
			config_decay[c][v]["efficiency"].Draw("same")
		legend.AddEntry(config_decay[c][v]["efficiency"], config_decay[c]["label"], "f")
	#CMSGraphics.printLumiOut(canvas, extraText="Simulation Preliminary")
	if b_low < 10.17 : legend.Draw()
	canvas.Update()
	'''
	canvas.cd(2)
	legend2 = CMSGraphics.makeLegend(nentries=2, left=0.75, scale=1.5)
	legend2.SetHeader("MC Events")
	config_decay["any"][v]["gen"].SetMinimum(0)
	#config_decay["any"][v]["gen"].SetMaximum(8e3)
	config_decay["any"][v]["gen"].SetLineWidth(3)
	config_decay["any"][v]["gen"].SetLineColor(ROOT.kAzure-4)
	config_decay["any"][v]["gen"].SetFillColorAlpha(ROOT.kAzure-4,0.35)
	config_decay["any"][v]["gen"].SetFillStyle(1001)
	config_decay["any"][v]["gen"].GetXaxis().SetTitle("#font[12]{m}_{#it{#mu#mu}} [GeV]")
	config_decay["any"][v]["gen"].GetYaxis().SetTitle("Events")
	config_decay["any"][v]["gen"].Draw("H same")
	legend2.AddEntry(config_decay["any"][v]["gen"], "All GEN-level", "f")
	config_decay["any"][v][gen_match_level].SetLineWidth(3)
	config_decay["any"][v][gen_match_level].SetLineColor(ROOT.kViolet)
	config_decay["any"][v][gen_match_level].SetFillColorAlpha(ROOT.kViolet,0.35)
	config_decay["any"][v][gen_match_level].SetFillStyle(1001)
	config_decay["any"][v][gen_match_level].Draw("H same")
	legend2.AddEntry(config_decay["any"][v][gen_match_level], "Partially Truthmatched", "f")
	#CMSGraphics.printLumiPrelLeft(canvas)
	legend2.Draw("same")
	'''

	canvas.Update()
	canvas.SaveAs("plots/plot_pt"+str(b_low).replace(".","p")+"to"+str(b_up).replace(".","p")+".png")







for v in config_var :
	plots = ( "gen", "gen_matched" )
	for c in config_decay :
		if not v in config_decay[c] :
			config_decay[c][v] = {}
		pt_condition = "gen_pt_mumu > 6.04 && gen_pt_trailing > 3 && gen_pt_leading > 3" #  gen_pt_mumu > 6.04 && gen_pt_mumu < 7.81 
		condition = pt_condition+" && decay=="+str(config_decay[c]["decay"]) if config_decay[c]["decay"] > -1 else pt_condition
		for p in plots :
			name = v+"_"+p+"_"+c
			if "bins" in config_var[v] :
				config_decay[c][v][p] = ROOT.TH1F( name, name, len(config_var[v]["bins"])-1, array('d', config_var[v]["bins"]) )
			elif "bins" in config_decay[c][v] :
				config_decay[c][v][p] = ROOT.TH1F( name, name, len(config_decay[c][v]["bins"])-1, array('d', config_decay[c][v]["bins"]) )
			else :
				config_decay[c][v][p] = ROOT.TH1F( name, name, config_var[v]["nbins"], config_var[v]["range"][0], config_var[v]["range"][1] )
			itree.Draw( p + "_" + v + ">>" + name, condition )
			config_decay[c][v][p].Sumw2()
		config_decay[c][v]["efficiency"] = config_decay[c][v]["gen_matched"].Clone("efficiency_"+c)
		config_decay[c][v]["efficiency"].Divide(config_decay[c][v]["gen"])
		print("gen_matched",c,v,config_decay[c][v]["gen_matched"].Integral())
		print("gen",c,v,config_decay[c][v]["gen"].Integral())

### MASS DIMU PLOTS
v = "mass_mumu"
# Dimuon mass efficiency
canvas = CMSGraphics.makeCMSCanvas(v,v,1400,1600)
canvas.SetLeftMargin(0.12)
canvas.SetBottomMargin(0.13)
canvas.cd()
legend = CMSGraphics.makeLegend(nentries=3, left=0.78, margin=0.4)
legend.SetHeader("GEN Decay")
to_plot = [ "any", "etaToMuMuGamma", "etaToMuMu" ]
for c in to_plot :
	config_decay[c][v]["efficiency"].SetMaximum(1.0)
	config_decay[c][v]["efficiency"].SetMinimum(0.0)
	config_decay[c][v]["efficiency"].SetLineWidth(3)
	config_decay[c][v]["efficiency"].SetLineColor(config_decay[c]["color"])
	config_decay[c][v]["efficiency"].GetXaxis().SetTitle(config_var[v]["label"])
	config_decay[c][v]["efficiency"].GetYaxis().SetTitle("Efficiency")
	if c=="any" or c=="etaToMuMuGamma" or c=="etaToMuMu" : 
		config_decay[c][v]["efficiency"].SetFillColorAlpha(config_decay[c]["color"],0.5)
		config_decay[c][v]["efficiency"].SetFillStyle(1001)
		config_decay[c][v]["efficiency"].Draw("H same")
	else :
		config_decay[c][v]["efficiency"].Draw("same")
	legend.AddEntry(config_decay[c][v]["efficiency"], config_decay[c]["label"], "f")
CMSGraphics.printLumiOut(canvas, extraText="Simulation Preliminary")
legend.Draw("same")
canvas.Update()
canvas.SaveAs("plots/efficiency_mass_dimu_ptmin6p04.png")


# etaToMuMuGamma Full GEN vs Truthmatched GEN
canvas = CMSGraphics.makeCMSCanvas("etaToMuMuGamma","etaToMuMuGamma",1800,1200)
canvas.cd()
#canvas.SetLogy()
legend = CMSGraphics.makeLegend(nentries=2, left=0.75, scale=1.5)
legend.SetHeader("#eta#rightarrow#mu#mu#gamma Events")
config_decay["etaToMuMuGamma"][v]["gen"].SetMinimum(0)
#config_decay["etaToMuMuGamma"][v]["gen"].SetMaximum(5 * 10 ** 5)
#config_decay["etaToMuMuGamma"][v]["gen"].SetAxisRange(0.2, 0.29)
config_decay["etaToMuMuGamma"][v]["gen"].SetLineWidth(3)
config_decay["etaToMuMuGamma"][v]["gen"].SetLineColor(ROOT.kAzure-4)
config_decay["etaToMuMuGamma"][v]["gen"].SetFillColorAlpha(ROOT.kAzure-4,0.35)
config_decay["etaToMuMuGamma"][v]["gen"].SetFillStyle(1001)
config_decay["etaToMuMuGamma"][v]["gen"].GetXaxis().SetTitle("#font[12]{m}_{#it{#mu#mu}} [GeV]")
config_decay["etaToMuMuGamma"][v]["gen"].GetYaxis().SetTitle("Events")
config_decay["etaToMuMuGamma"][v]["gen"].Draw("H same")
legend.AddEntry(config_decay["etaToMuMuGamma"][v]["gen"], "All GEN-level", "f")
config_decay["etaToMuMuGamma"][v]["gen_matched"].SetLineWidth(3)
config_decay["etaToMuMuGamma"][v]["gen_matched"].SetLineColor(ROOT.kViolet)
config_decay["etaToMuMuGamma"][v]["gen_matched"].SetFillColorAlpha(ROOT.kViolet,0.35)
config_decay["etaToMuMuGamma"][v]["gen_matched"].SetFillStyle(1001)
config_decay["etaToMuMuGamma"][v]["gen_matched"].Draw("H same")
legend.AddEntry(config_decay["etaToMuMuGamma"][v]["gen_matched"], "Truthmatched", "f")
CMSGraphics.printLumiPrelLeft(canvas)
legend.Draw("same")
canvas.Update()
canvas.SaveAs("plots/etaToMuMuGammaGEN.png")





### PT PLOTS
for v in ["pt_mumu", "pt_leading", "pt_trailing"] : 
	canvas = CMSGraphics.makeCMSCanvas(v,v,1400,1200)
	canvas.SetLeftMargin(0.12)
	canvas.SetBottomMargin(0.13)
	canvas.cd()
	legend = CMSGraphics.makeLegend(nentries=5, left=0.15, right=0.35, margin=0.4)
	legend.SetHeader("GEN Decay")
	to_plot = ( "any", "etaToMuMuGamma", "etaToMuMu", "omegaToMuMu", "omegaToPiMuMu" )
	for c in to_plot :
		config_decay[c][v]["efficiency"].SetMaximum(1.0)
		config_decay[c][v]["efficiency"].SetMinimum(0)
		config_decay[c][v]["efficiency"].SetLineWidth(3)
		config_decay[c][v]["efficiency"].SetLineColor(config_decay[c]["color"])
		config_decay[c][v]["efficiency"].GetXaxis().SetTitle(config_var[v]["label"])
		config_decay[c][v]["efficiency"].GetYaxis().SetTitle("Efficiency")
		if c == "any" : 
			config_decay[c][v]["efficiency"].SetFillColorAlpha(config_decay[c]["color"],0.5)
			config_decay[c][v]["efficiency"].SetFillStyle(1001)
			config_decay[c][v]["efficiency"].Draw("H same")
		else :
			config_decay[c][v]["efficiency"].Draw("same")
		legend.AddEntry(config_decay[c][v]["efficiency"], config_decay[c]["label"], "f")
	CMSGraphics.printLumiOut(canvas, extraText="Simulation Preliminary")
	legend.Draw("same")
	canvas.Update()
	canvas.SaveAs("plots/efficiency_" + v + ".png")



### DR PLOTS
v = "dr"
canvas = CMSGraphics.makeCMSCanvas(v,v,1400,1200)
canvas.SetLeftMargin(0.12)
canvas.SetBottomMargin(0.13)
canvas.cd()
legend = CMSGraphics.makeLegend(nentries=5, left=0.8, margin=0.4)
legend.SetHeader("GEN Decay")
to_plot = ( "etaToMuMu", "any", "etaToMuMuGamma", "omegaToMuMu", "omegaToPiMuMu" )
for c in to_plot :
	config_decay[c][v]["efficiency"].SetMaximum(1.0)
	config_decay[c][v]["efficiency"].SetMinimum(0)
	config_decay[c][v]["efficiency"].SetLineWidth(3)
	config_decay[c][v]["efficiency"].SetLineColor(config_decay[c]["color"])
	config_decay[c][v]["efficiency"].GetXaxis().SetTitle(config_var[v]["label"])
	config_decay[c][v]["efficiency"].GetYaxis().SetTitle("Efficiency")
	if c=="any" or c=="etaToMuMuGamma" or c=="etaToMuMu" : 
		config_decay[c][v]["efficiency"].SetFillColorAlpha(config_decay[c]["color"],0.5)
		if c=="etaToMuMu" :
			config_decay[c][v]["efficiency"].SetFillColorAlpha(config_decay[c]["color"],0.35) 
			config_decay[c][v]["efficiency"].SetFillStyle(3011)
		config_decay[c][v]["efficiency"].Draw("H same")
	else :
		config_decay[c][v]["efficiency"].Draw("same")
for c in ( "any", "etaToMuMuGamma", "etaToMuMu", "omegaToMuMu", "omegaToPiMuMu" ) :
	legend.AddEntry(config_decay[c][v]["efficiency"], config_decay[c]["label"], "f")
CMSGraphics.printLumiOut(canvas, extraText="Simulation Preliminary")
legend.Draw("same")
canvas.Update()
canvas.SaveAs("plots/efficiency_" + v + ".png")

