#!/usr/bin/env python3

import ROOT, math, os
from array import array
import utils.CMSGraphics as CMSGraphics

ROOT.gStyle.SetTitleOffset(2.5, "X")
ROOT.gStyle.SetTitleOffset(1.8, "Y")
ROOT.gROOT.ForceStyle()

ETA_MASS = 0.547862

ifile = ROOT.TFile.Open("data/histosDimuSelection.root")

muon_IDs = [ "anyID", "isHighPtMuon", "isLooseMuon", "isMediumMuon", "isSoftMuon", "isTightMuon" ]

muon_IDs2 = { "anyID" : "$p_{T} > 3$ GeV", "isHighPtMuon" : "High $p_{T}$", "isLooseMuon" : "Loose", "isMediumMuon" : "Medium", "isSoftMuon" : "Soft", "isTightMuon" :"Tight" }


muon_IDs = [ "anyID", "isLooseMuon", "isMediumMuon", "isSoftMuon" ]
#muon_IDs = [ "anyID" ]
vtxProb_bins = [ 0.005, 0.0075, 0.01, 0.015, 0.02, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
#vtxProb_bins = vtxProb_bins[5:]
#vtxProb_bins = [ 0.5 ]
histos = {}

s1 = {}
s2 = {}

for ID in muon_IDs :
	histos[ID] = {}
	for i_vtxProb in reversed(range(len(vtxProb_bins))) :
		if vtxProb_bins[i_vtxProb] == vtxProb_bins[-1] :
			name = "mumu_mass_" + ID + "_vtxProbMin" + str(vtxProb_bins[i_vtxProb]).replace(".","p")
		else :
			name = "mumu_mass_" + ID + "_vtxProb" + str(vtxProb_bins[i_vtxProb]).replace(".","p") + "to" + str(vtxProb_bins[i_vtxProb+1]).replace(".","p")
		print(name)
		histos[ID][i_vtxProb] = ifile.Get(name)
		h = histos[ID][i_vtxProb]
		h.Sumw2()
		rebinning = 2
		h.Rebin(rebinning)
		if vtxProb_bins[i_vtxProb] != vtxProb_bins[-1] :
			h.Add( histos[ID][i_vtxProb+1] )
		print(histos[ID][i_vtxProb].GetEntries())


		if not vtxProb_bins[i_vtxProb] in [ 0.005, 0.01, 0.025, 0.05, 0.075, 0.1 ] : continue
		bin_width = 0.0001 * rebinning
		m_range = [ round(0.51/bin_width)*bin_width, round(0.59/bin_width)*bin_width ]

		params = [
			(0.5478, 0.545, 0.55), # m
			(0.006, 0.000, 0.1), # s1
			(0.003, 0.000, 0.01), # s2
		]

		sig_params = 6
		bkg_params = 3
		sig_function = "gaus(0) + gaus(3)"
		bkg_function = "cheb{npbkg}({npsig})"
		bkg_function = bkg_function.replace("{npsig}",str(sig_params)).replace("{npbkg}",str(bkg_params-1))

		#peak = h.GetBinContent(h.GetBin(ETA_MASS))
		peak = h.GetMaximum()
		#print(h.GetBin(ETA_MASS), peak)
		print(peak)
		params1 = [
			2e+04,
			ETA_MASS,
			0.007,
			3e+03,
			ETA_MASS,
			0.004,
			5e+05,
			5e+05,
			-6e+05,
			8e+05,
		]
		par = array( 'd', params1 )
		fit_function = ROOT.TF1("fit_function", sig_function + "+" + bkg_function, m_range[0], m_range[1])
		fit_function.SetLineColor(ROOT.kViolet-7)
		fit_function.SetLineWidth(3)
		#fit_function.SetParameters(par)
		fit_function.SetParameter(1, ETA_MASS)
		fit_function.SetParLimits(1, ETA_MASS-0.0005, ETA_MASS+0.0005)
		fit_function.SetParameter(4, ETA_MASS)
		fit_function.SetParLimits(4, ETA_MASS-0.0005, ETA_MASS+0.0005)
		fit_function.SetParameter(2, 0.004)
		fit_function.SetParLimits(2, 0.003,0.005)
		fit_function.SetParameter(5, 0.007)
		fit_function.SetParLimits(5, 0.005,0.008)
		'''
		fit_function.SetParameter(1, params[0][0])
		fit_function.SetParLimits(1, params[0][1], params[0][2])
		fit_function.SetParameter(4, params[0][0])
		fit_function.SetParLimits(4, params[0][1], params[0][2])
		fit_function.SetParameter(2, params[1][0])
		fit_function.SetParLimits(2, params[1][1], params[1][2])
		fit_function.SetParameter(5, params[2][0])
		fit_function.SetParLimits(5, params[2][1], params[2][2])
		'''

		canvas = CMSGraphics.makeCMSCanvas(name,name,1400,1200)
		canvas.SetLeftMargin(0.15)
		canvas.SetBottomMargin(0.14)
		canvas.cd()

		result = h.Fit("fit_function", "RES")
		fit_params = fit_function.GetParameters()
		sig = ROOT.TF1("sig", sig_function, m_range[0], m_range[1])
		bkg = ROOT.TF1("bkg", bkg_function, m_range[0], m_range[1])
		for p in range(bkg_params+sig_params) :
			print(p, fit_params[p])
			sig.SetParameter(p, fit_params[p])
			bkg.SetParameter(p, fit_params[p])

		sig.SetLineColorAlpha(ROOT.kRed+1,0.5)
		sig.SetLineWidth(3)
		sig.Draw("same")
		bkg.SetLineColor(ROOT.kAzure+9)
		bkg.SetLineWidth(3)
		bkg.Draw("same")

		h.SetMinimum(0)
		h.SetMaximum(0.5 * peak)
		h.SetLineWidth(3)
		h.SetMarkerSize(1)
		h.SetLineColor(ROOT.kViolet-7)
		#h.SetFillColorAlpha(ROOT.kViolet-7,0.5)
		#h.SetFillStyle(1001)
		h.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
		h.GetYaxis().SetTitle("Events / " + str(bin_width) + " GeV")
		h.SetAxisRange(m_range[0],m_range[1])
		h.Draw("e same")

		sig_yield = sig.Integral(m_range[0],m_range[1]) / bin_width
		bkg_yield = bkg.Integral(m_range[0],m_range[1]) / bin_width
		print(sig_yield, bkg_yield)
		s1[name] = ( sig_yield, bkg_yield, sig_yield / math.sqrt(bkg_yield), sig_yield / math.sqrt(sig_yield+bkg_yield) )

		legend = CMSGraphics.makeLegend(nentries=4, left=0.7, margin=0.4)
		ID_name = "" if ID == "anyID" else ID + " "
		if vtxProb_bins[i_vtxProb] == vtxProb_bins[-1] :
			header1 = "2 " + ID_name + "#mu, " + str(vtxProb_bins[i_vtxProb]) + "<P_{vtx}"
		else :
			header1 = "2 " + ID_name + "#mu, " + str(vtxProb_bins[i_vtxProb]) + "< P_{vtx}<" + str(vtxProb_bins[i_vtxProb+1])
		#header = "#splitline{" + header1 + "}{Fit #chi^{2}/ndof = " + str( round( result.Chi2()/result.Ndf(), 2 ) ) + "}"
		header = "Fit #chi^{2}/ndof = " + str( round( result.Chi2()/result.Ndf(), 2 ) )
		legend.SetHeader(header)
		#legend.SetHeader("signal events: " + str(round(yld,5)) + "   background events: " + str(round(bkg_events,5)))
		legend.AddEntry(h, "Data", "lep")
		legend.AddEntry(fit_function, "Fit","L")
		legend.AddEntry(sig, "Signal","L")
		legend.AddEntry(bkg, "Background","L")
		
		CMSGraphics.printLumiLeft(canvas, extraText="Preliminary", lumitext="62.6 fb^{-fb} (13.6 TeV)")
		legend.Draw("same")
		canvas.Update()
		canvas.SaveAs("plots/dimuSelectionStudies/" + "mumu_mass_" + ID + "_vtxProbMin" + str(vtxProb_bins[i_vtxProb]).replace(".","p") + ".png")


print(s1)

