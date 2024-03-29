#!/usr/bin/env python3

import ROOT, array, math, os
import utils.CMSGraphics as CMSGraphics

ERA = "62.4 fb^{-1} (13.6 TeV)"

file_input = ROOT.TFile.Open("data/histoMMGFullRun3.root","READ")

canvas = CMSGraphics.makeCMSCanvas("canvas","canvas",1800,1200)
h_mass_dimu = file_input.Get("massDimu")
h_mass_dimu.SetAxisRange(0.2,1.4)
h_mass_dimu.GetXaxis().SetTitle("#font[12]{m}_{#it{#mu#mu}} [GeV]")
h_mass_dimu.GetYaxis().SetTitle("Events / MeV")
h_mass_dimu.SetLineWidth(2)
h_mass_dimu.SetLineColor(ROOT.kAzure-4)
h_mass_dimu.SetFillColorAlpha(ROOT.kAzure-9,0.35)
h_mass_dimu.SetFillStyle(3011)
h_mass_dimu.Draw()

legend = CMSGraphics.makeLegend()
legend.AddEntry(h_mass_dimu, "Parking Data", "f")

CMSGraphics.printLumiPrelLeft(canvas, lumitext=ERA)

legend.Draw("same")
canvas.Update()
canvas.SaveAs("template.png")