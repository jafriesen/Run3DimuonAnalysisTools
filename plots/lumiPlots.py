import ROOT
import math
from array import array
from pathlib import Path

import utils.CMSGraphics as CMSGraphics

#ROOT.gStyle.SetTitleOffset(2.5, "X")
#ROOT.gStyle.SetTitleOffset(1.8, "Y")
ROOT.gROOT.ForceStyle()

lumidata = Path("data/lumidata.txt").read_text().splitlines()

run_range = ( 355794, 372415 )

configs = {
	"Run2022C" : {
		"color" : ROOT.kRed,
	},
	"Run2022D" : {
		"color" : ROOT.kMagenta,
	},
	"Run2022E" : {
		"color" : ROOT.kBlue,
	},
	"Run2022F" : {
		"color" : ROOT.kCyan,
	},
	"Run2022G" : {
		"color" : ROOT.kGreen+4,
	},
	"Run2023B" : {
		"color" : ROOT.kSpring,
	},
	"Run2023C" : {
		"color" : ROOT.kYellow-3,
	},
	"Run2023D" : {
		"color" : ROOT.kOrange,
	},
}

run_range = ( 355800, 370800 )
fill_range = ( 8000, 9100 )
for c in configs :
	configs[c]["histo_runs"] = ROOT.TH1I( c + "_runs", c + "_runs", run_range[1]-run_range[0], run_range[0], run_range[1] )
	configs[c]["histo_fills"] = ROOT.TH1I( c + "_fills", c + "_fills", fill_range[1]-fill_range[0], fill_range[0], fill_range[1] )
	configs[c]["histo_runs2"] = ROOT.TH1F( c + "_runs2", c + "_runs2", 100, 0, 1 )
	configs[c]["histo_fills2"] = ROOT.TH1F( c + "_fills2", c + "_fills2", 100, 0, 1 )

histo_runs2 = ROOT.TH1F( "histo_runs2", "histo_runs2", 60, 0, 0.9 )
histo_fills2 = ROOT.TH1F( "histo_fills2", "histo_fills2", 60, 0, 1.2 )

fill = 0
lumi_in_fill = 0
for l in lumidata :
	if l[0:3] == "Run" :
		current = l
	else :
		data = l.split(",")
		if fill != int(data[0].split(":")[1]) :
			print(fill, lumi_in_fill)
			histo_fills2.Fill( lumi_in_fill )
			lumi_in_fill = 0
			fill = int(data[0].split(":")[1])
		print(current,data[0].split(":")[0],data[0].split(":")[1],data[-1])
		configs[current]["histo_runs"].Fill( int(data[0].split(":")[0]), float(data[-1])/10**9 )
		configs[current]["histo_fills"].Fill( int(data[0].split(":")[1]), float(data[-1])/10**9 )
		configs[current]["histo_runs2"].Fill( float(data[-1])/10**9 )
		lumi_in_fill += float(data[-1])/10**9
		histo_runs2.Fill( float(data[-1])/10**9 )

canvas = CMSGraphics.makeCMSCanvas("runs","runs",1400,1200)
canvas.SetRightMargin(0.08)
canvas.SetBottomMargin(0.13)
canvas.SetLeftMargin(0.12)
#canvas.SetLogy()
#canvas.SetLeftMargin(0.15)
#canvas.SetBottomMargin(0.14)
canvas.cd()
for c in configs :
	configs[c]["histo_runs"].GetYaxis().SetTitle("Integrated Luminosity / fb / 100 Runs")
	configs[c]["histo_runs"].GetXaxis().SetTitle("Run")
	configs[c]["histo_runs"].Rebin(100)
	configs[c]["histo_runs"].SetMaximum(configs[c]["histo_runs"].GetMaximum()*2.75)
	configs[c]["histo_runs"].SetLineColor(configs[c]["color"])
	configs[c]["histo_runs"].SetFillColorAlpha(configs[c]["color"],0.5)
	configs[c]["histo_runs"].SetFillStyle(1001)
	configs[c]["histo_runs"].Draw("hist same")
CMSGraphics.printLumiLeft(canvas, extraText="Preliminary", lumitext="62.6 fb^{-1} (13.6 TeV)")
canvas.Update()
canvas.SaveAs("plots/lumiPerRun.png")

canvas = CMSGraphics.makeCMSCanvas("fills","fills",1400,1200)
canvas.SetRightMargin(0.08)
canvas.SetBottomMargin(0.13)
canvas.SetLeftMargin(0.12)
#canvas.SetLeftMargin(0.15)
#canvas.SetBottomMargin(0.14)
canvas.cd()
for c in configs :
	configs[c]["histo_fills"].GetYaxis().SetTitle("Integrated Luminosity / fb / 10 Fills")
	configs[c]["histo_fills"].GetXaxis().SetTitle("Fill")
	configs[c]["histo_fills"].Rebin(10)
	configs[c]["histo_fills"].SetMaximum(configs[c]["histo_fills"].GetMaximum()*3.5)
	configs[c]["histo_fills"].SetLineColor(configs[c]["color"])
	configs[c]["histo_fills"].SetFillColorAlpha(configs[c]["color"],0.5)
	configs[c]["histo_fills"].SetFillStyle(1001)
	configs[c]["histo_fills"].Draw("hist same")
CMSGraphics.printLumiLeft(canvas, extraText="Preliminary", lumitext="62.6 fb^{-1} (13.6 TeV)")
canvas.Update()
canvas.SaveAs("plots/lumiPerFill.png")

canvas = CMSGraphics.makeCMSCanvas("runs2","runs2",1400,1200)
canvas.SetRightMargin(0.08)
canvas.SetBottomMargin(0.13)
canvas.SetLeftMargin(0.12)
#canvas.SetLogy()
#canvas.SetLeftMargin(0.15)
#canvas.SetBottomMargin(0.14)
canvas.cd()
histo_runs2.GetYaxis().SetTitle("Runs")
histo_runs2.GetXaxis().SetTitle("Integrated Luminosity [/fb]")
#configs[c]["histo_runs2"].Rebin(100)
#configs[c]["histo_runs2"].SetMaximum(configs[c]["histo_runs2"].GetMaximum()*2.75)
histo_runs2.SetLineColor(ROOT.kAzure-9)
histo_runs2.SetFillColorAlpha(ROOT.kAzure-9,0.75)
histo_runs2.SetFillStyle(1001)
histo_runs2.Draw("hist same")
CMSGraphics.printLumiLeft(canvas, extraText="Preliminary", lumitext="62.6 fb^{-1} (13.6 TeV)")
canvas.Update()
canvas.SaveAs("plots/lumiPerRun2.png")

canvas = CMSGraphics.makeCMSCanvas("fills2","fills2",1400,1200)
canvas.SetRightMargin(0.08)
canvas.SetBottomMargin(0.13)
canvas.SetLeftMargin(0.12)
#canvas.SetLeftMargin(0.15)
#canvas.SetBottomMargin(0.14)
canvas.cd()
histo_fills2.GetYaxis().SetTitle("Fills")
histo_fills2.GetXaxis().SetTitle("Integrated Luminosity [/fb]")
#configs[c]["histo_fills2"].Rebin(10)
#configs[c]["histo_fills2"].SetMaximum(configs[c]["histo_fills2"].GetMaximum()*3.5)
histo_fills2.SetLineColor(ROOT.kAzure-9)
histo_fills2.SetFillColorAlpha(ROOT.kAzure-9,0.75)
histo_fills2.SetFillStyle(1001)
histo_fills2.Draw("hist same")
CMSGraphics.printLumiLeft(canvas, extraText="Preliminary", lumitext="62.6 fb^{-1} (13.6 TeV)")
canvas.Update()
canvas.SaveAs("plots/lumiPerFill2.png")

