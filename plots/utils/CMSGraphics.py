import os, sys, imp
import ROOT
from array import array
import utils.CMS_lumi as CMS_lumi

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelOffset(0.02, "X")
ROOT.gStyle.SetTitleOffset(1.5, "X")
ROOT.gStyle.SetTitleOffset(1.3, "Y")
ROOT.gStyle.SetTitleSize(0.04, "XY")
ROOT.gROOT.ForceStyle()

T = 0.08
B = 0.1
L = 0.18
R = 0.04

def makeCMSCanvas(name="canvas",title="canvas",width=1800,height=1200):
    canvas = ROOT.TCanvas(name,name,50,50,width,height)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin(L * height / width)
    canvas.SetRightMargin(R)
    canvas.SetTopMargin(T)
    canvas.SetBottomMargin(B * width / height)
    if width == 800: ROOT.gStyle.SetTitleYOffset(1)
    return canvas

def makeLegend(nentries=1, left=0.75, right=0.98, margin=0.2, scale=1, top=0):
    height = nentries*0.06*scale
    leg = ROOT.TLegend(left,(top + 0.88 - height), right - R, top + 0.98 - T)
    leg.SetMargin(margin)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)
    return leg

def printLumiPrelLeft(canvas, lumitext="13.6 TeV", iPeriod=0):
    #change the CMS_lumi variables (see CMS_lumi.py)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 10 # inside frame left, default
    
    if ( iPos==0 ): CMS_lumi.relPosX = 0.12
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)


def printLumiPrelOut(canvas, lumitext="13 TeV", iPeriod=0):
    #change the CMS_lumi variables (see CMS_lumi.py)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = ""
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 0 # outside frame left
    
    if ( iPos==0 ): CMS_lumi.relPosX = 0.12
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)

def printLumiLeft(canvas, lumitext="13.6 TeV", iPeriod=0, extraText=""):
    #change the CMS_lumi variables (see CMS_lumi.py)
    if not extraText : CMS_lumi.writeExtraText = 0
    CMS_lumi.extraText = extraText
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 10 # inside frame left, default
    
    if ( iPos==0 ): CMS_lumi.relPosX = 0.12
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)


def printLumiOut(canvas, lumitext="13.6 TeV", iPeriod=0, extraText="", relPosX=0.12):
    #change the CMS_lumi variables (see CMS_lumi.py)
    if not extraText : CMS_lumi.writeExtraText = 0
    CMS_lumi.extraText = extraText
    CMS_lumi.lumi_sqrtS = lumitext # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    
    iPos = 0 # outside frame left
    
    if ( iPos==0 ): CMS_lumi.relPosX = relPosX
    
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)

