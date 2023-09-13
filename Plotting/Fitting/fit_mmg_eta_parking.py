import ROOT as rt
from utils.fitter import fitter_2mu as fitter
import utils.CMSGraphics
import CMS_lumi

rt.gROOT.SetBatch()
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)

bin_width = 0.001

m_range = [round(0.45/bin_width)*bin_width,round(0.65/bin_width)*bin_width]
n_mass_bins = int((m_range[1]-m_range[0]) / bin_width)
print (m_range, n_mass_bins)
f = rt.TFile.Open("Data/histos_MMG.root")
h = f.pfCandPhotons_minDr_massMMG_all

#f = rt.TFile.Open("histos_MMG_mass_parking.root")
#h = f.mass_mmg_minDR
h.Rebin(int(bin_width/0.001))

# Declare observable mass
mass = rt.RooRealVar("mass", "m_{#mu#mu#gamma} [GeV]", m_range[0], m_range[1])

gauss = 1
fit = fitter(mass, bkg_model='ThreshPol', sig_model='CB_Gauss')
if gauss : fit = fitter(mass, bkg_model='ThreshPol', sig_model='CB_Gauss')

# For each pT bin, fit the mass and add the yield and errors to the histos
tot_yields = 0.
yields, errors = {}, {}
massframe = mass.frame(rt.RooFit.Title(f"Signal + background fit"))
dh = rt.RooDataHist("dh", "dh", mass, rt.RooFit.Import(h))
#result = fit.model.fitTo(dh, rt.RooFit.Range(m_range[0],m_range[1]), rt.RooFit.Save())
result = fit.model.fitTo(dh, rt.RooFit.Save())
dh.plotOn(massframe, DrawOption="E", LineColor=rt.kRed, LineWidth=2, MarkerSize=0.5, MarkerColor=rt.kRed, Name="data")
fit.model.plotOn(massframe, LineColor=rt.kBlack, LineWidth=2, Name="model")
fit.model.plotOn(massframe, Components={fit.bkg()}, LineColor=rt.kGreen-1, LineStyle='-', LineWidth=2, Name="bkg")
fit.model.plotOn(massframe, Components={fit.sig()}, LineColor=rt.kViolet, LineStyle='-', LineWidth=2, Name="sig")
if gauss :
	fit.model.plotOn(massframe, Components={fit.sig.CB}, LineColor=rt.kRed, LineStyle='--', LineWidth=2, Name="CB")
	fit.model.plotOn(massframe, Components={fit.sig.Gauss}, LineColor=rt.kAzure, LineStyle='--', LineWidth=2, Name="Gauss")

c1 = rt.TCanvas("c1", "c1", 1200, 800)
c1.SetTopMargin(0.1)
c1.SetLeftMargin(0.12)
c1.SetRightMargin(0.08)
c1.SetBottomMargin(0.17)

massframe.GetXaxis().SetLabelOffset(0.02)
massframe.GetXaxis().SetTitleOffset(1.9)
massframe.GetXaxis().SetTitle("m_{#mu#mu#gamma} [GeV]")
massframe.GetYaxis().SetTitleOffset(1.5)
massframe.GetYaxis().SetTitle("Events / " + str(bin_width*1000) + " MeV")
massframe.SetLineWidth(1)
massframe.Draw()

CMS_lumi.writeExtraText = True
CMS_lumi.cmsTextSize    = 0.6
CMS_lumi.lumiTextSize   = 0.46
CMS_lumi.extraOverCmsTextSize = 0.75
CMS_lumi.relPosX = 0.1
CMS_lumi.lumi_sqrtS = "2022F (13.6 TeV)"

# Extract the signal yield and error from the fit
yld = fit.nsig.getVal()
yields[6] = yld
tot_yields += yld
error = fit.nsig.getError()
errors[6] = error

print(result.Print("v"))

chi2_ndof = massframe.chiSquare(result.floatParsFinal().getSize())
print(result.floatParsFinal().getSize())
print("Chi2/ndof:",chi2_ndof)

print("total signal yield: ", tot_yields)

chi2_sig_model = rt.RooChi2Var("chi2_sig_model","chi2", fit.sig(), dh)
chi2_over_ndof = chi2_sig_model.getVal() / (n_mass_bins - result.floatParsFinal().getSize())
#print(rt.RooChi2Var("chi2_sig_model","chi2", fit.model, dh, rt.RooFit.DataError(rt.RooAbsData.Expected), rt.RooFit.Range(m_range[0],m_range[1])))
print("chi2_sig_model:",chi2_sig_model.getVal(), "chi2_over_ndof:",chi2_over_ndof)

chi2_model = rt.RooChi2Var("chi2_model","chi2", fit.model, dh)
chi2_over_ndof = chi2_model.getVal() / (n_mass_bins - result.floatParsFinal().getSize())
print(rt.RooChi2Var("chi2_sig_model","chi2", fit.model, dh, rt.RooFit.DataError(rt.RooAbsData.Expected), rt.RooFit.Range(m_range[0],m_range[1])))
print("chi2_model:",chi2_model.getVal(), "chi2_over_ndof:",chi2_over_ndof)

legend = rt.TLegend (0.15, 0.6, .3, .89)
legend.SetTextSize (0.03)
legend.SetHeader("#chi^{2}/ndof: " + str(round(chi2_over_ndof,5)) + "   signal yield: " + str(round(tot_yields,5)))
legend.AddEntry (c1.FindObject("data"), "data", "LP")
legend.AddEntry (c1.FindObject("model"), "model", "LP")
legend.AddEntry (c1.FindObject("bkg"), "bkg", "LP")
legend.AddEntry (c1.FindObject("sig"), "sig", "LP")
if gauss :
	legend.AddEntry (c1.FindObject("CB"), "CB", "LP")
	legend.AddEntry (c1.FindObject("Gauss"), "Gauss", "LP")
legend.SetLineWidth (0)

CMS_lumi.CMS_lumi(c1, 0, 0)
legend.Draw("same")
c1.Update()

c1.SaveAs("Plots/p_mass_mmg_parking.png")

