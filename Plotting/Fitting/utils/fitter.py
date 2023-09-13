import ROOT as rt
from abc import ABC, abstractmethod
import utils.fit_function_library as library

# Declare base fitter class with interface
class fitter(ABC):
    def __init__(self):
        sig, bkg = None, None
        self.make_signal_pdf()
        self.make_bkg_pdf()
        self.make_model()

    @abstractmethod
    def make_signal_pdf(self):
        pass

    @abstractmethod
    def make_bkg_pdf(self):
        pass

    @abstractmethod
    def make_model(self):
        pass

    def set_sig_params(self, *args, **kwargs):
        self.sig.set_params(*args, **kwargs)

    def set_bkg_params(self, *args, **kwargs):
        self.bkg.set_params(*args, **kwargs)



# Derived class for fitting of the 4-mu mass spectrum
class fitter_4mu(fitter):
    def __init__(self, mass, bkg_model='Cheb2', sig_model='CB'):
        self.mass = mass
        self.bkg_model = bkg_model
        self.sig_model = sig_model
        super().__init__()
    
    def save_workspace(self, data, filename):
        w = rt.RooWorkspace("w", "workspace")
        w.Import(self.bkg)
        w.Import(self.sig)
        w.Import(self.nbkg)
        w.Import(data)
        w.writeToFile(filename)
    
    def make_signal_pdf(self):
        self.sig = library.get_fit_function(self.sig_model, self.mass)

    def make_bkg_pdf(self):
        if self.bkg_model == '':
            return
        self.bkg = library.get_fit_function(self.bkg_model, self.mass)
 
    def make_model(self):
        # Construct fit
        self.nsig = rt.RooRealVar("nsig", "nsig", 50, 1, 10000)
        self.esig = rt.RooExtendPdf("esig", "extended sig pdf", self.sig(), self.nsig)
        if self.bkg_model != '':
            self.nbkg = rt.RooRealVar("bkg_norm", "bkg_norm", 1000, 0, 10000)
            self.ebkg = rt.RooExtendPdf("ebkg", "extended bkg pdf", self.bkg(), self.nbkg)
            self.model = rt.RooAddPdf("s+b", "s+b", rt.RooArgList(self.esig, self.ebkg))
        else:
            self.model = self.esig





# Derived class for fitting of the 2-mu mass spectrum by pT slice
class fitter_2mu(fitter):
    def __init__(self, mass, bkg_model='Cheb3', sig_model='DoubleGauss'):
        self.mass = mass
        self.bkg_model = bkg_model
        self.sig_model = sig_model
        super().__init__()
    
    def mean_val(self):
        return self.sig.mean.getVal()
    
    def mean_err(self):
        return self.sig.mean.getError()

    def make_signal_pdf(self):
        if self.sig_model == '':
            return
        self.sig = library.get_fit_function(self.sig_model, self.mass)
        '''
        if self.sig_model == 'CB':
           self.sig.set_params(
                mcb=library.Param(0.5465, 0.544, 0.550),
                scb=library.Param(0.005, 0.003, 0.015),
                acb=library.Param(-0.01, -10, -0.01),
                ncb=library.Param(10, 5, 50)
            )
        '''

    def make_bkg_pdf(self):
        if self.bkg_model == '':
            return
        self.bkg = library.get_fit_function(self.bkg_model, self.mass)
        if self.bkg_model == 'Cheb2':
            self.bkg.set_params(
                a1=library.Param(0.5, -1., 1.),
                a2=library.Param(0.2, -1., 1.)
            )

    def make_model(self):
        # Sum the composite signal and background into an extended pdf nsig*sig+nbkg*bkg
        self.nsig = rt.RooRealVar("nsig", "number of signal events", 1e5, 0., 1e9)
        self.nbkg = rt.RooRealVar("nbkg", "number of background events", 1e7, 0, 1e9)
        if self.sig_model == '' :
            self.model = self.bkg()
        if self.bkg_model == '' :
            self.model = self.sig()
        else : 
            self.model = rt.RooAddPdf(f"{self.sig_model}_{self.bkg_model}", "sig+bkg", rt.RooArgList(self.bkg(), self.sig()), rt.RooArgList(self.nbkg, self.nsig))
