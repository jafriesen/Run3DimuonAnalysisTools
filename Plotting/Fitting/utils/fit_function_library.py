from typing import OrderedDict, Type
from webbrowser import get

from numpy import isin
import ROOT as rt
from abc import ABC, abstractmethod
from collections import namedtuple

Param = namedtuple('Param', ['val', 'min', 'max'])
ConstParam = namedtuple('ConstParam', ['val'])

CB_p = namedtuple('CB', ['mcb', 'scb', 'acb', 'ncb'])
DoubleGauss_p = namedtuple('DoubleGauss', ['mg', 'sg1', 'sg2', 'sig1frac'])
Cheb2_p = namedtuple('Cheb2', ['a1', 'a2'])
Cheb3_p = namedtuple('Cheb2', ['a1', 'a2', 'a3'])
Voigt_p = namedtuple('Voigt', ['mv', 'wv', 'sv'])
Cheb4_p = namedtuple('Cheb3', ['a1', 'a2', 'a3', 'a4'])


# def make_RRV(name, title, par):
#     if isinstance(par, Param):
#         var = rt.RooRealVar(name, title, par[0], par[1], par[2])
#     elif isinstance(par, ConstParam):
#         var = rt.RooRealVar(name, title, par[0])
#     else:
#         raise TypeError("Argument needs to be a (Const)Param!")
#     return var

# Declare base fit function class with interface
class FitFunction(ABC):
    def __init__(self, x):
        self.pars = OrderedDict()
        self.args = OrderedDict()
        self.x = x

    @abstractmethod
    def make_function(self):
        pass

    def __call__(self):
        return self.func

    def get_params(self):
        return self.pars
    
    def init_params(self):
        for k, v in self.pars.items():
            self.args[k] = self.make_RRV(k, v)
        self.make_function()

    def set_params(self, *args, **kwargs):
        # Take full Param obj or just the value, or namedtuple wrapper
        for arg in args:
            if isinstance(arg, tuple) and hasattr(arg, '_asdict'):
                for k, v in arg._asdict().items():
                    if k not in self.pars.keys():
                        raise NameError(f'Parameter {k} not initialized!')
                    self.pars[k] = ConstParam(v)
            else:
                raise TypeError('List arguments must be nametuples!')
        for k, v in kwargs.items():
            if k not in self.pars.keys():
                raise NameError(f'Parameter {k} not initialized!')
            if isinstance(v, Param) or isinstance(v, ConstParam):
                self.pars[k] = v
            elif isinstance(v, float):
                self.pars[k].val = v
            else:
                raise TypeError(f'Type of {v} not recognized!')

        for k, v in self.pars.items():
            if isinstance(v, Param):
                self.args[k].setVal(v.val)
                self.args[k].setMin(v.min)
                self.args[k].setMax(v.max)
            elif isinstance(v, ConstParam):
                self.args[k].setVal(v.val)

        # self.get_function()

    def make_RRV(self, name, val):
        if isinstance(val, Param):
            var = rt.RooRealVar(name, name, val.val, val.min, val.max)
        elif isinstance(val, ConstParam):
            var = rt.RooRealVar(name, name, val.val)
        else:
            raise TypeError("Argument needs to be a (Const)Param!")
        return var
    
    def create_args(self):
        for k, v in self.pars.items():
            self.args[k] = self.make_RRV(k, v)

    def get_arg_list(self):
        return list(self.args.values())
            

############# List of functions below

### Background functions 

class Threshold(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['x0'] = ConstParam(0.2113)
        self.pars['alpha'] = Param(0.01, 0.0, 5.0)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.func = rt.RooGenericPdf("bkg", "Threshold", "(@0-@1)**(@2)", rt.RooArgList(s.x, *s.get_arg_list()))

class ThreshExp(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['x0'] = ConstParam(0.2113)
        self.pars['alpha'] = Param(0.1, 0.0, 5.0)
        self.pars['C'] = Param(2, -10., 10.)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.func = rt.RooGenericPdf("bkg", "ThreshExp", "(@0-@1)**(@2)*(1-TMath::Exp(-@0/@3))", rt.RooArgList(s.x, *s.get_arg_list()))

class ThreshExpDimu(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['x0'] = ConstParam(0.2113)
        self.pars['alpha'] = Param(1.0, 0.001, 5.0)
        self.pars['C'] = Param(1., -10., 10.)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.func = rt.RooGenericPdf("bkg", "ThreshExpDimu", "(@0-@1)**(@2)*(1-TMath::Exp(-@0/@3))", rt.RooArgList(s.x, *s.get_arg_list()))

class ThreshPol(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['x0'] = ConstParam(0.2113)
        self.pars['alpha'] = Param(2.5, 0.1, 5.)
        self.pars['a0'] = Param(8., 1.0, 10.)
        self.pars['a1'] = Param(-5, -10., 10.)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.func = rt.RooGenericPdf("bkg", "ThreshPol", "(@0-@1)**(@2)*(@3 + @0 * @4)", rt.RooArgList(s.x, *s.get_arg_list()))

class Cheb1(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['a1'] = Param(0.987, -10., 10.)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.func = rt.RooChebychev("bkg", "Cheb1", s.x, rt.RooArgList(*s.get_arg_list()))

class Cheb2(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['a1'] = Param(0.987, -10., 10,)
        self.pars['a2'] = Param(0.145, -10., 10.)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.func = rt.RooChebychev("bkg", "Cheb2", s.x, rt.RooArgList(*s.get_arg_list()))

class Cheb3(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['a1'] = Param(-0.63, -1., 1,)
        self.pars['a2'] = Param(0.92, -1., 1.)
        self.pars['a3'] = Param(0.00009, -0.001, 0.001)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.func = rt.RooChebychev("bkg", "Cheb3", s.x, rt.RooArgList(*s.get_arg_list()))

class Cheb4(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['a1'] = Param(0.0148, -1., 1,)
        self.pars['a2'] = Param(0.0041, -1., 1.)
        self.pars['a3'] = Param(-0.001, -1., 1.)
        self.pars['a4'] = Param(0.001, -1., 1.)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.func = rt.RooChebychev("bkg", "Cheb4", s.x, rt.RooArgList(*s.get_arg_list()))

class Pol1(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['p0'] = Param(5500, 1e3, 1e7,)
        self.pars['p1'] = Param(-1000, -1e4, 0.0)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.func = rt.RooChebychev("bkg", "Pol1", s.x, rt.RooArgList(*s.get_arg_list()))



### Signal functions

class CB(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mcb'] = Param(0.548, 0.35, 0.75)
        self.pars['scb'] = Param(0.00575468, 0.003, 0.1)
        self.pars['acb'] = Param(-0.9, -1.0, -0.5)
        self.pars['ncb'] = ConstParam(10)
        #self.pars['mcb'] = Param(0.548, 0.45, 0.65)
        #self.pars['scb'] = Param(0.006, 0.005, 0.1)
        #self.pars['acb'] = Param(-0.3, -0.9, -0.1)
        #self.pars['ncb'] = ConstParam(100)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.func = rt.RooCBShape("CB", "CB", s.x, *s.get_arg_list())
        s.mean = s.args['mcb']
        # Best-fit from MC:
        # ("mcb", "mcb", 0.5498)
        # ("acb", "acb", -0.947)
        # ("ncb", "ncb", 21.04)
        # ("scb", "scb", 0.0057, 0.005, 0.006)

class Johnson(FitFunction):
    # DOES NOT SEEM TO CONVERGE!
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mu'] = Param(0.550, 0.545, 0.555)
        self.pars['lambda'] = Param(-0.1, -1, 0.5)
        self.pars['gamma'] = Param(0.05, 0.03, 0.08)
        self.pars['delta'] = Param(1., 0.5, 1.5)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.func = rt.RooJohnson("sig", "Johnson", s.x, *s.get_arg_list())
        s.mean = s.args['mu']

class SingleGauss(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mg'] = Param(0.5478, 0.35, 0.75)
        self.pars['sg'] = Param(0.005, 0.001, 0.1)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.func = rt.RooGaussian("sig", "Single Gaussian", s.x, *s.get_arg_list())
        s.mean = s.args['mg']

class DoubleGauss(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mg'] = Param(0.5478, 0.35, 0.75)
        self.pars['sg1'] = Param(0.001, 0.0005, 0.1)
        self.pars['sg2'] = Param(0.005, 0.001, 0.1)
        self.pars['sig1frac'] = Param(0.7, 0.0, 1.0)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.sig1 = rt.RooGaussian("sig1", "Gaussian 1", s.x, s.args['mg'], s.args['sg1'])
        s.sig2 = rt.RooGaussian("sig2", "Gaussian 2", s.x, s.args['mg'], s.args['sg2'])
        s.func = rt.RooAddPdf("sig", "Signal", rt.RooArgList(s.sig1, s.sig2), s.args['sig1frac'])
        s.mean = s.args['mg']

class DoubleGaussDimu(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mg'] = Param(0.5478, 0.5, 0.6)
        self.pars['sg1'] = Param(0.2, 0.01, 0.5)
        self.pars['sg2'] = Param(0.08, 0.003, 0.15)
        self.pars['sig1frac'] = Param(0.8, 0.1, 0.9)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.sig1 = rt.RooGaussian("sig1", "Gaussian 1", s.x, s.args['mg'], s.args['sg1'])
        s.sig2 = rt.RooGaussian("sig2", "Gaussian 2", s.x, s.args['mg'], s.args['sg2'])
        s.func = rt.RooAddPdf("sig", "Signal", rt.RooArgList(s.sig1, s.sig2), s.args['sig1frac'])
        s.mean = s.args['mg']

class TripleGauss(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mg'] = Param(0.5478, 0.5, 0.6)
        self.pars['sg1'] = Param(0.01, 0.001, 0.012)
        self.pars['sg2'] = Param(0.005, 0.001, 0.008)
        self.pars['sg3'] = Param(0.001, 0.001, 0.012)
        self.pars['sig1frac'] = Param(0.8, 0.1, 0.9)
        self.pars['sig2frac'] = Param(0.1, 0.0, 1.0)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.sig1 = rt.RooGaussian("sig1", "Gaussian 1", s.x, s.args['mg'], s.args['sg1'])
        s.sig2 = rt.RooGaussian("sig2", "Gaussian 2", s.x, s.args['mg'], s.args['sg2'])
        s.sig3 = rt.RooGaussian("sig3", "Gaussian 3", s.x, s.args['mg'], s.args['sg3'])
        s.func = rt.RooAddPdf("sig", "Signal", rt.RooArgList(s.sig1, s.sig2, s.sig3), rt.RooArgList(s.args['sig1frac'], s.args['sig2frac']))
        s.mean = s.args['mg']

class Voigtian(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mv'] = Param(0.5478, 0.545, 0.550)
        self.pars['wv'] = Param(0.005, 0.001, 0.008)
        self.pars['sv'] = Param(0.005, 0.001, 0.008)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.func = rt.RooVoigtian("sig", "Voigtian", s.x, *s.get_arg_list())
        s.mean = s.args['mv']

class Landau(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['ml'] = Param(0.5478, 0.545, 0.550)
        self.pars['sl'] = Param(0.005, 0.001, 0.008)
        self.init_params()
        self.set_params(**kwargs)
    
    def make_function(s):
        s.create_args()
        s.func = rt.RooLandau("sig", "Landau", s.x, *s.get_arg_list())
        s.mean = s.args['ml']

class CB_Gauss(FitFunction):
    def __init__(self, x, **kwargs):
        super().__init__(x)
        self.pars['mcb'] = Param(0.548, 0.45, 0.65)
        self.pars['scb'] = Param(0.01, 0.001, 0.025)
        self.pars['acb'] = Param(-0.4, -1.0, -0.2)
        self.pars['ncb'] = Param(80, 1, 1000)
        self.pars['sg'] = Param(0.001, 0.005, 0.08)
        self.pars['CB_frac'] = Param(0.8, 0.5, 0.99)
        self.init_params()
        self.set_params(**kwargs)

    def make_function(s):
        s.create_args()
        s.CB = rt.RooCBShape("CB", "CB", s.x, s.args['mcb'], s.args['scb'], s.args['acb'], s.args['ncb'])
        s.Gauss = rt.RooGaussian("Gauss", "Gauss", s.x, s.args['mcb'], s.args['sg'])
        s.func = rt.RooAddPdf("CB_Gauss", "CB+Gauss", rt.RooArgList(s.CB, s.Gauss), rt.RooArgList(s.args['CB_frac']))
        s.mean = s.args['mcb']

### Map str -> fit function

fit_functions = {
    'SingleGauss': SingleGauss,
    'DoubleGauss': DoubleGauss,
    'DoubleGaussDimuD': DoubleGaussDimu,
    'TripleGauss': TripleGauss,
    'CB': CB,
    'Johnson': Johnson,
    'Voigtian': Voigtian,
    'Landau': Landau,
    'CB_Gauss': CB_Gauss,
    
    'Cheb1': Cheb1,
    'Cheb2': Cheb2,
    'Cheb3': Cheb3,
    'Cheb4': Cheb4,
    'Pol1': Pol1,
    'Threshold': Threshold,
    'ThreshExp': ThreshExp,
    'ThreshExpDimu': ThreshExpDimu,
    'ThreshPol': ThreshPol
}
def get_fit_function(name, var):
    if name not in fit_functions.keys():
        raise NameError('Fit function not known!')
    return fit_functions[name](var)

