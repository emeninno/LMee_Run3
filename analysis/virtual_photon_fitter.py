import ROOT
from ROOT import TH1D, TF1, TMath

class VPF(): #virtual photon function
    def __init__(self, h1LF, h1HF, h1vph):
        self.h1LF    = h1LF.Clone("h1LF_tmp");
        self.h1HF    = h1HF.Clone("h1HF_tmp");
        self.h1vph   = h1vph.Clone("h1vph_tmp");
    def __call__(self,x,par):
        r = par[0];#direct photon fraction
        mbinLF  = self.h1LF .FindBin(x[0]);
        mbinHF  = self.h1HF .FindBin(x[0]);
        mbinVPH = self.h1vph.FindBin(x[0]);
        dNdm = r * self.h1vph.GetBinContent(mbinVPH) + (1 - r) * self.h1LF.GetBinContent(mbinLF) + self.h1HF.GetBinContent(mbinHF);
        return dNdm;

