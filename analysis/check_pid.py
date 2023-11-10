import os, sys, shutil
import numpy as np
import ROOT
from ROOT import TFile, TH1D, TH2D, TF1
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange
from file_manager import FileManager
from histo_manager import slice_histogram, rebin_histogram

def check_nsigma_tpc(list_hierarchies, arr_pin, arr_eta, suffix=""):
    print(sys._getframe().f_code.co_name);
    print(list_hierarchies);
    fm = FileManager(list_hierarchies);
    parnames = ["Ele","Pi","Ka","Pr"];
    parsymbols = ["e","#pi","K","p"];
    np = len(parnames);

    outname = "pid_nsigma_TPC{0}.root".format(suffix);
    outfile = TFile.Open(outname,"RECREATE");

    for ip in range(0,np):
        h2_nsigma_pin = fm.get("TPCnSig{0}_pIN".format(parnames[ip]));
        h2_nsigma_pin.Sumw2();
        outfile.WriteTObject(h2_nsigma_pin);
        for ipin in range(0,len(arr_pin)-1):
            pin1 = arr_pin[ipin];
            pin2 = arr_pin[ipin+1];
            h1 = slice_histogram(h2_nsigma_pin,pin1,pin2,"Y",False);
            h1.SetName("h1nsigma_{0}_pin{1}".format(parnames[ip],ipin));
            h1.SetTitle("nsigma TPC {0} , {1} < p_{{in}} < {2} GeV/c".format(parnames[ip],pin1,pin2));
            h1.SetXTitle("n#sigma^{{TPC}}_{{{0}}}".format(parsymbols[ip]));
            h1.SetYTitle("Number of Tracks");

            f1 = TF1("f1_{0}_pin{1}".format(parnames[ip],ipin),"gaus",-10,+10);
            f1.SetNpx(2000);
            f1.SetParameter(1,0);
            f1.SetParameter(2,1);
            h1.Fit(f1,"SMEL","",-5,+5);

            outfile.WriteTObject(h1);
            outfile.WriteTObject(f1);
    outfile.Close();
def check_nsigma_tof(list_hierarchies, arr_pin, arr_eta, suffix=""):
    print(sys._getframe().f_code.co_name);
    print(list_hierarchies);
    fm = FileManager(list_hierarchies);
    parnames = ["Ele","Pi","Ka","Pr"];
    parsymbols = ["e","#pi","K","p"];
    np = len(parnames);

    outname = "pid_nsigma_TOF{0}.root".format(suffix);
    outfile = TFile.Open(outname,"RECREATE");

    for ip in range(0,np):
        h2_nsigma_pin = fm.get("TOFnSig{0}_pIN".format(parnames[ip]));
        h2_nsigma_pin.Sumw2();
        outfile.WriteTObject(h2_nsigma_pin);
        for ipin in range(0,len(arr_pin)-1):
            pin1 = arr_pin[ipin];
            pin2 = arr_pin[ipin+1];
            h1 = slice_histogram(h2_nsigma_pin,pin1,pin2,"Y",False);
            h1.SetName("h1nsigma_{0}_pin{1}".format(parnames[ip],ipin));
            h1.SetTitle("nsigma TOF {0} , {1} < p_{{in}} < {2} GeV/c".format(parnames[ip],pin1,pin2));
            h1.SetXTitle("n#sigma^{{TOF}}_{{{0}}}".format(parsymbols[ip]));
            h1.SetYTitle("Number of Tracks");

            f1 = TF1("f1_{0}_pin{1}".format(parnames[ip],ipin),"gaus",-10,+10);
            f1.SetNpx(2000);
            f1.SetParameter(1,0);
            f1.SetParameter(2,1);
            h1.Fit(f1,"SMEL","",-5,+5);

            outfile.WriteTObject(h1);
            outfile.WriteTObject(f1);
    outfile.Close();

if __name__ == "__main__":
    list_hierarchies = [
        "../AnalysisResults.root"
        ,"analysis-track-selection"
	,"output"
	,"TrackBarrel_BeforeCuts"
    ];

    #arr_pin = np.array([0.01,0.02,0.04,0.06,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,4,6,8,10],dtype=float);
    arr_pin = np.array([0.05,0.1,0.5,1.0,10],dtype=float);
    arr_eta = np.array([-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],dtype=float);

    check_nsigma_tpc(list_hierarchies, arr_pin, arr_eta,"");
    check_nsigma_tof(list_hierarchies, arr_pin, arr_eta,"");

