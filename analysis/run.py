import os, sys, shutil
import math
import argparse
import numpy as np
import ctypes
import yaml
import ROOT
from ROOT import TFile
from analyze_pair import analyze_mee_ptee, analyze_mee_ptee_dcaee, analyze_mee_ptee_efficiency

parser = argparse.ArgumentParser('Example program');
parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
args = parser.parse_args();

filename = args.input;
with open(args.config, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
#_________________________________________________________________________________________
def run(filename,config,ismc):
    print(sys._getframe().f_code.co_name);
    arr_mee = np.array(config["common"]["mee_bin"],dtype=float);
    arr_ptee = np.array(config["common"]["ptee_bin"],dtype=float);
    arr_dcaee = np.array(config["common"]["dcaee_bin"],dtype=float);
    print("mass binning = ",arr_mee);
    print("pTee binning = ",arr_ptee);
    print("DCAee binning",arr_dcaee);
    print("ismc = ",ismc);
    print("input = ",filename);
    rootfile = TFile.Open(filename,"READ");

    cutnames = config[args.type]["cutnames"]
    print(cutnames); 
    nc = len(cutnames);

    if config["common"]["do_mee_ptee_dcaee"] == True:
        outname = "output_{0}_mee_ptee_dcaee_{1}_{2}TeV_{3}.root".format(args.type,config["common"]["system"],config["common"]["energy"],config["common"]["period"]);
        print("output file name = ",outname);
        outfile = TFile(outname,"RECREATE");

        if ismc:
            for ic in range(0,nc):
                outlist = analyze_mee_ptee_efficiency(rootfile,cutnames[ic],arr_mee,arr_ptee);
                outlist.SetOwner(True);
                outfile.WriteTObject(outlist);
                outlist.Clear();
        else:
            for ic in range(0,nc):
                outlist = analyze_mee_ptee_dcaee(rootfile,cutnames[ic],arr_mee,arr_ptee,arr_dcaee);
                outlist.SetOwner(True);
                outfile.WriteTObject(outlist);
                outlist.Clear();
    elif config["common"]["do_mee_ptee"] == True:
        outname = "output_{0}_mee_ptee_{1}_{2}TeV_{3}.root".format(args.type,config["common"]["system"],config["common"]["energy"],config["common"]["period"]);
        print("output file name = ",outname);
        outfile = TFile(outname,"RECREATE");

        if ismc:
            for ic in range(0,nc):
                outlist = analyze_mee_ptee_efficiency(rootfile,cutnames[ic],arr_mee,arr_ptee);
                outlist.SetOwner(True);
                outfile.WriteTObject(outlist);
                outlist.Clear();
        else:
            for ic in range(0,nc):
                outlist = analyze_mee_ptee(rootfile,cutnames[ic],arr_mee,arr_ptee,arr_dcaee);
                outlist.SetOwner(True);
                outfile.WriteTObject(outlist);
                outlist.Clear();
    else:
        print("please check what to do in",args.config);

    outfile.Close();
    rootfile.Close();
#_________________________________________________________________________________________
ismc = False;
if args.type == "data":
    ismc = False;
elif args.type == "mc":
    ismc = True;
else:
    print("unknown type.sys.exit()");
    sys.exit();
run(filename,config,ismc);
