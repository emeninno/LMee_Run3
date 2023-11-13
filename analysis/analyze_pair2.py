import os, sys, shutil
import math
import numpy as np
import ctypes
import ROOT
from ROOT import TFile, TDirectory, TList, THashList, TH1D, TH2D, TH3D, TH1F, TH2F
from ROOT import gROOT, gSystem, gStyle, gPad
from ROOT import TCanvas
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange
#from file_manager import FileManager
from histo_manager import slice_histogram, rebin_histogram, get_R_factor, get_R_factor_2D, get_bkg2D, get_bkg_subtracted, get_SBratio

class PairAnalysis:
# pair analysis class
	def __init__(self):
		print("default constructor is called");
	def __init__(self,filename,cutname):
		
		outlist = TList();
		outlist.SetName("output histograms, cut {0}".format(cutname));

		print("file = {0}, cuts = {1}".format(filename, cutname));
		rootfile = TFile.Open(filename,"READ");
		rootfile.ls();
		self.rootdir = rootfile.Get("analysis-event-mixing");
		self.rootdir.ls();
		self.rootdir.cd();
		self.list_mix = self.rootdir.Get("output");
		self.list_mix.ls();
		#lists for pos.neg.+neg.pos, pos.pos., and neg.neg. pairs 
		parnames = ["PM","PP","MM"];
		parsymbols = ["ULS","LSpp","LSmm"];
		hMass = {'hMass_0': TH1F(), 'hMass_1': TH1F, 'hMass_2': TH1F}; #need to be improved
		hPt = {'hPt_0': TH1F(), 'hPt_1': TH1F, 'hPt_2':TH1F};
		hMass_Pt = {'hMass_Pt_0': TH2F(), 'hMass_Pt_1': TH2F, 'hMass_Pt_2':TH2F};
		np = len(parnames);
		for ip in range(0,np):
			self.list2 = self.list_mix.FindObject("PairsBarrelME{0}_lmee_{1}".format(parnames[ip],cutname));
			self.list2.Print();
			hMass[ip] = self.list2.FindObject("Mass");
			hPt[ip] = self.list2.FindObject("Pt");
			hMass[ip].SetName("hMass_mix_{0}". format(parnames[ip]));
			hMass[ip].SetTitle("hMass_mix_{0}". format(parnames[ip]));
			hPt[ip].SetName("hPt_mix_{0}". format(parnames[ip]));
			hPt[ip].SetTitle("hPt_mix_{0}". format(parnames[ip]));
			hMass_Pt[ip] = self.list2.FindObject("Mass_Pt");
			hMass_Pt[ip].SetName("hMassVsPt_mix_{0}". format(parnames[ip]));
			hMass_Pt[ip].SetTitle("hMassVsPt_mix_{0}". format(parnames[ip]));
		Rfac = get_R_factor(hMass[0],None,hMass[1], hMass[2]);
		print("Rfactor (1D) histo: {}".format(Rfac));

		Rfac2D = get_R_factor_2D(hMass_Pt[0],None,hMass_Pt[1], hMass_Pt[2]);
		print("Rfactor (2D) histo: {0}".format(Rfac2D));
		
		outlist.Add(Rfac);
		outlist.Add(Rfac2D);

		#get signal and background pairs from same-event pairs
		self.rootdir_se = rootfile.Get("analysis-same-event-pairing");
		self.rootdir_se.ls();
		self.rootdir_se.cd();
		self.list_se = self.rootdir_se.Get("output");
		self.list_se.ls();
		hMassSE = {'hMass_PM': TH1F(), 'hMass_PP': TH1F, 'hMass_MM':TH1F};
		hPtSE = {'hPt_PM': TH1F(), 'hPt_PP': TH1F, 'hPt_MM':TH1F};
		hMass_PtSE = {'hMass_Pt_PM': TH2F(), 'hMass_Pt_PP': TH2F, 'hMass_Pt_MM':TH2F};

		for ip in range(0,np):
			self.list_se2 = self.list_se.FindObject("PairsBarrelSE{0}_lmee_{1}".format(parnames[ip],cutname));
			hMassSE[ip] = self.list_se2.FindObject("Mass");
			hPtSE[ip] = self.list_se2.FindObject("Pt");
			hMassSE[ip].SetName("hMass_se_{0}". format(parnames[ip]));
			hMassSE[ip].SetTitle("hMass_se_{0}". format(parnames[ip]));									                        
			hPtSE[ip].SetName("hPt_se_{0}". format(parnames[ip]));
			hPtSE[ip].SetTitle("hPt_se_{0}". format(parnames[ip]));
			hMass_PtSE[ip] = self.list_se2.FindObject("Mass_Pt");
			hMass_PtSE[ip].SetName("hMassVsPt_se_{0}". format(parnames[ip]));														
			hMass_PtSE[ip].SetTitle("hMassVsPt_se_{0}". format(parnames[ip]));

		#get 2D background from LSpp and LSmm pairs
		hBkg = get_bkg2D(hMass_PtSE[1], hMass_PtSE[2]);
		outlist.Add(hBkg);	
		
		#get 2D corrected background (corrected Rfactor) and spectrum
		hSpectrum = hMass_PtSE[0].Clone("hSpectrum");
		hBkg_corr = hBkg.Clone("hBkg_corr");
		
		hBkg_corr.Multiply(Rfac2D);
		hSpectrum.Add(hBkg_corr,-1);
		
		outlist.Add(hMassSE[0]);
		outlist.Add(hMassSE[1]);
		outlist.Add(hMassSE[2]);
		outlist.Add(hPtSE[0]);
		outlist.Add(hPtSE[1]);
		outlist.Add(hPtSE[2]);
		outlist.Add(hMass_PtSE[0]);
		outlist.Add(hMass_PtSE[1]);
		outlist.Add(hMass_PtSE[2]);
		outlist.Add(hBkg_corr);
		outlist.Add(hSpectrum);
		#drawing
		c1 = TCanvas("c1","2D mass,pT spectrum",200,10,700,500) #make nice Canvas
		c1.SetGrid();	
		hSpectrum.Draw();		

	def __new__(self):
		return outlist;

#	def analyze_mee_ptee(rootfile,cutname,arr_mee, arr_ptee):
 #  	 print(sys._getframe().f_code.co_name);
 #  	 outlist = TList();
 #   	outlist.SetName(cutname);
 #   	return outlist;



#_________________________________________________________________________________________
#to be implemented
#	def analyze_mee_ptee_dcaee(rootfile,cutname,arr_mee, arr_ptee, arr_dcaee):
 #   	print(sys._getframe().f_code.co_name);
  #  	outlist = TList();
   # 	outlist.SetName(cutname);
    #	return outlist;
#_________________________________________________________________________________________
#to be implemented!
#	def analyze_mee_ptee_efficiency(rootfile,cutname,arr_mee, arr_ptee):
 #   	print(sys._getframe().f_code.co_name);
  #  	outlist = TList();
   # 	outlist.SetName(cutname);
    #	return outlist;
#____#_____________________________________________________________________________________

	    
#need to define the mee, ptee and dcaee bins
if __name__ == "__main__":
	ana = PairAnalysis("../AnalysisResults.root","eNSigmaRun3_strongNSigE");
	print(ana);

arr_mee = np.array([0.05,0.1,0.5,1.0,10],dtype=float);
arr_ptee = np.array([-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],dtype=float);
arr_dcaee = np.array([-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],dtype=float);

# request user action before ending (and deleting graphic window)
#raw_input ( ' Press <ret > to end -> ' ) 
