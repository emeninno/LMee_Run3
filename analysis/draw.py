import numpy as np
import datetime
import ROOT
from  file_manager import FileManager
from ROOT import TFile, TDirectory, THashList, TH1F, TH2F, TCanvas, TLegend, TPaveText, TPython, TMath
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#_____________________________________________________________________
def draw_dedx_pin(cutnames):

    list_hierarchies = [
        "AnalysisResults_table_maker.root"
        ,"table-maker"
        ,"output"
    ];
    fm = FileManager(list_hierarchies);

    list_ev = fm.get("Event_AfterCuts");
    h1vtx = list_ev.FindObject("VtxZ");
    nev = h1vtx.GetEntries();
    print("nev = {0} M events".format(nev/1e+6));
    gStyle.SetPalette(55);
    ncuts = len(cutnames);

    date = datetime.date.today().strftime("%Y%m%d");

    for i in range(0,ncuts):
        c1 = TCanvas("c{0}".format(i),"cv_{0}".format(cutnames[i]),0,0,800,800);
        c1.SetLogx(1);
        c1.SetLogz(1);
        c1.SetTicks(1,1);
        c1.SetMargin(0.12,0.12,0.12,0.03);

        frame = c1.DrawFrame(0.01,0,10,200);
        list_track_barrel = fm.get("TrackBarrel_{0}".format(cutnames[i]));
        h2 = list_track_barrel.FindObject("TPCdedx_pIN");
        h2.SetDirectory(0);
        h2.SetContour(100);
        h2.GetXaxis().SetRangeUser(0.01,10);
        h2.GetYaxis().SetRangeUser(0,200);
        h2.GetXaxis().SetTitleOffset(1.5);
        h2.GetYaxis().SetTitleOffset(1.5);
        h2.GetXaxis().SetTitleSize(0.035);
        h2.GetYaxis().SetTitleSize(0.035);
        h2.GetXaxis().SetLabelSize(0.035);
        h2.GetYaxis().SetLabelSize(0.035);
        h2.GetXaxis().SetMoreLogLabels(1);
        h2.Draw("colz");
        ROOT.SetOwnership(h2,False);

        c1.Modified();
        c1.Update();
        c1.SaveAs("{0}_pilotbeam_pp_900GeV_TPCdedx_{1}.eps".format(date,cutnames[i]));
        c1.SaveAs("{0}_pilotbeam_pp_900GeV_TPCdedx_{1}.pdf".format(date,cutnames[i]));
        ROOT.SetOwnership(c1,False);

    fm.close();
#_____________________________________________________________________
def draw_ap():
    list_hierarchies = [
        "AnalysisResults_apass4_noQC_V0.root"
        ,"v0-selector"
    ];
    fm = FileManager(list_hierarchies);

    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,800,800);
    #c1.SetLogz(1);
    c1.SetTicks(1,1);
    c1.SetMargin(0.12,0.12,0.15,0.03);

    hV0AP = fm.get("hV0APplot");
    ROOT.SetOwnership(hV0AP,False);
    hV0AP.SetXTitle("#alpha = #frac{p_{L}^{+} - p_{L}^{-}}{p_{L}^{+} + p_{L}^{-}}");
    hV0AP.SetYTitle("q_{T} (GeV/c)");
    hV0AP.Draw("colz");
    hV0AP.GetXaxis().SetTitleOffset(1.5);
    hV0AP.GetYaxis().SetTitleOffset(1.5);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0AP.eps".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0AP.pdf".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0AP.png".format(date));
    ROOT.SetOwnership(c1,False);
#_____________________________________________________________________
def draw_v0_dedx():
    list_hierarchies = [
        "AnalysisResults_apass4_noQC_V0.root"
        ,"track-pid-qa"
    ];
    fm = FileManager(list_hierarchies);

    parnames = ["all","El","Pi","Ka","Pr"];
    np = len(parnames);
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,1200,800);
    c1.Divide(3,2);

    for ip in range(0,np):
        p1 = c1.cd(ip+1);
        p1.SetMargin(0.12,0.12,0.15,0.03);
        p1.SetLogx(1);
        p1.SetLogz(1);
        h2_dedx_pin = fm.get("h2TPCdEdx_Pin_{0}".format(parnames[ip]));
        h2_dedx_pin.Draw("colz");
        ROOT.SetOwnership(h2_dedx_pin,False);
        h2_dedx_pin.SetXTitle("p_{in} (GeV/c)");
        h2_dedx_pin.SetYTitle("TPC dE/dx (a.u.)");
        h2_dedx_pin.GetXaxis().SetTitleOffset(1.5);
        h2_dedx_pin.GetYaxis().SetTitleOffset(1.5);
        h2_dedx_pin.GetXaxis().SetTitleSize(0.04);
        h2_dedx_pin.GetYaxis().SetTitleSize(0.04);
        h2_dedx_pin.Draw("colz");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_TPCdEdx.eps".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_TPCdEdx.pdf".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_TPCdEdx.png".format(date));
    ROOT.SetOwnership(c1,False);
#_____________________________________________________________________
def draw_v0_tofbeta():
    list_hierarchies = [
        "AnalysisResults_apass4_noQC_V0.root"
        ,"track-pid-qa"
    ];
    fm = FileManager(list_hierarchies);

    parnames = ["all","El","Pi","Ka","Pr"];
    np = len(parnames);
    gStyle.SetPalette(55);
    c1 = TCanvas("c0","c0",0,0,1200,800);
    c1.Divide(3,2);

    for ip in range(0,np):
        p1 = c1.cd(ip+1);
        p1.SetMargin(0.12,0.12,0.15,0.03);
        p1.SetLogx(1);
        p1.SetLogz(1);
        h2_beta_pin = fm.get("h2TOFbeta_Pin_{0}".format(parnames[ip]));
        h2_beta_pin.Draw("colz");
        ROOT.SetOwnership(h2_beta_pin,False);
        h2_beta_pin.SetXTitle("p_{in} (GeV/c)");
        h2_beta_pin.SetYTitle("TOF #beta");
        h2_beta_pin.GetXaxis().SetTitleOffset(1.5);
        h2_beta_pin.GetYaxis().SetTitleOffset(1.5);
        h2_beta_pin.GetXaxis().SetTitleSize(0.04);
        h2_beta_pin.GetYaxis().SetTitleSize(0.04);
        h2_beta_pin.Draw("colz");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_TOFbeta.eps".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_TOFbeta.pdf".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_TOFbeta.png".format(date));
    ROOT.SetOwnership(c1,False);
#_____________________________________________________________________
def draw_v0_mass():
    list_hierarchies = [
        "AnalysisResults_apass4_noQC_V0.root"
        ,"v0-selector"
    ];
    fm = FileManager(list_hierarchies);

    parnames_symbol = ["#gamma","K_{S}^{0}","#Lambda","#bar{#Lambda}"];
    parnames = ["Gamma","K0S","Lambda","AntiLambda"];
    titles = ["e^{+}e^{-}","#pi^{+}#pi^{-}","p#pi^{-}","#bar{p}#pi^{+}"];
    np = len(parnames);
    gStyle.SetPalette(55);
    c1 = TCanvas("c1","c1",0,0,1200,800);
    c1.Divide(2,2);

    for ip in range(0,np):
        p1 = c1.cd(ip+1);
        p1.SetMargin(0.12,0.03,0.12,0.03);
        p1.SetTicks(1,1);
        ROOT.SetOwnership(p1,False);
        h1 = fm.get("hMass{0}".format(parnames[ip]));
        h1.SetDirectory(0);
        h1.SetXTitle("M_{{{0}}} (GeV/c^{{2}})".format(titles[ip]));
        h1.SetYTitle("Number of Entries");
        h1.GetXaxis().SetTitleOffset(1.1);
        h1.GetYaxis().SetTitleOffset(1.1);
        h1.GetXaxis().SetTitleSize(0.05);
        h1.GetYaxis().SetTitleSize(0.05);
        h1.GetXaxis().SetLabelSize(0.04);
        h1.GetYaxis().SetLabelSize(0.04);
        ROOT.SetOwnership(h1,False);
        h1.Sumw2();
        h1.Draw("HE0");

        txt = TPaveText(0.15,0.75,0.4,0.92,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(12);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.05);
        txt.AddText("pp at #sqrt{s} = 900 GeV");
        txt.AddText("pilot beam data");
        txt.AddText("{0} #rightarrow {1}".format(parnames_symbol[ip],titles[ip]));
        txt.Draw();
        ROOT.SetOwnership(txt,False);


    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_mass.eps".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_mass.pdf".format(date));
    c1.SaveAs("{0}_pilotbeam_pp_900GeV_V0_mass.png".format(date));
    ROOT.SetOwnership(c1,False);

#_____________________________________________________________________
#_____________________________________________________________________
if __name__ == "__main__":
    #cutnames = ["BeforeCuts","jpsiO2MCdebugCuts"];
    #draw_dedx_pin(cutnames);
    draw_ap();
    #draw_v0_dedx();
    #draw_v0_tofbeta();
    #draw_v0_mass();

