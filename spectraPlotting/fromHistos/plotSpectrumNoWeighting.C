#include "../../generalFunctions/random_functions.C"
#include "../../generalFunctions/cuts.C"
#include "../../generalFunctions/plotting.C"
#include "../../generalFunctions/diElecFuncs.C"
#include "../../generalFunctions/readInOut.C"
#include "../../generalFunctions/runArrays.C"
#include "../../generalFunctions/texLabels.C"
#include "../../paths/path_locations.C"

void setMultClass(Int_t whichMultClass, TString& phiVnoCutsHistName, TString& phiVCutsHistName, Float_t* zAxisRange);
std::vector<Int_t> getBinRanges(Float_t xRange[], Float_t yRange[], Float_t zRange[], const TH3F* inputSpectrum);
TH1F* getRfactorProjections(const std::vector<TH2F*> vec_inputs, const std::vector<Int_t> binRanges, Int_t axis = 0);
TH1F* getEffProjections(const std::vector<TH2F*> vec_inputs, const std::vector<Int_t> binRanges, Int_t axis = 0);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Takes spectrum and efficiency correction inputs from the dielectron      //
//  framework. Produces the following plots in 2D, and projected onto the    //
//  mass and pairPt axes:                                                    //
//    - unlike-sign and combinatorial background                             //
//    - statistical significance and signal-to-background ratio              //
//    - acceptance correction factor (R factor)                              //
//    - spectra with and with correction and R factors                       //
//    - efficiencies from all sources (including phiV)                       //
//                                                                           //
//  Function arguments:                                                      //
//    - whichCutSet:    "kCutSet1", etc                                      //
//                        - specified in, e.g., Config_acapon.C              //
//    - filesPath:      output from dielectron framework.                    //
//                        - e.g. "cent/16q/2018_08_27_kCutSet1_00100"        //
//    - effCorrPath:    e.g "cent/18f3/2018_08_27_kCutSet1_00100"            //
//    - saveName:       Specify the name of output ROOT file where the       //
//                      results will be saved to. The full path to the       //
//                      correct directory, based on <filesPath>, will be     //
//                      prepended.                                           //
//                        - If empty, plots will not be saved.               //
//                                                                           //
//    - whichMultClass: Specify the desired multiplicity class. Script       //
//                       assumes data input histograms are mass/pairPt/mult. //
//                         - 0 = 00100, 1 = 0020, 2 = 2040,                  //
//                           3 = 4060,  4 = 0100, 5 = 00100(course mult bins)//
//                       mass                                                //
//                                                                           //
// Analysis specific customisations can be found in 'script settings' below. //
///////////////////////////////////////////////////////////////////////////////


void plotSpectraNoWeighting(TString filesPath,
                            TString effCorrPath  = "",
                            Int_t whichMultClass = 0,
                            TString saveLocation = "")
{

  //------------------------- Script settings ---------------------------//
  // Write "ALICE Preliminary" on all plots
  Bool_t prelimPlots = kFALSE;

  // - Dielectron output folder, cut setting,  and 3D histogram name
  std::vector<TString> diElecOutputNames  = {"acapon_out_0", "kCutSet1", "pInvMass_PairPt_CentralityV0A",};
  // - Dielectron efficiency output folder, cut setting
  std::vector<TString> diElecEffNames = {"efficiency0", "kCutSet1"};
  // Which source of electrons to be used
  Int_t whichSource = 0; ///0=all, 1=LF, 2=HF

  // - - - - PhiV efficiency settings
  // Pair cuts are not applied to the dielectron efficiency task
  Bool_t applyPhiVcorrection = kTRUE;
  TString phiVfileName       = "cent/18f3/multDepPhiVcutEfficiency";
  /* TString phiVfileName       = "cent/18f3/multDepPhiVcutEfficiency_noITSsharedCut"; */

  // - - - - Efficiency settings
  // Rebin efficiency. Efficiency bins do not need to match the data, so a
  // rebinning scheme can be used to suppress statistical fluctuations.
  Bool_t rebinEff = kTRUE;
  std::vector<Double_t> newMassBinsEff   = plotting::vec_massBins;
  std::vector<Double_t> newPairPtBinsEff = plotting::vec_pairPtBins;
  Int_t numEffRebin[2] = {2,10};
  // After rebinning, one can apply bin filling scheme with following numbers as
  // minimum allowed tracks per given bin. If not satisfied, an average over the
  // neighbouring bins is used.
  // The flag will also enforce that the efficiency is less than 1.
  Bool_t fillEmptyEffBins = kFALSE;
  Int_t minGenPairsInBin  = 35;
  Int_t minAccPairsInBin  = 25;

  // Point at which the R factor is set to 1 (mass axis)
  Float_t RfacOne = 0.3;
  Bool_t scaleByPhysSel  = kTRUE;
  Bool_t scaleByBinWidth = kTRUE;

  // Restrict mass and pairPt for projections
  Float_t xAxisRange[2] = {0, 5};
  Float_t yAxisRange[2] = {0, 8};
  // Set with wichMultClass argument
  Float_t zAxisRange[2] = {};

  //TODO: assign correct axis labels based on this choice
  // Pick two axes to be plotted
  TString projectionAxes = "yxe";

  // Option to do a standard rebin on the x-y axes
  // Will be applied to data histograms before calculation begins
  // If flag set and vectors empty, standard rebin will be used
  Bool_t rebin2D = kTRUE;
  std::vector<Float_t> newMassBins   = plotting::vec_massBins;
  std::vector<Float_t> newPairPtBins = plotting::vec_pairPtBins;
  Int_t numRebin[2] = {1,4}; // Standard rebin values

  // - - - - Plotting settings
  // Efficiency plotting options
  Bool_t restrictEffPlots = kFALSE;
  Float_t minPairPtEffPlots = 0;
  Float_t maxPairPtEffPlots = 6;
  Float_t minMassEffPlots   = 0;
  Float_t maxMassEffPlots   = 5;
  TString effPlotStyle = "COLZ";
  // Set if only using LF correction
  Bool_t zoomXaxis = kFALSE;

  // 1D projections will be saved from the two chosen axes under the following
  // names. Axis name appended (currently assuming mass and pairPt)
  std::vector<TString> saveName = {"ULS", "combBackgr", "rawSignal", "rawSigNoR",
                                   "correctedSignal", "SB", "signficance", "Rfac",
                                  "diElecEff", "phiVeff", "effectiveEff"};
  //---------------------------------------------------------------------//

  // Use much wider binning if looking mult classes
  if(whichMultClass != 0){
    std::cout << "Setting multiplcitiy bins!" << std::endl;
    newMassBins      = plotting::vec_massBinsMultDep;
    newPairPtBins    = plotting::vec_pairPtBinsMultDep;
    newMassBinsEff   = plotting::vec_massBinsMultDep;
    newPairPtBinsEff = plotting::vec_pairPtBinsMultDep;
  }

  TStopwatch* watch = new TStopwatch();
  watch->Start();

  // Function returns following variables with correct multiplicity labels appended
  // and correct multiplicity range set
  TString phiVnoCutsHistName = "noCut";
  TString phiVcutsHistName   = "cut";
  setMultClass(whichMultClass, phiVnoCutsHistName, phiVcutsHistName, zAxisRange);

  filesPath.Prepend(paths::results_loc + paths::diElec_dir + "data/");
  filesPath.Append(".root");

  // Number of events after cuts stored here
  Double_t numEvents = 0.0;

  // Get input spectra (same event and mixed event(me) )
  // Vector contents: {ULS, PosLS, NegLS, ULSme, PosLSme, NegLSme}
  Bool_t hasMixedEvents = kTRUE;
  std::vector<TH3F*> vec_inputSpectra = readInOut::getInputSpectra(filesPath, diElecOutputNames, numEvents);
  if(vec_inputSpectra.size() == 0){
    return;
  }
  // Sometimes mixing doesn't work because.....
  // If so, set flag
  if(vec_inputSpectra[3]->GetEntries() == 0){
    printError("No mixed event plots found!");
    hasMixedEvents = kFALSE;
  }

  // Retrieve and store min/max bins for x,y,z
  std::cout << "Original bin ranges" << std::endl;
  std::vector<Int_t> binRanges = getBinRanges(xAxisRange, yAxisRange, zAxisRange, vec_inputSpectra[0]);

  // Restrict Z axis before projection and store projections
  std::vector<TH2F*> vec_inputSpectra2D;
  for(Int_t i = 0; i < vec_inputSpectra.size(); ++i){
    vec_inputSpectra[i]->GetZaxis()->SetRange(binRanges[4], binRanges[5]);
    TH2F* histProj = (TH2F*)vec_inputSpectra[i]->Project3D(projectionAxes);
    vec_inputSpectra2D.push_back(histProj);
  }

  if(rebin2D){
    if(newMassBins.size() == 0){
      std::cout << "Performing simple 2D rebin on all input spectra" << std::endl;
      std::cout << "X axis rebin: " << numRebin[0] << ", y-axis rebin: " << numRebin[1] << std::endl;
      for(Int_t i = 0; i < vec_inputSpectra2D.size(); ++i){
        vec_inputSpectra2D[i]->Rebin2D(numRebin[0], numRebin[1]);
      }
    }
    else if(newMassBins.size() > 1){
      for(Int_t i = 0; i < vec_inputSpectra2D.size(); ++i){
        vec_inputSpectra2D[i] = plotting::rebin2Dhist(vec_inputSpectra2D[i], newMassBins, newPairPtBins);
        if(vec_inputSpectra2D[i] == 0x0){
          return;
        }
      }
    }
    // If rebinning, find new bins for projections
    std::cout << "Bin ranges after rebinning" << std::endl;
    binRanges = getBinRanges(xAxisRange, yAxisRange, zAxisRange, vec_inputSpectra[0]);
  }

  //----------- Get generated smeared and accepted histograms
  //Vector contents: {generated, accepted, diElecEfficiency}
  std::vector<TH2F*> vec_effInputs;
  if(effCorrPath != ""){
    effCorrPath.Prepend(paths::results_loc + paths::diElec_dir + "effCorrs/");
    effCorrPath.Append(".root");

    vec_effInputs = readInOut::getEfficiencySpectra(effCorrPath, diElecEffNames, whichSource);
    if(vec_effInputs.size() < 2){
      return;
    }

    if(rebinEff){
      if(newMassBinsEff.size() == 0){
        std::cout << "Standard rebinning applied to efficiency histograms" << std::endl;
        vec_effInputs[0]->Rebin2D(numEffRebin[0], numEffRebin[1]);
        vec_effInputs[1]->Rebin2D(numEffRebin[0], numEffRebin[1]);
      }
      else if(newMassBinsEff.size() > 1){
        for(Int_t i = 0; i < vec_effInputs.size(); ++i){
          vec_effInputs[i] = plotting::rebin2Dhist(vec_effInputs[i], newMassBinsEff, newPairPtBinsEff);
          if(vec_effInputs[i] == 0x0){
            return;
          }
        }
      }
    }
    if(fillEmptyEffBins){
      std::cout << "Filling empty bins" << std::endl;
      fixEfficiency(vec_effInputs[0], vec_effInputs[1]);
    }

    TH2F* histEffDielec2D = (TH2F*)vec_effInputs[0]->Clone("histEffDielec2D");
    histEffDielec2D->Reset();
    histEffDielec2D->Divide(vec_effInputs[1], vec_effInputs[0], 1, 1, "B");
    vec_effInputs.push_back(histEffDielec2D);

  }else{ // Push NULL pointer so plotting functions work
    vec_effInputs.push_back(0x0);
  }

  //##########################################################
  // Get phiV cut efficiency
  TH2F* phiVeff = 0x0;
  // Projection plots
  TH1F* phiVmassEff   = 0x0;
  TH1F* phiVpairPtEff = 0x0;
  if(applyPhiVcorrection && vec_effInputs.size() != 0){
    TFile* inPhiVfile = TFile::Open(paths::results_loc + paths::conversionRej_dir+ phiVfileName + ".root", "READ");
    if(!inPhiVfile){
      printError("PhiV eff correction file not found.");
      return;
    }
    TH2F* phiVnoCut = (TH2F*)inPhiVfile->Get(phiVnoCutsHistName);
    TH2F* phiVcut   = (TH2F*)inPhiVfile->Get(phiVcutsHistName);
    if(!phiVnoCut | !phiVcut){
      printError("PhiV input plot not found");
      return;
    }

    // Only three bins in mass before 100 MeV.
    if(rebin2D && newMassBins.size () == 0){
      phiVcut->Rebin2D(numRebin[0] > 3 ? 3 : numRebin[0], numRebin[1]);
      phiVnoCut->Rebin2D(numRebin[0] > 3 ? 3 : numRebin[0], numRebin[1]);
    }
    phiVeff = (TH2F*)phiVcut->Clone("phiVeff");
    phiVeff->Reset();
    phiVeff->Divide(phiVcut, phiVnoCut, 1, 1, "B");

    // PhiV projections
    TH1F* phiVmassCutTemp   = (TH1F*)phiVcut->ProjectionX("phVmassNumTemp",     binRanges[2], binRanges[3], "e");
    TH1F* phiVmassNoCutTemp = (TH1F*)phiVnoCut->ProjectionX("phVmassNoCutTemp", binRanges[2], binRanges[3], "e");

    TH1F* phiVpairPtCutTemp   = (TH1F*)phiVcut->ProjectionY("phVpairPtNumTemp",     binRanges[0], binRanges[1], "e");
    TH1F* phiVpairPtNoCutTemp = (TH1F*)phiVnoCut->ProjectionY("phVpairPtNoCutTemp", binRanges[0], binRanges[1], "e");

    phiVmassEff = (TH1F*)phiVmassCutTemp->Clone("phiVmassEff");
    phiVmassEff->Clear();
    phiVmassEff->Divide(phiVmassCutTemp, phiVmassNoCutTemp, 1, 1, "e");

    phiVpairPtEff = (TH1F*)phiVpairPtCutTemp->Clone("phiVpairPtEff");
    phiVpairPtEff->Reset();
    phiVpairPtEff->Divide(phiVpairPtCutTemp, phiVpairPtNoCutTemp, 1, 1, "B");
    // If looking at mass ranges above 100 GeV, set phiV on pairPt axis to 1
    // without error
    if(xAxisRange[0] > 0.1){
      for(Int_t i = 0; i < phiVpairPtEff->GetNbinsX(); ++i){
        phiVpairPtEff->SetBinContent(i, 1);
        phiVpairPtEff->SetBinError(i, 0);
      }
    }

    // Set all bins above 100 MeV to 1 with 0 error
    Int_t lastPhiVbin = phiVeff->GetXaxis()->FindBin(0.1000000000001);
    Int_t lastPtBin = phiVeff->GetNbinsY();
    for(Int_t i = lastPhiVbin; i <= phiVeff->GetNbinsX(); i++){

      phiVmassEff->SetBinContent(i, 1);
      phiVmassEff->SetBinError(i, 0.00000000000001);

      for(Int_t j = 0; j <= lastPtBin; j++){
        phiVeff->SetBinContent(i, j, 1);
        phiVeff->SetBinError(i, j, 0.000000000000001);
      }
    }

    delete phiVmassCutTemp;
    delete phiVmassNoCutTemp;
    delete phiVpairPtCutTemp;
    delete phiVpairPtNoCutTemp;
    delete inPhiVfile;
  }

  //##############################################################
  // Begin creating final plots
  TH2F* R              = 0x0;
  TH2F* RforComp       = 0x0;
  if(hasMixedEvents){
    R        = (TH2F*)calcDiElecRfactor(vec_inputSpectra2D[3], vec_inputSpectra2D[4], vec_inputSpectra2D[5], RfacOne);
    RforComp = (TH2F*)calcDiElecRfactor(vec_inputSpectra2D[3], vec_inputSpectra2D[4], vec_inputSpectra2D[5], 5);
  }

  TH2F* combBackgr     = (TH2F*)calcDiElecBackgr(vec_inputSpectra2D[1], vec_inputSpectra2D[2]);
  TH2F* spectrum       = (TH2F*)calcDiElecSpectrum(vec_inputSpectra2D[0], combBackgr, R, vec_effInputs.back());
  if(applyPhiVcorrection && effCorrPath != ""){
    applyEffCorr(spectrum, phiVeff);
    vec_effInputs.push_back(phiVeff);
  }


  TH2F* rawSpectrumNoR = (TH2F*)calcDiElecSpectrum(vec_inputSpectra2D[0], combBackgr, 0x0, 0x0);
  TH2F* rawSpectrum    = (TH2F*)calcDiElecSpectrum(vec_inputSpectra2D[0], combBackgr, RforComp, 0x0);

  // Create vector containing final projections
  // Vector contents: {ULS, combBackgr, raw signal, raw signal no R factor, raw
  // signal, SB ratio, stat. significance, R factor}
  // Projections onto X axis
  std::vector<TH1F*> vec_xAxisProjections;
  vec_xAxisProjections.push_back((TH1F*)vec_inputSpectra2D[0]->ProjectionX("ULSxAxisProj", binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1F*)combBackgr->ProjectionX("combBackgrXaxisProj",     binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1F*)rawSpectrum->ProjectionX("rawSigXaxisProj",        binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1F*)rawSpectrumNoR->ProjectionX("rawSigNoRxAxisProj",  binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1F*)spectrum->ProjectionX("signalXaxisProj",           binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back(calcDiElecSB(vec_xAxisProjections[2], vec_xAxisProjections[1]));
  vec_xAxisProjections.push_back(calcDiElecSignificance(vec_xAxisProjections[2], vec_xAxisProjections[1]));
  if(hasMixedEvents){
    vec_xAxisProjections.push_back(getRfactorProjections(vec_inputSpectra2D, binRanges, 0)); // 7th object in vector
  }else{
    vec_xAxisProjections.push_back(0x0);
  }


  // Projections onto Y axis
  std::vector<TH1F*> vec_yAxisProjections;
  vec_yAxisProjections.push_back((TH1F*)vec_inputSpectra2D[0]->ProjectionY("ULSyAxisProj", binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1F*)combBackgr->ProjectionY("combBackgrYaxisProj",     binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1F*)rawSpectrum->ProjectionY("rawSigYaxisProj",        binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1F*)rawSpectrumNoR->ProjectionY("rawSigNoRyAxisProj",  binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1F*)spectrum->ProjectionY("signalYaxisProj",           binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back(calcDiElecSB(vec_yAxisProjections[2], vec_yAxisProjections[1]));
  vec_yAxisProjections.push_back(calcDiElecSignificance(vec_yAxisProjections[2], vec_yAxisProjections[1]));
  if(hasMixedEvents){
    vec_yAxisProjections.push_back(getRfactorProjections(vec_inputSpectra2D, binRanges, 1)); // 7th object in vector
  }
  else{
    vec_yAxisProjections.push_back(0x0);
  }

  //###########################################################3
  // Push efficiency projections into existing projection vectors
  // Contents pushed into vector (in final order):
  //          {diElec eff, phiVeff, effective efficiency}
  if(effCorrPath != ""){
    vec_xAxisProjections.push_back(getEffProjections(vec_effInputs, binRanges, 0));
    vec_yAxisProjections.push_back(getEffProjections(vec_effInputs, binRanges, 1));
    if(applyPhiVcorrection){
      vec_xAxisProjections.push_back(phiVmassEff);
      vec_yAxisProjections.push_back(phiVpairPtEff);
    }
    // Create "effective efficiency plots"
    // Raw signal/fully corrected signal
    TH1F* effEff_Xaxis = (TH1F*)vec_xAxisProjections[2]->Clone("effEffXaxis");
    if(!effEff_Xaxis->Divide(effEff_Xaxis, vec_xAxisProjections[4], 1, 1, "B")){
      printError("Xaxis effective efficiency division failed!");
    }
    vec_xAxisProjections.push_back(effEff_Xaxis);
    TH1F* effEff_Yaxis = (TH1F*)vec_yAxisProjections[2]->Clone("effEffYaxis");
    if(!effEff_Yaxis->Divide(effEff_Yaxis, vec_yAxisProjections[4], 1, 1, "B")){
      printError("Yaxis effective efficiency division failed!");
    }
    vec_yAxisProjections.push_back(effEff_Yaxis);
  }

  //###########################################################
  // ----- Event scaling and rebinning -------
  if(scaleByPhysSel && numEvents != 0){
    numEvents *= (zAxisRange[1] - zAxisRange[0])/100.;
    std::cout << "Scaling spectra by " << numEvents << " events" << std::endl;
    vec_xAxisProjections[4]->Scale(1./numEvents);
    vec_yAxisProjections[4]->Scale(1./numEvents);
    spectrum->Scale(1./numEvents);
  }

  // Scale projected plots by bin width
  if(scaleByBinWidth){
    for(Int_t i = 0; i < 5; ++i){
      vec_xAxisProjections[i]->Scale(1., "width");
      vec_yAxisProjections[i]->Scale(1., "width");
    }
  }

  //#########################################################
  //------ Create labels for plots
  TString multSel = TString::Format("%.0f - %.0f", zAxisRange[0], zAxisRange[1]);

  // Labels for plots.
  // Two sets:
  // 1) plots small/no y axis label ("left")
  // 2) plots with large Y axis labels ("shifted")
  Float_t yPositionTexLabel = prelimPlots ? 0.82 : 0.87;
  TLatex* texPrelimLeft = texLabels::getTitle(0.2, yPositionTexLabel, "ALICE Preliminary", 0.05);
  std::vector<TLatex*> massLabelsLeft   = texLabels::getMassTexLabels("MB", yAxisRange[0], yAxisRange[1], 0.2, yPositionTexLabel-0.05);
  std::vector<TLatex*> pairPtLabelsLeft = texLabels::getPairPtTexLabels("MB", xAxisRange[0], xAxisRange[1], 0.2, yPositionTexLabel-0.05);

  // Labels for the final spectra
  TLatex* texPrelimShifted = texLabels::getTitle(0.25, yPositionTexLabel, "ALICE Preliminary", 0.05);
  std::vector<TLatex*> massLabelsShifted   = texLabels::getMassTexLabels("MB", yAxisRange[0], yAxisRange[1], 0.25, yPositionTexLabel-0.05);
  std::vector<TLatex*> pairPtLabelsShifted = texLabels::getPairPtTexLabels("MB", xAxisRange[0], xAxisRange[1], 0.25, yPositionTexLabel-0.05);


  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  //#########################################################
  // Plot 2D raw inputs
  TCanvas* canvasRawInputs2D = new TCanvas("canvasRawInputs2D", "canvasRawInputs2D");
  canvasRawInputs2D->SetWindowSize(1700, 900);
  canvasRawInputs2D->SetCanvasSize(1650, 850);
  canvasRawInputs2D->Divide(3,1);
  canvasRawInputs2D->cd(1);
  gPad->SetLogz();
  plotting::Set2DTitles(vec_inputSpectra2D[0], "Unlike-sign Dielectron Pairs");
  vec_inputSpectra2D[0]->Draw("COLZ");
  canvasRawInputs2D->cd(2);
  gPad->SetLogz();
  plotting::Set2DTitles(vec_inputSpectra2D[1], "Positive like-sign dielectron pairs");
  vec_inputSpectra2D[1]->Draw("COLZ");
  canvasRawInputs2D->cd(3);
  gPad->SetLogz();
  plotting::Set2DTitles(vec_inputSpectra2D[2], "Negative like-sign dielectron pairs");
  vec_inputSpectra2D[2]->Draw("COLZ");


  // Define margins for projections plot
  Float_t marginRight = 0.01;
  Float_t marginLeft = 0.2;
  Float_t marginBottom = 0.15;
  Float_t marginTop = 0.01;
  //#########################################################
  //--------- ULS and LS distributions
  plotting::format1Dhist(vec_xAxisProjections[0], kBlue);
  plotting::format1Dhist(vec_xAxisProjections[1], kRed);
  plotting::format1Dhist(vec_xAxisProjections[2], plotting::kGREEN);
  plotting::format1Dhist(vec_yAxisProjections[0], kBlue);
  plotting::format1Dhist(vec_yAxisProjections[1], kRed);
  plotting::format1Dhist(vec_yAxisProjections[2], plotting::kGREEN);

  TLegend* legUSLS = new TLegend(0.53, 0.58, 0.98, 0.69);
  legUSLS->SetFillStyle(0);
  legUSLS->AddEntry(vec_xAxisProjections[0], "Unlike-sign pairs",       "lep");
  legUSLS->AddEntry(vec_xAxisProjections[1], "Like-sign combinatorial background", "lep");
  legUSLS->AddEntry(vec_xAxisProjections[2], "Raw dielectron signal",   "lep");

  TCanvas* canvXaxisInputs = new TCanvas("canvXaxisInputs", "Raw Spectra");
  canvXaxisInputs->SetWindowSize(900, 900);
  canvXaxisInputs->SetCanvasSize(850, 850);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  vec_xAxisProjections[0]->SetMaximum(vec_xAxisProjections[0]->GetMaximum()*300);
  vec_xAxisProjections[0]->SetTitle("");
  vec_xAxisProjections[0]->SetYTitle(plotting::yieldTotalMassLabel);
  vec_xAxisProjections[0]->SetXTitle(plotting::massAxisLabel);
  vec_xAxisProjections[0]->GetYaxis()->SetTitleOffset(1.7);
  vec_xAxisProjections[0]->Draw("E0");
  vec_xAxisProjections[1]->Draw("SAME E0");
  vec_xAxisProjections[2]->Draw("SAME E0");
  legUSLS->Draw();
  if(prelimPlots){
    texPrelimShifted->Draw();
  }
  for(Int_t i = 0; i < massLabelsShifted.size()-1; ++i){
    massLabelsShifted[i]->Draw();
  }
  TCanvas* canvYaxisInputs = new TCanvas("canvYaxisInputs", "Raw Spectra");
  canvYaxisInputs->SetWindowSize(900, 900);
  canvYaxisInputs->SetCanvasSize(850, 850);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  vec_yAxisProjections[0]->SetMaximum(vec_yAxisProjections[0]->GetMaximum()*3000);
  vec_yAxisProjections[0]->SetTitle("");
  vec_yAxisProjections[0]->SetYTitle(plotting::yieldTotalPairPtLabel);
  vec_yAxisProjections[0]->SetXTitle(plotting::pairPtAxisLabel);
  vec_yAxisProjections[0]->GetYaxis()->SetTitleOffset(1.7);
  vec_yAxisProjections[0]->Draw("E0");
  vec_yAxisProjections[1]->Draw("SAME E0");
  vec_yAxisProjections[2]->Draw("SAME E0");
  legUSLS->Draw();
  if(prelimPlots){
    texPrelimShifted->Draw();
  }
  for(Int_t i = 0; i < pairPtLabelsShifted.size(); ++i){
    pairPtLabelsShifted[i]->Draw();
  }

  //#########################################################
  //-------- Raw signal R factor comparison
  plotting::format1Dhist(vec_xAxisProjections[3], kBlue);
  plotting::format1Dhist(vec_yAxisProjections[3], kBlue);

  TLegend* legRatio = new TLegend(0.5, 0.5, 0.7, 0.7);
  legRatio->AddEntry(vec_xAxisProjections[2], "Full R factor", "lep");
  legRatio->AddEntry(vec_xAxisProjections[3], "No R factor", "lep");

  TCanvas* canvXaxisRfacComp = 0x0;
  TCanvas* canvYaxisRfacComp = 0x0;
  if(hasMixedEvents){
    canvXaxisRfacComp = new TCanvas("canvXaxisRfacComp", "R factor comparison");
    canvXaxisRfacComp->SetWindowSize(900, 900);
    canvXaxisRfacComp->SetCanvasSize(850, 850);
    TPad* padMass1 = new TPad("padMass1", "padMass1", 0, 0.3, 1, 1.0);
    padMass1->SetBottomMargin(0);
    padMass1->SetGridx();
    padMass1->Draw();
    padMass1->cd();
    gPad->SetLogy();
    if(zoomXaxis){
      vec_xAxisProjections[3]->GetXaxis()->SetRangeUser(0, 1.2);
    }
    vec_xAxisProjections[3]->SetMaximum(vec_xAxisProjections[3]->GetMaximum()*100);
    vec_xAxisProjections[3]->SetTitle("R factor effect");
    vec_xAxisProjections[3]->SetYTitle("d#it{N}_{ee}/d#it{m}_{ee} (GeV/#it{c}^{2})^{-1}");
    vec_xAxisProjections[3]->Draw("E0");
    vec_xAxisProjections[2]->Draw("SAME E0");
    legRatio->Draw();
    for(Int_t i = 0; i < massLabelsLeft.size(); ++i){
      massLabelsLeft[i]->Draw();
    }
    TGaxis* axisMass = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    axisMass->SetLabelFont(43);
    axisMass->SetLabelSize(15);
    axisMass->Draw();
    canvXaxisRfacComp->cd();
    TPad* padMass2 = new TPad("padMass2", "padMass2", 0, 0.05, 1, 0.3);
    padMass2->SetTopMargin(0);
    padMass2->SetBottomMargin(0.2);
    padMass2->SetGridx();
    padMass2->Draw();
    padMass2->cd();
    TH1F* rFacMassRatio = (TH1F*)vec_xAxisProjections[2]->Clone("rFacMassRatio");
    rFacMassRatio->Divide(vec_xAxisProjections[3]);
    plotting::formatRatioPlot(rFacMassRatio, "without/with", kBlack);
    rFacMassRatio->SetMinimum(0.95);
    rFacMassRatio->SetMaximum(1.05);
    rFacMassRatio->Draw("HIST");

    canvYaxisRfacComp = new TCanvas("canvYaxisRfacComp", "R factor comparison");
    canvYaxisRfacComp->SetWindowSize(900, 900);
    canvYaxisRfacComp->SetCanvasSize(850, 850);
    TPad* padPairPt1 = new TPad("padPairPt1", "padPairPt1", 0, 0.3, 1, 1.0);
    padPairPt1->SetBottomMargin(0);
    padPairPt1->SetGridx();
    padPairPt1->Draw();
    padPairPt1->cd();
    gPad->SetLogy();
    if(zoomXaxis){
      vec_yAxisProjections[3]->GetXaxis()->SetRangeUser(0, 1.2);
    }
    vec_yAxisProjections[3]->SetMaximum(vec_yAxisProjections[3]->GetMaximum()*100);
    vec_yAxisProjections[3]->SetTitle("R factor effect");
    vec_yAxisProjections[3]->SetYTitle("d#it{N}_{ee}/d#it{m}_{ee} (GeV/#it{c}^{2})^{-1}");
    vec_yAxisProjections[3]->Draw("E0");
    vec_yAxisProjections[2]->Draw("SAME E0");
    legRatio->Draw();
    for(Int_t i = 0; i < pairPtLabelsLeft.size(); ++i){
      pairPtLabelsLeft[i]->Draw();
    }
    TGaxis* axisPairPt = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    axisPairPt->SetLabelFont(43);
    axisPairPt->SetLabelSize(15);
    axisPairPt->Draw();
    canvYaxisRfacComp->cd();
    TPad* padPairPt2 = new TPad("padPairPt2", "padPairPt2", 0, 0.05, 1, 0.3);
    padPairPt2->SetTopMargin(0);
    padPairPt2->SetBottomMargin(0.2);
    padPairPt2->SetGridx();
    padPairPt2->Draw();
    padPairPt2->cd();
    TH1F* rFacPairPtRatio = (TH1F*)vec_yAxisProjections[2]->Clone("rFacPairPtRatio");
    rFacPairPtRatio->Divide(vec_yAxisProjections[3]);
    plotting::formatRatioPlot(rFacPairPtRatio, "without/with", kBlack);
    rFacPairPtRatio->SetMinimum(0.95);
    rFacPairPtRatio->SetMaximum(1.05);
    rFacPairPtRatio->Draw("HIST");
  }

  //#########################################################
  //-------------  Final Spectra
  plotting::format1Dhist(vec_xAxisProjections[4], plotting::kGREEN, plotting::massAxisLabel);
  plotting::format1Dhist(vec_yAxisProjections[4], plotting::kGREEN, plotting::pairPtAxisLabel);

  if(scaleByPhysSel){
    vec_xAxisProjections[4]->SetYTitle(plotting::yieldMassLabel);
    vec_yAxisProjections[4]->SetYTitle(plotting::yieldPairPtLabel);
  }

  TCanvas* canvXaxisSignal = new TCanvas("canvXaxisSignal", "Dielectron Spectrum");
  canvXaxisSignal->SetWindowSize(900, 900);
  canvXaxisSignal->SetCanvasSize(850, 850);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  if(zoomXaxis){
    vec_xAxisProjections[4]->GetXaxis()->SetRangeUser(0, 1.2);
  }
  vec_xAxisProjections[4]->SetTitle("");
  vec_xAxisProjections[4]->GetYaxis()->SetTitleOffset(1.7);
  vec_xAxisProjections[4]->SetMaximum(vec_xAxisProjections[4]->GetMaximum()*10);
  vec_xAxisProjections[4]->Draw("E0");
  if(prelimPlots){
    texPrelimShifted->Draw();
  }
  for(Int_t i = 0; i < massLabelsShifted.size()-1; ++i){
    massLabelsShifted[i]->Draw();
  }
  TCanvas* canvYaxisSignal = new TCanvas("canvYaxisSignal", "Dielectron Spectrum");
  canvYaxisSignal->SetWindowSize(900, 900);
  canvYaxisSignal->SetCanvasSize(850, 850);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  vec_yAxisProjections[4]->SetTitle("");
  vec_yAxisProjections[4]->GetYaxis()->SetTitleOffset(1.7);
  vec_yAxisProjections[4]->SetMaximum(vec_yAxisProjections[4]->GetMaximum()*500);
  vec_yAxisProjections[4]->Draw("E0");
  if(prelimPlots){
    texPrelimShifted->Draw();
  }

  for(Int_t i = 0; i < pairPtLabelsShifted.size()-1; ++i){
    pairPtLabelsShifted[i]->Draw();
  }

  //#########################################################
  //------ SB and significance plots
  plotting::format1Dhist(vec_xAxisProjections[5], kBlue);
  plotting::format1Dhist(vec_xAxisProjections[6], plotting::kGREEN);
  plotting::format1Dhist(vec_yAxisProjections[5], kBlue);
  plotting::format1Dhist(vec_yAxisProjections[6], plotting::kGREEN);

  TLegend* legStats = new TLegend(0.47, 0.58, 0.92, 0.69);
  legStats->SetFillStyle(0);
  legStats->AddEntry(vec_xAxisProjections[5], "S/B", "lep");
  legStats->AddEntry(vec_xAxisProjections[6], "Statistical significance", "lep");

  TCanvas* canvStatsXaxis = new TCanvas("canvStatsXaxis", "canvStatsXaxis");
  canvStatsXaxis->SetWindowSize(900, 900);
  canvStatsXaxis->SetCanvasSize(850, 850);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetBottomMargin(marginBottom);
  vec_xAxisProjections[5]->SetXTitle(plotting::massAxisLabel);
  vec_xAxisProjections[5]->SetTitle("");
  vec_xAxisProjections[5]->SetMinimum(0.001);
  vec_xAxisProjections[5]->SetMaximum(100000);
  vec_xAxisProjections[5]->Draw("E0");
  vec_xAxisProjections[6]->Draw("SAME E0");
  legStats->Draw();
  if(prelimPlots){
    texPrelimLeft->Draw();
  }
  for(Int_t i = 0; i < massLabelsShifted.size()-1; ++i){
    massLabelsLeft[i]->Draw();
  }

  TCanvas* canvStatsYaxis = new TCanvas("canvStatsYaxis", "canvStatsYaxis");
  canvStatsYaxis->SetWindowSize(900, 900);
  canvStatsYaxis->SetCanvasSize(850, 850);
  gPad->SetRightMargin(marginRight);
  gPad->SetBottomMargin(marginBottom);
  gPad->SetLogy();
  vec_yAxisProjections[5]->SetXTitle(plotting::pairPtAxisLabel);
  vec_yAxisProjections[5]->SetTitle("");
  vec_yAxisProjections[5]->SetMinimum(0.0001);
  vec_yAxisProjections[5]->SetMaximum(vec_yAxisProjections[6]->GetMaximum()*10000);
  vec_yAxisProjections[5]->Draw("E0");
  vec_yAxisProjections[6]->Draw("SAME E0");
  legStats->Draw();
  if(prelimPlots){
    texPrelimLeft->Draw();
  }
  for(Int_t i = 0; i < pairPtLabelsShifted.size(); ++i){
    pairPtLabelsLeft[i]->Draw();
  }


  //#########################################################
  //------- R factor plot
  if(hasMixedEvents){
    plotting::format1Dhist(vec_xAxisProjections[7], kBlue);
    plotting::format1Dhist(vec_yAxisProjections[7], kBlue);
  }

  TCanvas* canvRfac = 0x0;
  if(hasMixedEvents){
    canvRfac = new TCanvas("canvRfac", "R factor");
    canvRfac->SetWindowSize(1800, 900);
    canvRfac->SetCanvasSize(1700, 850);
    canvRfac->Divide(3,1);
    canvRfac->cd(1);
    RforComp->SetTitle("R factor");
    RforComp->SetYTitle(plotting::pairPtAxisLabel);
    RforComp->GetZaxis()->SetRangeUser(0.7,1.3);
    RforComp->Draw("COLZ");
    canvRfac->cd(2);
    vec_xAxisProjections[7]->SetTitle("R factor");
    vec_xAxisProjections[7]->SetXTitle(plotting::massAxisLabel);
    vec_xAxisProjections[7]->SetMaximum(1.1);
    vec_xAxisProjections[7]->SetMinimum(0.9);
    vec_xAxisProjections[7]->Draw();
    for(Int_t i = 0; i < massLabelsShifted.size(); ++i){
      massLabelsShifted[i]->Draw();
    }
    canvRfac->cd(3);
    vec_yAxisProjections[7]->SetTitle("R factor");
    vec_yAxisProjections[7]->SetXTitle(plotting::pairPtAxisLabel);
    vec_yAxisProjections[7]->SetMaximum(1.1);
    vec_yAxisProjections[7]->SetMinimum(0.9);
    vec_yAxisProjections[7]->Draw();
    for(Int_t i = 0; i < pairPtLabelsShifted.size(); ++i){
      pairPtLabelsShifted[i]->Draw();
    }
  }

  //#########################################################
  //------------ Efficiency plots
  if(effCorrPath != ""){
    plotting::format1Dhist(vec_xAxisProjections[8], kBlue, plotting::massAxisLabel, "#epsilon_{signal}");
    plotting::format1Dhist(vec_yAxisProjections[8], kBlue, plotting::pairPtAxisLabel, "#epsilon_{signal}");
    if(applyPhiVcorrection){
      plotting::format1Dhist(vec_xAxisProjections[9], kBlack);
      plotting::format1Dhist(vec_yAxisProjections[9], kBlack);
    }
    plotting::format1Dhist(vec_xAxisProjections.back(), plotting::kORANGE, 4);
    plotting::format1Dhist(vec_yAxisProjections.back(), plotting::kORANGE, 4);
  }

  // Create TLegend for different efficiency inputs
  TLegend* legEff = new TLegend(0, 0.18, 0.7, 0.88);
  if(effCorrPath != ""){
    legEff->AddEntry(vec_xAxisProjections[8], "Track Cut Efficiency (Dielectron Framework)", "lep");
    if(applyPhiVcorrection){
      legEff->AddEntry(vec_xAxisProjections[9], "Phi V Cut Efficiency", "lep");
    }
    legEff->AddEntry(vec_xAxisProjections.back(), "Effective Efficiency", "lep");
    legEff->SetTextSize(0.04);
  }

  TCanvas* canvasEff2D = 0x0;
  TCanvas* canvasEff1D = 0x0;
  if(effCorrPath != ""){
     // 1D efficiency projections
    canvasEff1D = new TCanvas("canvasEff1D", "canvasEff1D");
    canvasEff1D->SetWindowSize(1800, 900);
    canvasEff1D->SetCanvasSize(1750, 850);
    canvasEff1D->Divide(3,1);
    canvasEff1D->cd(1);
    vec_yAxisProjections[8]->SetTitle("Efficiencies");
    vec_yAxisProjections[8]->SetMinimum(0);
    vec_yAxisProjections[8]->SetMaximum(1.1);
    vec_yAxisProjections[8]->Draw("E0");
    if(applyPhiVcorrection){
      vec_yAxisProjections[9]->Draw("SAME E0");
    }
    vec_yAxisProjections.back()->Draw("SAME E0");
    canvasEff1D->cd(2);
    vec_xAxisProjections[8]->SetTitle("Efficiencies");
    vec_xAxisProjections[8]->SetMinimum(0);
    vec_xAxisProjections[8]->SetMaximum(1.1);
    vec_xAxisProjections[8]->Draw("E0");
    if(applyPhiVcorrection){
      vec_xAxisProjections[9]->Draw("SAME E0");
    }
    vec_xAxisProjections.back()->Draw("SAME E0");
    canvasEff1D->cd(3);
    legEff->Draw();

    // 2D efficiency plots
    canvasEff2D = new TCanvas("canvasEff2D", "canvasEff2D");
    canvasEff2D->SetWindowSize(1800, 1100);
    canvasEff2D->SetCanvasSize(1750, 1000);
    canvasEff2D->Divide(3,2);
    canvasEff2D->cd(1);
    gPad->SetLogz();
    plotting::Set2DTitles(vec_effInputs[0], "Generated Pairs");
    vec_effInputs[0]->SetMaximum(vec_effInputs[0]->GetMaximum()*1);
    if(restrictEffPlots){
      vec_effInputs[0]->GetXaxis()->SetRangeUser(minMassEffPlots, maxMassEffPlots);
      vec_effInputs[0]->GetYaxis()->SetRangeUser(minPairPtEffPlots, maxPairPtEffPlots);
    }
    vec_effInputs[0]->Draw(effPlotStyle);
    canvasEff2D->cd(2);
    gPad->SetLogz();
    plotting::Set2DTitles(vec_effInputs[1], "Pairs After Cuts");
    vec_effInputs[1]->SetMaximum(vec_effInputs[0]->GetMaximum()*1);
    if(restrictEffPlots){
      vec_effInputs[1]->GetXaxis()->SetRangeUser(minMassEffPlots, maxMassEffPlots);
      vec_effInputs[1]->GetYaxis()->SetRangeUser(minPairPtEffPlots, maxPairPtEffPlots);
    }
    vec_effInputs[1]->Draw(effPlotStyle);
    canvasEff2D->cd(3);
    if(vec_effInputs[2]->GetMaximum() > 1){
      printError("WARNING: Efficiency greater than 1!!");
    }
    else{
      vec_effInputs[2]->SetMaximum(1.);
    }
    plotting::Set2DTitles(vec_effInputs[2], "Pair Efficiency From DiElec Framework");
    vec_effInputs[2]->SetMaximum(1);
    vec_effInputs[2]->SetMinimum(0);
    if(restrictEffPlots){
      vec_effInputs[2]->GetXaxis()->SetRangeUser(minMassEffPlots, maxMassEffPlots);
      vec_effInputs[2]->GetYaxis()->SetRangeUser(minPairPtEffPlots, maxPairPtEffPlots);
    }
    vec_effInputs[2]->Draw(effPlotStyle);
    if(applyPhiVcorrection){
      canvasEff2D->cd(4);
      plotting::Set2DTitles(vec_effInputs[3], "#varphi_{V} efficiency");
      vec_effInputs[3]->Draw(effPlotStyle);
    }
  }

  watch->Stop();
  watch->Print();

  // Alter some save names based on presence of efficiency correction and mixed
  // events
  if(!hasMixedEvents && effCorrPath != ""){
    saveName[2] = "rawSignalNoR";
    saveName[4] = "correctedSignalNoR";
  }
  else if(!hasMixedEvents && effCorrPath == ""){
    saveName[2] = "rawSignalNoR";
    saveName[4] = "rawSignalScaledNoR";
  }
  else if(hasMixedEvents && !effCorrPath){
    saveName[4]= "rawSignalScaled";
  }

  TFile* outFile = 0x0;
  if(saveLocation != ""){
    saveLocation.Prepend(paths::results_loc + paths::spectra_dir + "corrected_spectra/diElec_framework/");
    saveLocation.Append(".root");
    outFile = TFile::Open(saveLocation ,"RECREATE");
      if(!outFile){
        printError("File could not be opened for writing!");
      }
      else{
        outFile->cd();
        // Save corrected spectra in 2D
        spectrum->Write("spectrum2D_corrected");
        // Save raw 2D inputs
        for(Int_t i = 0; i < vec_inputSpectra2D.size(); ++i){
          vec_inputSpectra2D[i]->Write();
        }
        // Save projection plots
        for(Int_t i = 0; i < vec_xAxisProjections.size(); ++i){
          // Skip R factor if not present (only contains NULL pointer)
          if(!hasMixedEvents && i == 7){ continue;}
          vec_xAxisProjections[i]->Write(saveName[i] + "_mass");
          vec_yAxisProjections[i]->Write(saveName[i] + "_pairPt");
        }
        // Save specific canvases
        if(hasMixedEvents){
          canvXaxisRfacComp->Write("canvRfacImpact_mass");
          canvYaxisRfacComp->Write("canvRfacImpact_mass");
        }
        if(effCorrPath != ""){
          canvasEff2D->Write("canvEff2D");
          canvasEff1D->Write("canvEff1D");
        }
        std::cout << "Plots saved to -> " << saveLocation << std::endl;
      }
    }

  return;
}

void setMultClass(Int_t whichMultClass, TString& phiVnoCutsHistName, TString& phiVCutsHistName, Float_t* zAxisRange){

  if(whichMultClass == 0 || whichMultClass == 5){
    std::cout << "Looking at multiplicity range: 00 - 100 %" << std::endl;
    phiVCutsHistName.Append("00100");
    phiVnoCutsHistName.Append("00100");
    zAxisRange[0] = 0;
    zAxisRange[1] = 100;
  }else if(whichMultClass == 1){
    std::cout << "Looking at multiplicity range: 00 - 20 %" << std::endl;
    phiVCutsHistName.Append("0020");
    phiVnoCutsHistName.Append("0020");
    zAxisRange[0] = 0;
    zAxisRange[1] = 20;
  }else if(whichMultClass == 2){
    std::cout << "Looking at multiplicity range: 20 - 40 %" << std::endl;
    phiVCutsHistName.Append("2040");
    phiVnoCutsHistName.Append("2040");
    zAxisRange[0] = 20;
    zAxisRange[1] = 40;
  }else if(whichMultClass == 3){
    std::cout << "Looking at multiplicity range: 40 - 60 %" << std::endl;
    phiVCutsHistName.Append("4060");
    phiVnoCutsHistName.Append("4060");
    zAxisRange[0] = 40;
    zAxisRange[1] = 60;
  }else if(whichMultClass == 4){
    std::cout << "Looking at multiplicity range: 60 - 100 %" << std::endl;
    phiVCutsHistName.Append("60100");
    phiVnoCutsHistName.Append("60100");
    zAxisRange[0] = 60;
    zAxisRange[1] = 100;
  }
    return;
  }

std::vector<Int_t> getBinRanges(Float_t xRange[], Float_t yRange[], Float_t zRange[], const TH3F* inputSpectrum){

  std::vector<Int_t> binRanges;

  //##########################################################
  // ------------  Get range restrictions
  binRanges.push_back(inputSpectrum->GetXaxis()->FindBin(xRange[0]+0.00001));
  binRanges.push_back(inputSpectrum->GetXaxis()->FindBin(xRange[1]-0.00001));
  binRanges.push_back(inputSpectrum->GetYaxis()->FindBin(yRange[0]+0.00001));
  binRanges.push_back(inputSpectrum->GetYaxis()->FindBin(yRange[1]-0.00001));
  binRanges.push_back(inputSpectrum->GetZaxis()->FindBin(zRange[0]+0.00001));
  binRanges.push_back(inputSpectrum->GetZaxis()->FindBin(zRange[1]-0.00001));

  std::cout << "X axis range: " << xRange[0] << ", " << xRange[1];
  std::cout << ", Bins: " << binRanges[0] << ", " << binRanges[1] << std::endl;
  std::cout << "Y axis range: " << yRange[0] << ", " << yRange[1];
  std::cout << ", Bins: " << binRanges[2] << ", " << binRanges[3] << std::endl;
  std::cout << "Z axis range: " << zRange[0] << ", " << zRange[1];
  std::cout << ", Bins: " << binRanges[4] << ", " << binRanges[5] << std::endl;

  return binRanges;
}

TH1F* getRfactorProjections(const std::vector<TH2F*> vec_inputs, const std::vector<Int_t> binRanges, Int_t axis = 0){

  //TODO: remove hard coded numbers
  TH1F* tempPN = 0x0;
  TH1F* tempNN = 0x0;
  TH1F* tempPP = 0x0;
  TH1F* rFacProjection = 0x0;
  if(axis == 0){
    tempPN = (TH1F*)vec_inputs[3]->ProjectionX("tempPN", binRanges[2], binRanges[3], "e");
    tempPP = (TH1F*)vec_inputs[4]->ProjectionX("tempPP", binRanges[2], binRanges[3], "e");
    tempNN = (TH1F*)vec_inputs[5]->ProjectionX("tempNN", binRanges[2], binRanges[3], "e");
  }
  else if(axis == 1){
    tempPN = (TH1F*)vec_inputs[3]->ProjectionY("tempPN", binRanges[0], binRanges[1], "e");
    tempPP = (TH1F*)vec_inputs[4]->ProjectionY("tempPP", binRanges[0], binRanges[1], "e");
    tempNN = (TH1F*)vec_inputs[5]->ProjectionY("tempNN", binRanges[0], binRanges[1], "e");
  }
  if(!tempPN || !tempNN || !tempPP){
    printError("Projections for R factor failed.");
    return 0x0;
  }else{
    rFacProjection = calcDiElecRfactor(tempPN, tempPP, tempNN, (axis == 0 ? 5 : 10));
    return rFacProjection;
  }

}

TH1F* getEffProjections(const std::vector<TH2F*> vec_inputs, const std::vector<Int_t> binRanges, Int_t axis = 0){

  TH1F* effGen = 0x0;
  TH1F* effAcc = 0x0;
  if(axis == 0){
    effGen = (TH1F*)vec_inputs[0]->ProjectionX("effGen", binRanges[2], binRanges[3], "e");
    effAcc = (TH1F*)vec_inputs[1]->ProjectionX("effAcc", binRanges[2], binRanges[3], "e");
  }else if(axis == 1){
    effGen = (TH1F*)vec_inputs[0]->ProjectionY("effGen", binRanges[0], binRanges[1], "e");
    effAcc = (TH1F*)vec_inputs[1]->ProjectionY("effAcc", binRanges[0], binRanges[1], "e");
  }

  TH1F* effDielec = (TH1F*)effGen->Clone("effDielec");
  effDielec->Divide(effAcc, effGen, 1, 1, "B");

  delete effGen;
  delete effAcc;

  return effDielec;
}
