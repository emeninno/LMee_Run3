#include "../../paths/path_locations.C"
#include "../../generalFunctions/random_functions.C"
#include "../../generalFunctions/readInOut.C"
#include "../../generalFunctions/cuts.C"
#include "../../generalFunctions/plotting.C"
#include "../../generalFunctions/diElecFuncs.C"
#include "../../generalFunctions/runArrays.C"
#include "../../generalFunctions/texLabels.C"

void setMultClass(Int_t whichMultClass, TString& phiVnoCutsHistName, TString& phiVCutsHistName, Float_t* zAxisRange);
std::vector<Int_t> getBinRanges(Float_t xRange[], Float_t yRange[], Float_t zRange[], const TH3D* inputSpectrum);
std::vector<Int_t> getBinRanges(Float_t xRange[], Float_t yRange[], const TH2D* inputSpectrum);
TH1D* getRfactorProjections(const std::vector<TH2D*> vec_inputs, const std::vector<Int_t> binRanges, Int_t axis = 0);
TH1D* getEffProjections(const std::vector<TH2D*> vec_inputs, const std::vector<Int_t>& binRanges, Int_t axis = 0, Int_t whichProj = 0);

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Takes spectrum and efficiency correction inputs from the dielectron      //
//  framework, as well as a phiV correction and weighting file which are     //
//  determined locally from TTrees. Produces the following plots in 2D, and  //
//  projected onto the mass and pairPt axes:                                 //
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
//    - effCorrLF:      e.g "cent/18f3/2018_08_27_kCutSet1_00100"            //
//    - effCorrHF:      e.g "cent/19h9h/2018_08_27_kCutSet1_00100"           //
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
//                                                                           //
//                                                                           //
//    Important inputs set in `script settings`                              //
//    - whichWeights:   an output from 'cocktail/createWeightingFile.C'      //
//                      which weights the LF and HF contributions for the    //
//                      efficiency correction according to the ratio from    //
//                      the cocktail.                                        //
//    - phiVpath:       phiV correction, obtained with                       //
//                      'convRejections/getPhiVcutEffMultDep.C'              //
//                         - leave empty and no phiV correction will be      //
//                           applied                                         //
//                                                                           //
// Analysis specific customisations can be found in 'script settings' below. //
///////////////////////////////////////////////////////////////////////////////


void plotSpectra(TString whichCutSet,
                 TString filesPath,
                 TString effCorrLF,
                 TString effCorrHF = "",
                 TString saveName = "",
                 Int_t whichMultClass = 0)
{

  // ------------------------- Script settings ------------------------------//
  // - Dielectron output folder, cut setting,  and 3D histogram name
  std::vector<TString> diElecOutputNames  = {"acapon_out_0", whichCutSet, "pInvMass_PairPt_CentralityV0A"};
  // - Dielectron efficiency output folder, cut setting
  std::vector<TString> diElecEffNames = {"efficiency0", whichCutSet};

  // Which cocktail weighting file to use
  TString whichWeights = "POWHEG";

  TString phiVpath        = "merged/18f3/2020_03_18_phiVeff";
  // Post correction to the light flavour efficiency
  // Will decrease LF efficiency by 3%
  Bool_t adHocLF = kTRUE;

  // PhiV cut settings
  // If phiV cut was applied during train trains (i.e. hard coded into histogram)
  Double_t massPhiVcut = 0.14;

  // -----  Rebinning options
  // Data and efficiency bins do not need to match.
  // All rebinning done before any operations.
  // If rebin flag(s) set and vectors empty, standard rebin will be used
  //  - Data
  const Bool_t rebin2D = kTRUE;
  std::vector<Double_t> newMassBins   = plotting::vec_massBinsPPprelim;
  std::vector<Double_t> newPairPtBins = plotting::vec_pairPtBinsPPprelim;
  const Int_t numRebin[2] = {1,4}; // Standard rebin values
  // - Efficiency
  const Bool_t rebinEff = kTRUE;
  std::vector<Double_t> newMassBinsEff   = plotting::vec_massBinsPPprelim;
  std::vector<Double_t> newPairPtBinsEff = plotting::vec_pairPtBinsPPprelim;
  const Int_t numEffRebin[2] = {2,10};

  Bool_t rebinPhiVeff = kFALSE;

  // Apply phiV cut
  // Set to zero if not to be used (check input histogram too!)
  const Float_t phiVmassCutNew = 0;
  const Float_t phiVcutNew = 0;

  // Point at which the R factor is set to 1 (mass axis)
  const Float_t RfacOne = 0.3;
  // Scale by number of events
  const Bool_t scaleByPhysSel  = kTRUE;
  // Scale final spectra (excluding stat. sig. SB and R factor) by bin width
  const Bool_t scaleByBinWidth = kTRUE;
  // Scale final statistical significance plot by 1/sqrt{N}
  Bool_t scaleStatSigByEvents = kFALSE;

  // Restrict X and Y axis ranges before projections
  Float_t xAxisRange[2] = {0, 3.5};
  Float_t yAxisRange[2] = {0, 8};

  // - - - - Plotting settings
  // Restrict 2D X and Y ranges (will match projection ranges)
  const Bool_t restrict2DplotRanges = kTRUE;
  const TString effPlotStyle = "COLZ";

  // Plot extra latex title on most plots (e.g. "ALICE Preliminary")
  const TString latexTitle = "This Thesis";
  // ------------------------- Script settings ------------------------------//

  TStopwatch* watch = new TStopwatch();
  watch->Start();


  if(!effCorrLF.IsNull() && !effCorrHF.IsNull() && whichWeights.IsNull()){
    printError("Must specify cocktail weighting file if using LF and HF");
    return;
  }

  // Ensure correct phiV cuts and correction are being applied
  if((phiVcutNew != 0 || phiVmassCutNew != 0) && !diElecOutputNames.at(2).Contains("Phiv")){
    printError("Applying post phiV cut but not importing histogram with phiV");
    return;
  }


  // Use much wider binning if looking multiplicity classes
  /* if(whichMultClass != 0){ */
  /*   std::cout << "Setting multiplicity bins!" << std::endl; */
  /*   newMassBins      = plotting::vec_massBinsMultDep; */
  /*   newPairPtBins    = plotting::vec_pairPtBinsMultDep; */
  /*   newMassBinsEff   = plotting::vec_massBinsMultDep; */
  /*   newPairPtBinsEff = plotting::vec_pairPtBinsMultDep; */
  /* } */

  // Set PhiV correction flag
  Bool_t usePhiVcorr = !phiVpath.IsNull();

  // Function returns following variables with correct multiplicity labels appended
  // and correct multiplicity range set
  TString phiVnoCutsHistName = "noCut";
  TString phiVcutsHistName   = "cut";
  Float_t zAxisRange[2] = {-99, -99};
  setMultClass(whichMultClass, phiVnoCutsHistName, phiVcutsHistName, zAxisRange);

  // Number of events after cuts stored here
  Double_t numEvents = 0.0;

  // Get input spectra (same event and mixed event(me) )
  // Vector contents: {ULS, PosLS, NegLS, ULSme, PosLSme, NegLSme}
  Bool_t hasMixedEvents = kTRUE;
  std::vector<TH3D*> vec_inputSpectra = readInOut::getInputSpectra(paths::results_loc+paths::diElec_dir+"data/"+filesPath+".root", diElecOutputNames, numEvents);
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
  std::vector<Int_t> binRanges = getBinRanges(xAxisRange, yAxisRange, zAxisRange, vec_inputSpectra.at(0));

  // Apply non-standard phiV cut
  if(phiVcutNew != 0 || phiVmassCutNew != 0){
    TString phiVmassString;
    phiVmassString.Form("%01.2f", phiVmassCutNew);
    TString phiVcutString;
    phiVcutString.Form("%01.2f", phiVcutNew);
    printInfo("Manually applying phiV cut: mass < "+ phiVmassString + " && phiV > " + phiVcutString);

    // Set content of bins in phiV cut range to zero
    for(Int_t i = 0; i < vec_inputSpectra.size(); ++i){

      Float_t maxMassBin = vec_inputSpectra.at(0)->GetXaxis()->FindBin(phiVmassCutNew-0.000001);
      Float_t maxPhiVbin = vec_inputSpectra.at(0)->GetZaxis()->FindBin(phiVcutNew+0.000001);

      for(Int_t binX = 0; binX <= maxMassBin; ++binX){
        for(Int_t binY = 0; binY <= vec_inputSpectra.at(0)->GetYaxis()->GetNbins(); ++binY){
          for(Int_t binZ = maxPhiVbin; binZ <= vec_inputSpectra.at(i)->GetZaxis()->GetNbins(); ++binZ){

            vec_inputSpectra.at(i)->SetBinContent(binX, binY, binZ, 0);
            vec_inputSpectra.at(i)->SetBinError(binX,  binY,  binZ, 0);
          }
        }
      }
    }

    // Set correct bin ranges (otherwise script assumes multiplicity axis)
    binRanges[4] = 0;
    binRanges[5] = -1;
  }// End post phiV cut

  // Restrict Z axis before projection and store projections
  std::vector<TH2D*> vec_inputSpectra2D;
  for(Int_t i = 0; i < vec_inputSpectra.size(); ++i){
    vec_inputSpectra.at(i)->GetZaxis()->SetRange(binRanges[4], binRanges[5]);
    TH2D* histProj = (TH2D*)vec_inputSpectra.at(i)->Project3D("yxe");
    vec_inputSpectra2D.push_back(histProj);
  }

  if(rebin2D){
    if(newMassBins.size() == 0){
      std::cout << "Performing simple 2D rebin on all input spectra" << std::endl;
      std::cout << "X axis rebin: " << numRebin[0] << ", y-axis rebin: " << numRebin[1] << std::endl;
      for(Int_t i = 0; i < vec_inputSpectra2D.size(); ++i){
        vec_inputSpectra2D.at(i)->Rebin2D(numRebin[0], numRebin[1]);
      }
    }
    else if(newMassBins.size() > 1){
      std::cout << "Performing custom 2D rebin on all input spectra" << std::endl;
      for(Int_t i = 0; i < vec_inputSpectra2D.size(); ++i){
        vec_inputSpectra2D.at(i) = plotting::rebin2Dhist(vec_inputSpectra2D.at(i), newMassBins, newPairPtBins);
        if(vec_inputSpectra2D.at(i) == 0x0){
          return;
        }
      }
    }
    // If rebinning, find new bins for projections
    std::cout << "Bin ranges after rebinning" << std::endl;
    binRanges = getBinRanges(xAxisRange, yAxisRange, vec_inputSpectra2D.at(0));
  }

  //----------- Get generated smeared and accepted histograms
  // Vector contents: {generated, accepted, diElecEfficiency}
  std::vector<TH2D*> vec_effInputs;
  std::vector<TH2D*> vec_effInputsLF;
  std::vector<TH2D*> vec_effInputsHF;
  std::vector<TH2D*> vec_weights;
  Bool_t hasEffCorr = kFALSE;
  if(!effCorrLF.IsNull() && !effCorrHF.IsNull()){

    hasEffCorr = kTRUE;

    TString effFilePath = paths::results_loc + paths::diElec_dir + "effCorrs/";

    vec_effInputsLF = readInOut::getEfficiencySpectra(effFilePath+effCorrLF+".root", diElecEffNames, 1);
    if(vec_effInputsLF.size() != 2){
      return;
    }
    vec_effInputsHF = readInOut::getEfficiencySpectra(effFilePath+effCorrHF+".root", diElecEffNames, 2);
    if(vec_effInputsHF.size() != 2){
      return;
    }

    if(rebinEff){
      if(newMassBinsEff.size() == 0){
        std::cout << "Standard rebinning applied to efficiency histograms" << std::endl;
        for(Int_t i = 0; i < vec_effInputsLF.size(); ++i){
          vec_effInputsLF.at(0)->Rebin2D(numEffRebin[0], numEffRebin[1]);
          vec_effInputsLF.at(1)->Rebin2D(numEffRebin[0], numEffRebin[1]);
          vec_effInputsHF.at(0)->Rebin2D(numEffRebin[0], numEffRebin[1]);
          vec_effInputsHF.at(1)->Rebin2D(numEffRebin[0], numEffRebin[1]);
        }
      }
      else if(newMassBinsEff.size() > 1){
        for(Int_t i = 0; i < vec_effInputsLF.size(); ++i){
          vec_effInputsLF.at(i) = plotting::rebin2Dhist(vec_effInputsLF.at(i), newMassBinsEff, newPairPtBinsEff);
          vec_effInputsHF.at(i) = plotting::rebin2Dhist(vec_effInputsHF.at(i), newMassBinsEff, newPairPtBinsEff);
          if(vec_effInputsLF.at(i) == 0x0 || vec_effInputsHF.at(i) == 0x0){
            printError("Custom rebin of efficiency plot failed");
            return;
          }
        }
      }
    }

    // Use weighting file to weight sameMother and HF contributions
    if(!whichWeights.IsNull() && !effCorrLF.IsNull() && !effCorrHF.IsNull()){

      TFile* inWeightFile = TFile::Open(paths::results_loc + paths::cocktail_dir +
                                        "weightingFiles/" + whichWeights + ".root", "READ");
      if(!inWeightFile){
        printError("Weighting file not found!");
        return;
      }

      TH2D* frac_LF = (TH2D*)inWeightFile->Get("fraction_sameMother");
      TH2D* frac_HF = (TH2D*)inWeightFile->Get("fraction_HF");
      if(!frac_LF || !frac_HF){
        printError("Weighting histogram not found!");
        return;
      }


      vec_weights.push_back(frac_LF);
      vec_weights.push_back(frac_HF);

      // General purpose Monte Carlos that were used for pPb diElec analysis
      // suffer from incorrect SPD dead maps and need adhoc 3% correction.
      if(adHocLF){
        printInfo("Scaling LF efficiency by 3% (before merging)");
        for(Int_t i = 0; i < vec_effInputsLF.at(1)->GetNbinsX(); ++i){
          for(Int_t j = 0; j < vec_effInputsLF.at(1)->GetNbinsY(); ++j){
            vec_effInputsLF.at(1)->SetBinContent(i, j, vec_effInputsLF.at(1)->GetBinContent(i, j)*0.97);
            vec_effInputsLF.at(1)->SetBinError(i, j, TMath::Sqrt(vec_effInputsLF.at(1)->GetBinContent(i, j)));
          }
        }
      }

      // Create temporary source efficiencies
      TH2D* tempEffLF = (TH2D*)vec_effInputsLF.at(0)->Clone("tempEffLF");
      tempEffLF->Divide(vec_effInputsLF.at(1), vec_effInputsLF.at(0), 1, 1, "B");
      TH2D* tempEffHF = (TH2D*)vec_effInputsHF.at(0)->Clone("tempEffHF");
      tempEffHF->Divide(vec_effInputsHF.at(1), vec_effInputsHF.at(0), 1, 1, "B");

      // Weight contributions
      for(Int_t i = 1; i <= tempEffLF->GetNbinsX(); ++i){

        Float_t xBinCenter = tempEffLF->GetXaxis()->GetBinCenter(i);

        for(Int_t j = 1; j <= tempEffLF->GetNbinsY(); ++j){

          Float_t yBinCenter = tempEffLF->GetYaxis()->GetBinCenter(j);

          // Find corresponding bins in weights histograms
          Int_t binWeightingX = frac_LF->GetXaxis()->FindBin(xBinCenter);
          Int_t binWeightingY = frac_LF->GetYaxis()->FindBin(yBinCenter);

          Float_t LFfrac = frac_LF->GetBinContent(binWeightingX, binWeightingY);
          Float_t HFfrac = frac_HF->GetBinContent(binWeightingX, binWeightingY);

          if(LFfrac == 0){
            HFfrac = 1;
          }else if(HFfrac == 0){
            LFfrac = 1;
          }

          tempEffLF->SetBinContent(i, j, tempEffLF->GetBinContent(i, j)*LFfrac);
          tempEffLF->SetBinError(i,   j, tempEffLF->GetBinError(i, j)*LFfrac);

          tempEffHF->SetBinContent(i, j, tempEffHF->GetBinContent(i, j)*HFfrac);
          tempEffHF->SetBinError(i,   j, tempEffHF->GetBinError(i, j)*HFfrac);
        }
      }
      TH2D* weightedEff = (TH2D*)tempEffLF->Clone("histEffDielec2D");
      weightedEff->Add(tempEffHF);

      vec_effInputs.push_back(weightedEff);
      vec_effInputs.push_back(vec_effInputsLF.at(0));
      vec_effInputs.push_back(vec_effInputsLF.at(1));
      vec_effInputs.push_back(vec_effInputsHF.at(0));
      vec_effInputs.push_back(vec_effInputsHF.at(1));

    }
  }else{ // Push NULL pointer so plotting functions work
    vec_effInputs.push_back(0x0);
    printInfo("No efficiency correction applied. All results raw.");
  }// End efficiency retrieval section

  //##########################################################
  // Get phiV cut efficiency
  TH2D* phiVeff = 0x0;
  // Projection plots
  TH1D* phiVmassEff   = 0x0;
  TH1D* phiVpairPtEff = 0x0;
  if(usePhiVcorr && vec_effInputs.size() != 0){
    TFile* inPhiVfile = TFile::Open(paths::results_loc + paths::conversionRej_dir+ phiVpath + ".root", "READ");
    if(!inPhiVfile){
      printError("PhiV eff correction file not found.");
      return;
    }
    TH2D* phiVnoCut = (TH2D*)inPhiVfile->Get(phiVnoCutsHistName);
    TH2D* phiVcut   = (TH2D*)inPhiVfile->Get(phiVcutsHistName);
    if(!phiVnoCut | !phiVcut){
      printError("PhiV input plot not found");
      return;
    }

    if(rebinPhiVeff){
      printInfo("Rebinning phiV efficiency. Desired behaviour?!?!?!?!?!");
      phiVcut   = plotting::rebin2Dhist(phiVcut, newMassBins, newPairPtBins);
      phiVnoCut = plotting::rebin2Dhist(phiVnoCut, newMassBins, newPairPtBins);
      if(phiVcut == 0x0 || phiVnoCut == 0x0){
        printError("Custom rebin of phiV failed!");
        return;
      }
    }
    phiVeff = (TH2D*)phiVcut->Clone("phiVeff");
    phiVeff->Reset();
    phiVeff->Divide(phiVcut, phiVnoCut, 1, 1, "B");

    // PhiV projections
    TH1D* phiVmassCutTemp   = (TH1D*)phiVcut->ProjectionX("phVmassNumTemp",     binRanges[2], binRanges[3], "e");
    TH1D* phiVmassNoCutTemp = (TH1D*)phiVnoCut->ProjectionX("phVmassNoCutTemp", binRanges[2], binRanges[3], "e");

    TH1D* phiVpairPtCutTemp   = (TH1D*)phiVcut->ProjectionY("phVpairPtNumTemp",     binRanges[0], binRanges[1], "e");
    TH1D* phiVpairPtNoCutTemp = (TH1D*)phiVnoCut->ProjectionY("phVpairPtNoCutTemp", binRanges[0], binRanges[1], "e");

    phiVmassEff = (TH1D*)phiVmassCutTemp->Clone("phiVmassEff");
    phiVmassEff->Clear();
    phiVmassEff->Divide(phiVmassCutTemp, phiVmassNoCutTemp, 1, 1, "e");

    phiVpairPtEff = (TH1D*)phiVpairPtCutTemp->Clone("phiVpairPtEff");
    phiVpairPtEff->Reset();
    phiVpairPtEff->Divide(phiVpairPtCutTemp, phiVpairPtNoCutTemp, 1, 1, "B");
    // If looking at mass ranges above 100 GeV, set phiV on pairPt axis to 1
    // without error
    if(xAxisRange[0] > massPhiVcut){
      for(Int_t i = 0; i < phiVpairPtEff->GetNbinsX(); ++i){
        phiVpairPtEff->SetBinContent(i, 1);
        phiVpairPtEff->SetBinError(i, 0);
      }
    }

    // Set all bins above 100 MeV to 1 with 0 error
    Int_t lastPhiVbin = phiVeff->GetXaxis()->FindBin(massPhiVcut + 0.0000000000001);
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
  TH2D* R              = 0x0;
  TH2D* RforComp       = 0x0;
  if(hasMixedEvents){
    R        = (TH2D*)calcDiElecRfactor(vec_inputSpectra2D.at(3), vec_inputSpectra2D.at(4), vec_inputSpectra2D.at(5), RfacOne);
    RforComp = (TH2D*)calcDiElecRfactor(vec_inputSpectra2D.at(3), vec_inputSpectra2D.at(4), vec_inputSpectra2D.at(5), 5);
  }

  TH2D* combBackgr     = (TH2D*)calcDiElecBackgr(vec_inputSpectra2D.at(1), vec_inputSpectra2D.at(2));
  TH2D* spectrum       = (TH2D*)calcDiElecSpectrum(vec_inputSpectra2D.at(0), combBackgr, R, vec_effInputs.front());
  if(usePhiVcorr && hasEffCorr){
    applyEffCorr(spectrum, phiVeff);
    vec_effInputs.push_back(phiVeff);
  }


  TH2D* rawSpectrumNoR = (TH2D*)calcDiElecSpectrum(vec_inputSpectra2D.at(0), combBackgr, 0x0, 0x0);
  TH2D* rawSpectrum    = (TH2D*)calcDiElecSpectrum(vec_inputSpectra2D.at(0), combBackgr, RforComp, 0x0);
  // 1D projections will be saved from the two chosen axes under the following
  // names. Axis name appended (currently assuming mass and pairPt)
   // Create vector containing final projections
  // Vector contents: {ULS, combBackgr, raw signal, raw signal no R factor, raw
  // signal, SB ratio, statistical significance, R factor}
  // Projections onto X axis
  std::vector<TH1D*> vec_xAxisProjections;
  vec_xAxisProjections.push_back((TH1D*)vec_inputSpectra2D.at(0)->ProjectionX("ULSxAxisProj", binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1D*)combBackgr->ProjectionX("combBackgrXaxisProj",     binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1D*)rawSpectrum->ProjectionX("rawSigXaxisProj",        binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1D*)rawSpectrumNoR->ProjectionX("rawSigNoRxAxisProj",  binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back((TH1D*)spectrum->ProjectionX("signalXaxisProj",           binRanges[2], binRanges[3], "e"));
  vec_xAxisProjections.push_back(calcDiElecSB(vec_xAxisProjections.at(2), vec_xAxisProjections.at(1)));
  vec_xAxisProjections.push_back(calcDiElecSignificance(vec_xAxisProjections.at(2), vec_xAxisProjections.at(1)));
  if(scaleStatSigByEvents){
    vec_xAxisProjections.back()->Scale(1./TMath::Sqrt(numEvents));
  }
  if(hasMixedEvents){
    vec_xAxisProjections.push_back(getRfactorProjections(vec_inputSpectra2D, binRanges, 0)); // 7th object in vector
  }else{
    vec_xAxisProjections.push_back(0x0);
  }


  // Projections onto Y axis
  std::vector<TH1D*> vec_yAxisProjections;
  vec_yAxisProjections.push_back((TH1D*)vec_inputSpectra2D.at(0)->ProjectionY("ULSyAxisProj", binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1D*)combBackgr->ProjectionY("combBackgrYaxisProj",     binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1D*)rawSpectrum->ProjectionY("rawSigYaxisProj",        binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1D*)rawSpectrumNoR->ProjectionY("rawSigNoRyAxisProj",  binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back((TH1D*)spectrum->ProjectionY("signalYaxisProj",           binRanges[0], binRanges[1], "e"));
  vec_yAxisProjections.push_back(calcDiElecSB(vec_yAxisProjections.at(2), vec_yAxisProjections.at(1)));
  vec_yAxisProjections.push_back(calcDiElecSignificance(vec_yAxisProjections.at(2), vec_yAxisProjections.at(1)));
  if(scaleStatSigByEvents){
    vec_yAxisProjections.back()->Scale(1./TMath::Sqrt(numEvents));
  }
  if(hasMixedEvents){
    vec_yAxisProjections.push_back(getRfactorProjections(vec_inputSpectra2D, binRanges, 1)); // 7th object in vector
  }
  else{
    vec_yAxisProjections.push_back(0x0);
  }

  // Vector to store plot names. Will be added to depending on which files are used.
  std::vector<TString> plotSaveNames= {"ULS", "combBackgr", "rawSignal", "rawSigNoR",
                                       "correctedSignal", "SB", "significance", "Rfac"};

  //###########################################################3
  // Push efficiency projections into existing projection vectors
  // Contents pushed into vector (in final order):
  //          {diElecEffLF, diElecEffHF, phiVeff, effective efficiency}
  if(hasEffCorr){
    if(!effCorrLF.IsNull()){
      vec_xAxisProjections.push_back(getEffProjections(vec_effInputs, binRanges, 0, 0));
      vec_yAxisProjections.push_back(getEffProjections(vec_effInputs, binRanges, 1, 0));
      plotSaveNames.push_back("diElecEff_LF");
    }
    if(!effCorrHF.IsNull()){
      vec_xAxisProjections.push_back(getEffProjections(vec_effInputs, binRanges, 0, 1));
      vec_yAxisProjections.push_back(getEffProjections(vec_effInputs, binRanges, 1, 1));
      plotSaveNames.push_back("diElecEff_HF");
    }
    if(usePhiVcorr){
      vec_xAxisProjections.push_back(phiVmassEff);
      vec_yAxisProjections.push_back(phiVpairPtEff);
      plotSaveNames.push_back("phiVeff");
    }
    // Create "effective efficiency plots"
    // Raw signal/fully corrected signal
    TH1D* effEff_Xaxis = (TH1D*)vec_xAxisProjections.at(2)->Clone("effEffXaxis");
    if(!effEff_Xaxis->Divide(effEff_Xaxis, vec_xAxisProjections.at(4), 1, 1, "B")){
      printError("Xaxis effective efficiency division failed!");
    }
    vec_xAxisProjections.push_back(effEff_Xaxis);
    TH1D* effEff_Yaxis = (TH1D*)vec_yAxisProjections.at(2)->Clone("effEffYaxis");
    if(!effEff_Yaxis->Divide(effEff_Yaxis, vec_yAxisProjections.at(4), 1, 1, "B")){
      printError("Yaxis effective efficiency division failed!");
    }
    vec_yAxisProjections.push_back(effEff_Yaxis);
    plotSaveNames.push_back("effectiveEff");
  }

  //###########################################################
  // ----- Event scaling and rebinning -------
  if(scaleByPhysSel && numEvents != 0){
    numEvents *= (zAxisRange[1] - zAxisRange[0])/100.;
    std::cout << "Scaling spectra by " << numEvents << " events" << std::endl;
    vec_xAxisProjections.at(4)->Scale(1./numEvents);
    vec_yAxisProjections.at(4)->Scale(1./numEvents);
    spectrum->Scale(1./numEvents);
  }

  // Scale projected plots by bin width
  if(scaleByBinWidth){
    for(Int_t i = 0; i < 5; ++i){
      vec_xAxisProjections.at(i)->Scale(1., "width");
      vec_yAxisProjections.at(i)->Scale(1., "width");
    }
  }


  //#########################################################
  //------ Create labels for plots
  TString multSel;
  if(zAxisRange[0] == 0 && zAxisRange[1] == 100){
    multSel = "MB";
  }else{
    multSel = TString::Format("%.0f - %.0f", zAxisRange[0], zAxisRange[1]);
  }

  // Labels for plots.
  // Two sets:
  // 1) plots small/no y axis label ("left")
  // 2) plots with large Y axis labels ("shifted")
  Float_t yPositionTexLabel = !latexTitle.IsNull() ? 0.9 : 0.94;
  TLatex* texTitleLeft = texLabels::getTitle(0.2, yPositionTexLabel, latexTitle, plotting::textSize);
  std::vector<TLatex*> massLabelsLeft   = texLabels::getMassTexLabels(multSel, yAxisRange[0], yAxisRange[1], 0.2, yPositionTexLabel-0.05, plotting::textSize);
  std::vector<TLatex*> pairPtLabelsLeft = texLabels::getPairPtTexLabels(multSel, xAxisRange[0], xAxisRange[1], 0.2, yPositionTexLabel-0.05, plotting::textSize);

  // Labels for the final spectra
  TLatex* texTitleShifted = texLabels::getTitle(0.25, yPositionTexLabel, latexTitle, plotting::textSize);
  std::vector<TLatex*> massLabelsShifted   = texLabels::getMassTexLabels(multSel, yAxisRange[0], yAxisRange[1], 0.25, yPositionTexLabel-0.05, plotting::textSize);
  std::vector<TLatex*> pairPtLabelsShifted = texLabels::getPairPtTexLabels(multSel, xAxisRange[0], xAxisRange[1], 0.25, yPositionTexLabel-0.05, plotting::textSize);

  // Set 2D spectra titles and limit ranges
  {
    TString spectraTitles[] = {"Unlike-sign Dielectron Pairs", "Positive like-sign dielectron pairs",
                      "Negative like-sign dielectron pairs",  "Raw Spectrum"};

    for(Int_t i = 0; i < 4; ++i){
      if(hasEffCorr){
        vec_effInputs.at(i)->GetXaxis()->SetRangeUser(xAxisRange[0]-0.00001, xAxisRange[1]-0.00001);
        vec_effInputs.at(i)->GetYaxis()->SetRangeUser(yAxisRange[0]-0.00001, yAxisRange[1]-0.00001);
      }
      plotting::Set2DTitles(vec_inputSpectra2D.at(i), spectraTitles[i]);
    }
  }

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  //#########################################################
  // Plot 2D spectra
  TCanvas* canvasSpectra = new TCanvas("canvasSpectra", "canvasSpectra");
  canvasSpectra->SetWindowSize(1700, 1100);
  canvasSpectra->SetCanvasSize(1650, 1000);
  canvasSpectra->Divide(3,2);
  for(Int_t i = 0; i < 4; ++i){
    canvasSpectra->cd(i+1);
    gPad->SetLogz();
    gPad->SetRightMargin(0.15);
    if(restrict2DplotRanges){
      vec_inputSpectra2D.at(i)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_inputSpectra2D.at(i)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
    }
    vec_inputSpectra2D.at(i)->Draw("COLZ");
  }
  canvasSpectra->cd(5);
  /* gPad->SetLogz(); */
  gPad->SetRightMargin(0.15);
  spectrum->SetTitle("Corrected Signal");
  spectrum->Draw("COLZ");


  // Define margins for projections plot
  Float_t marginRight = 0.03;
  Float_t marginLeft = 0.2;
  Float_t marginBottom = 0.15;
  Float_t marginTop = 0.01;
  //#########################################################
  //--------- ULS and LS distributions
  plotting::format1Dhist(vec_xAxisProjections.at(0), kBlue);
  plotting::format1Dhist(vec_xAxisProjections.at(1), kRed);
  plotting::format1Dhist(vec_xAxisProjections.at(2), plotting::kGREEN);
  plotting::format1Dhist(vec_yAxisProjections.at(0), kBlue);
  plotting::format1Dhist(vec_yAxisProjections.at(1), kRed);
  plotting::format1Dhist(vec_yAxisProjections.at(2), plotting::kGREEN);

  TLegend* legUSLS = new TLegend(0.53, 0.58, 0.98, 0.69);
  legUSLS->SetFillStyle(0);
  legUSLS->AddEntry(vec_xAxisProjections.at(0), "Unlike-sign pairs",       "lep");
  legUSLS->AddEntry(vec_xAxisProjections.at(1), "Like-sign combinatorial background", "lep");
  legUSLS->AddEntry(vec_xAxisProjections.at(2), "Raw dielectron signal",   "lep");

  TCanvas* canvXaxisInputs = new TCanvas("canvXaxisInputs", "Raw Spectra", plotting::canvWidth, plotting::canvHeight);
  gPad->SetLogy();
  gPad->SetTopMargin(marginTop);
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  vec_xAxisProjections.at(0)->SetMaximum(vec_xAxisProjections.at(0)->GetMaximum()*300);
  vec_xAxisProjections.at(0)->SetTitle("");
  vec_xAxisProjections.at(0)->SetYTitle(plotting::yieldTotalMassLabel);
  vec_xAxisProjections.at(0)->SetXTitle(plotting::massAxisLabel);
  vec_xAxisProjections.at(0)->GetYaxis()->SetTitleOffset(1.7);
  vec_xAxisProjections.at(0)->Draw("E0");
  vec_xAxisProjections.at(1)->Draw("SAME E0");
  vec_xAxisProjections.at(2)->Draw("SAME E0");
  legUSLS->Draw();
  if(!latexTitle.IsNull()){
    texTitleShifted->Draw();
  }
  for(Int_t i = 0; i < massLabelsShifted.size()-1; ++i){
    massLabelsShifted.at(i)->Draw();
  }
  TCanvas* canvYaxisInputs = new TCanvas("canvYaxisInputs", "Raw Spectra", plotting::canvWidth, plotting::canvHeight);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  gPad->SetBottomMargin(marginBottom);
  vec_yAxisProjections.at(0)->SetMaximum(vec_yAxisProjections.at(0)->GetMaximum()*3000);
  vec_yAxisProjections.at(0)->SetTitle("");
  vec_yAxisProjections.at(0)->SetYTitle(plotting::yieldTotalPairPtLabel);
  vec_yAxisProjections.at(0)->SetXTitle(plotting::pairPtAxisLabel);
  vec_yAxisProjections.at(0)->GetYaxis()->SetTitleOffset(1.7);
  vec_yAxisProjections.at(0)->Draw("E0");
  vec_yAxisProjections.at(1)->Draw("SAME E0");
  vec_yAxisProjections.at(2)->Draw("SAME E0");
  legUSLS->Draw();
  if(!latexTitle.IsNull()){
    texTitleShifted->Draw();
  }
  for(Int_t i = 0; i < pairPtLabelsShifted.size()-1; ++i){
    pairPtLabelsShifted.at(i)->Draw();
  }

  //#########################################################
  //-------- Raw signal R factor comparison
  plotting::format1Dhist(vec_xAxisProjections.at(3), kBlue);
  plotting::format1Dhist(vec_yAxisProjections.at(3), kBlue);

  TLegend* legRatio = new TLegend(0.59, 0.43, 0.83, 0.79);
  legRatio->AddEntry(vec_xAxisProjections.at(2), "Full R factor", "lep");
  legRatio->AddEntry(vec_xAxisProjections.at(3), "No R factor", "lep");

  TCanvas* canvXaxisRfacComp = 0x0;
  TCanvas* canvYaxisRfacComp = 0x0;
  if(hasMixedEvents){
    canvXaxisRfacComp = new TCanvas("canvXaxisRfacComp", "R factor comparison");
    canvXaxisRfacComp->SetWindowSize(900, 900);
    canvXaxisRfacComp->SetCanvasSize(850, 850);
    TPad* padMass1 = new TPad("padMass1", "padMass1", 0, 0.3, 1, 1.0);
    padMass1->SetBottomMargin(0);
    padMass1->Draw();
    padMass1->cd();
    gPad->SetLogy();
    vec_xAxisProjections.at(3)->SetMaximum(vec_xAxisProjections.at(3)->GetMaximum()*100);
    vec_xAxisProjections.at(3)->SetTitle("R factor effect");
    vec_xAxisProjections.at(3)->SetYTitle("d#it{N}_{ee}/d#it{m}_{ee} (GeV/#it{c}^{2})^{-1}");
    vec_xAxisProjections.at(3)->Draw("E0");
    vec_xAxisProjections.at(2)->Draw("SAME E0");
    legRatio->Draw();
    for(Int_t i = 0; i < massLabelsLeft.size(); ++i){
      massLabelsLeft.at(i)->Draw();
    }
    TGaxis* axisMass = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    axisMass->SetLabelFont(43);
    axisMass->SetLabelSize(15);
    axisMass->Draw();
    canvXaxisRfacComp->cd();
    TPad* padMass2 = new TPad("padMass2", "padMass2", 0, 0.05, 1, 0.3);
    padMass2->SetTopMargin(0);
    padMass2->SetBottomMargin(0.2);
    padMass2->Draw();
    padMass2->cd();
    TH1D* rFacMassRatio = (TH1D*)vec_xAxisProjections.at(2)->Clone("rFacMassRatio");
    rFacMassRatio->Divide(vec_xAxisProjections.at(3));
    plotting::formatRatioPlot(rFacMassRatio, "with/without", kBlack);
    rFacMassRatio->SetMinimum(0.95);
    rFacMassRatio->SetMaximum(1.05);
    rFacMassRatio->Draw("HIST");

    canvYaxisRfacComp = new TCanvas("canvYaxisRfacComp", "R factor comparison");
    canvYaxisRfacComp->SetWindowSize(900, 900);
    canvYaxisRfacComp->SetCanvasSize(850, 850);
    TPad* padPairPt1 = new TPad("padPairPt1", "padPairPt1", 0, 0.3, 1, 1.0);
    padPairPt1->SetBottomMargin(0);
    padPairPt1->Draw();
    padPairPt1->cd();
    gPad->SetLogy();
    vec_yAxisProjections.at(3)->SetMaximum(vec_yAxisProjections.at(3)->GetMaximum()*100);
    vec_yAxisProjections.at(3)->SetTitle("R factor effect");
    vec_yAxisProjections.at(3)->SetYTitle("d#it{N}_{ee}/d#it{m}_{ee} (GeV/#it{c}^{2})^{-1}");
    vec_yAxisProjections.at(3)->Draw("E0");
    vec_yAxisProjections.at(2)->Draw("SAME E0");
    legRatio->Draw();
    for(Int_t i = 0; i < pairPtLabelsLeft.size(); ++i){
      pairPtLabelsLeft.at(i)->Draw();
    }
    TGaxis* axisPairPt = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    axisPairPt->SetLabelFont(43);
    axisPairPt->SetLabelSize(15);
    axisPairPt->Draw();
    canvYaxisRfacComp->cd();
    TPad* padPairPt2 = new TPad("padPairPt2", "padPairPt2", 0, 0.05, 1, 0.3);
    padPairPt2->SetTopMargin(0);
    padPairPt2->SetBottomMargin(0.2);
    padPairPt2->Draw();
    padPairPt2->cd();
    TH1D* rFacPairPtRatio = (TH1D*)vec_yAxisProjections.at(2)->Clone("rFacPairPtRatio");
    rFacPairPtRatio->Divide(vec_yAxisProjections.at(3));
    plotting::formatRatioPlot(rFacPairPtRatio, "with/without", kBlack);
    rFacPairPtRatio->SetMinimum(0.95);
    rFacPairPtRatio->SetMaximum(1.05);
    rFacPairPtRatio->Draw("HIST");
  }

  //#########################################################
  //-------------  Final Spectra
  plotting::format1Dhist(vec_xAxisProjections.at(4), plotting::kGREEN, plotting::massAxisLabel);
  plotting::format1Dhist(vec_yAxisProjections.at(4), plotting::kGREEN, plotting::pairPtAxisLabel);

  if(scaleByPhysSel){
    vec_xAxisProjections.at(4)->SetYTitle(plotting::yieldMassLabel);
    vec_yAxisProjections.at(4)->SetYTitle(plotting::yieldPairPtLabel);
  }

  TCanvas* canvXaxisSignal = new TCanvas("canvXaxisSignal", "Dielectron Spectrum", plotting::canvWidth, plotting::canvHeight);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  vec_xAxisProjections.at(4)->SetTitle("");
  vec_xAxisProjections.at(4)->GetYaxis()->SetTitleOffset(1.7);
  vec_xAxisProjections.at(4)->SetMaximum(vec_xAxisProjections[4]->GetMaximum()*10);
  vec_xAxisProjections.at(4)->Draw("E0");
  if(!latexTitle.IsNull()){
    texTitleShifted->Draw();
  }
  for(Int_t i = 0; i < massLabelsShifted.size()-1; ++i){
    massLabelsShifted.at(i)->Draw();
  }
  TCanvas* canvYaxisSignal = new TCanvas("canvYaxisSignal", "Dielectron Spectrum", plotting::canvWidth, plotting::canvHeight);
  gPad->SetLogy();
  gPad->SetRightMargin(marginRight);
  gPad->SetLeftMargin(marginLeft);
  gPad->SetBottomMargin(marginBottom);
  vec_yAxisProjections.at(4)->SetTitle("");
  vec_yAxisProjections.at(4)->GetYaxis()->SetTitleOffset(1.7);
  vec_yAxisProjections.at(4)->SetMaximum(vec_yAxisProjections[4]->GetMaximum()*500);
  vec_yAxisProjections.at(4)->Draw("E0");
  if(!latexTitle.IsNull()){
    texTitleShifted->Draw();
  }
  for(Int_t i = 0; i < pairPtLabelsShifted.size()-1; ++i){
    pairPtLabelsShifted.at(i)->Draw();
  }

  //#########################################################
  //------ SB and significance plots
  plotting::format1Dhist(vec_xAxisProjections.at(5), kBlue);
  plotting::format1Dhist(vec_xAxisProjections.at(6), kBlue);//plotting::kGREEN);

  plotting::format1Dhist(vec_yAxisProjections.at(5), kBlue);
  plotting::format1Dhist(vec_yAxisProjections.at(6), plotting::kGREEN);

  TLegend* legStats = new TLegend(0.47, 0.58, 0.92, 0.69);
  legStats->SetFillStyle(0);
  legStats->AddEntry(vec_xAxisProjections.at(5), "S/B", "lep");
  legStats->AddEntry(vec_xAxisProjections.at(6), "Statistical significance", "lep");

  TCanvas* canvStatsXaxis = new TCanvas("canvStatsXaxis", "canvStatsXaxis", plotting::canvWidth, plotting::canvHeight);
  gPad->SetLogy();
  gPad->SetRightMargin(plotting::marginRight);
  gPad->SetBottomMargin(marginBottom);
  vec_xAxisProjections.at(5)->SetXTitle(plotting::massAxisLabel);
  vec_xAxisProjections.at(5)->SetTitle("");
  vec_xAxisProjections.at(5)->SetMinimum(0.001);
  vec_xAxisProjections.at(5)->SetMaximum(100000);
  vec_xAxisProjections.at(5)->Draw("E0");
  vec_xAxisProjections.at(6)->Draw("E0");
  legStats->Draw();
  if(!latexTitle.IsNull()){
    texTitleLeft->Draw();
  }
  for(Int_t i = 0; i < massLabelsShifted.size()-1; ++i){
    massLabelsLeft.at(i)->Draw();
  }

  TCanvas* canvStatsYaxis = new TCanvas("canvStatsYaxis", "canvStatsYaxis");
  canvStatsYaxis->SetWindowSize(900, 900);
  canvStatsYaxis->SetCanvasSize(850, 850);
  gPad->SetRightMargin(marginRight);
  gPad->SetBottomMargin(marginBottom);
  gPad->SetLogy();
  vec_yAxisProjections.at(5)->SetXTitle(plotting::pairPtAxisLabel);
  vec_yAxisProjections.at(5)->SetTitle("");
  vec_yAxisProjections.at(5)->SetMinimum(0.0001);
  vec_yAxisProjections.at(5)->SetMaximum(vec_yAxisProjections.at(6)->GetMaximum()*10000);
  vec_yAxisProjections.at(5)->Draw("E0");
  vec_yAxisProjections.at(6)->Draw("SAME E0");
  legStats->Draw();
  if(!latexTitle.IsNull()){
    texTitleLeft->Draw();
  }
  for(Int_t i = 0; i < pairPtLabelsShifted.size(); ++i){
    pairPtLabelsLeft.at(i)->Draw();
  }


  //#########################################################
  //------- R factor plot
  if(hasMixedEvents){
    plotting::format1Dhist(vec_xAxisProjections.at(7), kBlue);
    plotting::format1Dhist(vec_yAxisProjections.at(7), kBlue);
  }

  TCanvas* canvRfac = 0x0;
  if(hasMixedEvents){
    canvRfac = new TCanvas("canvRfac", "R factor", plotting::canvWidth, plotting::canvHeight);
    canvRfac->Divide(3,1);
    canvRfac->cd(1);
    RforComp->SetTitle("R factor");
    RforComp->SetYTitle(plotting::pairPtAxisLabel);
    RforComp->GetZaxis()->SetRangeUser(0.7,1.3);
    gPad->SetRightMargin(0.15);
    RforComp->Draw("COLZ");
    canvRfac->cd(2);
    vec_xAxisProjections.at(7)->SetTitle("");
    vec_xAxisProjections.at(7)->SetXTitle(plotting::massAxisLabel);
    vec_xAxisProjections.at(7)->SetYTitle("R");
    vec_yAxisProjections.at(7)->GetYaxis()->SetTitleOffset(1.7);
    vec_xAxisProjections.at(7)->SetMaximum(1.1);
    vec_xAxisProjections.at(7)->SetMinimum(0.9);
    vec_xAxisProjections.at(7)->Draw();
    if(!latexTitle.IsNull()){
      texTitleLeft->Draw();
    }
    for(Int_t i = 0; i < massLabelsLeft.size()-1; ++i){
      massLabelsLeft.at(i)->Draw();
    }
    canvRfac->cd(3);
    vec_yAxisProjections.at(7)->SetTitle("R factor");
    vec_yAxisProjections.at(7)->SetXTitle(plotting::pairPtAxisLabel);
    vec_yAxisProjections.at(7)->SetYTitle(plotting::rFactorLabel);
    vec_yAxisProjections.at(7)->GetYaxis()->SetTitleOffset(1.7);
    vec_yAxisProjections.at(7)->SetMaximum(1.1);
    vec_yAxisProjections.at(7)->SetMinimum(0.9);
    vec_yAxisProjections.at(7)->Draw();
    for(Int_t i = 0; i < pairPtLabelsShifted.size(); ++i){
      pairPtLabelsShifted.at(i)->Draw();
    }
  }

  //#########################################################
  //------------ Efficiency plots
  if(hasEffCorr){
    if(!effCorrLF.IsNull()){
      plotting::format1Dhist(vec_xAxisProjections.at(8), kBlue, plotting::massAxisLabel, "#epsilon_{signal}");
      plotting::format1Dhist(vec_yAxisProjections.at(8), kBlue, plotting::pairPtAxisLabel, "#epsilon_{signal}");
    }
    if(!effCorrHF.IsNull()){
      plotting::format1Dhist(vec_xAxisProjections.at(9), kBlack, plotting::massAxisLabel, "#epsilon_{signal}");
      plotting::format1Dhist(vec_yAxisProjections.at(9), kBlack, plotting::pairPtAxisLabel, "#epsilon_{signal}");
    }
    if(usePhiVcorr){
      plotting::format1Dhist(vec_xAxisProjections.at(10), kRed);
      plotting::format1Dhist(vec_yAxisProjections.at(10), kRed);
    }
    plotting::format1Dhist(vec_xAxisProjections.back(), plotting::kORANGE, 4);
    plotting::format1Dhist(vec_yAxisProjections.back(), plotting::kORANGE, 4);
  }

  // Create TLegend for different efficiency inputs
  TLegend* legEff = new TLegend(0, 0.18, 0.7, 0.88);
  if(hasEffCorr){
    if(!effCorrLF.IsNull()){
      legEff->AddEntry(vec_xAxisProjections.at(8), "Track Cut Efficiency (LF)", "lep");
    }
    if(!effCorrHF.IsNull()){
      legEff->AddEntry(vec_xAxisProjections.at(9), "Track Cut Efficiency (HF)", "lep");
    }
    if(usePhiVcorr){
      legEff->AddEntry(vec_xAxisProjections.at(10), "Phi V Cut Efficiency", "lep");
    }
    legEff->AddEntry(vec_xAxisProjections.back(), "Effective Efficiency", "lep");
    legEff->SetTextSize(0.04);
  }

  TCanvas* canvasEff2D = 0x0;
  TCanvas* canvasEff1D = 0x0;
  if(hasEffCorr){
     // 1D efficiency projections
    canvasEff1D = new TCanvas("canvasEff1D", "canvasEff1D");
    canvasEff1D->SetWindowSize(1800, 900);
    canvasEff1D->SetCanvasSize(1750, 850);
    canvasEff1D->Divide(3,1);
    canvasEff1D->cd(1);
    vec_yAxisProjections.at(8)->SetTitle("Efficiencies");
    vec_yAxisProjections.at(8)->SetMinimum(0);
    vec_yAxisProjections.at(8)->SetMaximum(1.1);
    vec_yAxisProjections.at(8)->Draw("E0");
    vec_yAxisProjections.at(9)->Draw("SAME E0");
    if(usePhiVcorr){
      vec_yAxisProjections.at(10)->Draw("SAME E0");
    }
    vec_yAxisProjections.back()->Draw("SAME E0");
    canvasEff1D->cd(2);
    vec_xAxisProjections.at(8)->SetTitle("Efficiencies");
    vec_xAxisProjections.at(8)->SetMinimum(0);
    vec_xAxisProjections.at(8)->SetMaximum(1.1);
    vec_xAxisProjections.at(8)->Draw("E0");
    vec_xAxisProjections.at(9)->Draw("SAME E0");
    if(usePhiVcorr){
      vec_xAxisProjections.at(10)->Draw("SAME E0");
    }
    vec_xAxisProjections.back()->Draw("SAME E0");
    canvasEff1D->cd(3);
    legEff->Draw();

    // 2D efficiency plots
    canvasEff2D = new TCanvas("canvasEff2D", "canvasEff2D", plotting::canvWidth, plotting::canvHeight);
    canvasEff2D->SetWindowSize(1800, 1100);
    canvasEff2D->SetCanvasSize(1750, 1000);
    canvasEff2D->Divide(3,2);
    // Light flavour plots
    canvasEff2D->cd(1);
    gPad->SetLogz();
    plotting::Set2DTitles(vec_effInputs.at(1), "Generated Pairs (LF)");
    vec_effInputs.at(1)->SetMaximum(vec_effInputs.at(1)->GetMaximum()*1);
    if(restrict2DplotRanges){
      vec_effInputs.at(1)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_effInputs.at(1)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
    }
    gPad->SetRightMargin(0.15);
    vec_effInputs.at(1)->Draw(effPlotStyle);
    canvasEff2D->cd(4);
    gPad->SetLogz();
    plotting::Set2DTitles(vec_effInputs.at(2), "Pairs After Cuts (LF)");
    vec_effInputs.at(2)->SetMaximum(vec_effInputs.at(1)->GetMaximum()*1);
    if(restrict2DplotRanges){
      vec_effInputs.at(2)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_effInputs.at(2)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
    }
    gPad->SetRightMargin(0.15);
    vec_effInputs.at(2)->Draw(effPlotStyle);
    // Heavy flavour
    canvasEff2D->cd(2);
    gPad->SetLogz();
    plotting::Set2DTitles(vec_effInputs.at(3), "Generated Pairs (HF)");
    vec_effInputs.at(3)->SetMaximum(vec_effInputs.at(3)->GetMaximum()*1);
    if(restrict2DplotRanges){
      vec_effInputs.at(3)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_effInputs.at(3)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
    }
    gPad->SetRightMargin(0.15);
    vec_effInputs.at(3)->Draw(effPlotStyle);
    canvasEff2D->cd(5);
    gPad->SetLogz();
    plotting::Set2DTitles(vec_effInputs.at(4), "Pairs After Cuts (HF)");
    vec_effInputs.at(4)->SetMaximum(vec_effInputs.at(4)->GetMaximum()*1);
    if(restrict2DplotRanges){
      vec_effInputs.at(4)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_effInputs.at(4)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
    }
    gPad->SetRightMargin(0.15);
    vec_effInputs.at(4)->Draw(effPlotStyle);
    // PhiV correction
    if(usePhiVcorr){
      canvasEff2D->cd(3);
      plotting::Set2DTitles(vec_effInputs.back(), "");
      vec_effInputs.back()->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
      vec_effInputs.back()->GetXaxis()->SetRangeUser(0, 0.5);
    gPad->SetRightMargin(0.15);
      vec_effInputs.back()->Draw(effPlotStyle);
    }
    // Final efficiency
    canvasEff2D->cd(6);
    // Set acceptance hole contents to zero (aesthetic)
    /* for(Int_t i = 0; i < vec_effInputs.at(0)->GetNbinsX(); ++i){ */
    /*   for(Int_t j = 0; j < vec_effInputs.at(0)->GetNbinsY(); ++j){ */
    /*     Double_t xVal = vec_effInputs[0]->GetXaxis()->GetBinCenter(i); */
    /*     Double_t yVal = vec_effInputs[0]->GetYaxis()->GetBinCenter(j); */

    /*     if(TMath::Sqrt(TMath::Power(xVal,2) + TMath::Power(yVal,2)) < 0.4){ */
    /*       vec_effInputs[0]->SetBinContent(i, j, 0); */
    /*     } */

    /*     if( i < 5 && j < 5) std::cout << "Bin: " << i << ", " << j << ": " << vec_effInputs[0]->GetBinContent(i, j) << std::endl; */
    /*   } */
    /* } */
    if(vec_effInputs.at(0)->GetMaximum() > 1){
      printError("WARNING: Efficiency greater than 1!!");
    }
    else{
      vec_effInputs.at(0)->SetMaximum(1.);
    }
    plotting::Set2DTitles(vec_effInputs.at(0), "");
    vec_effInputs.at(0)->SetMaximum(1);
    vec_effInputs.at(0)->SetMinimum(0);
    if(restrict2DplotRanges){
      vec_effInputs.at(0)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_effInputs.at(0)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
    }
    gPad->SetRightMargin(0.15);
    vec_effInputs.at(0)->Draw(effPlotStyle);
  }

    


  TCanvas* canvasWeights = 0x0;
  if(!whichWeights.IsNull()){
    plotting::Set2DTitles(vec_weights.at(0), "Weights (sameMother)");
    plotting::Set2DTitles(vec_weights.at(1), "Weights (Heavy Flavour)");

    if(restrict2DplotRanges){
      vec_weights.at(0)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_weights.at(0)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
      vec_weights.at(1)->GetXaxis()->SetRangeUser(xAxisRange[0]+0.00001, xAxisRange[1]-0.00001);
      vec_weights.at(1)->GetYaxis()->SetRangeUser(yAxisRange[0]+0.00001, yAxisRange[1]-0.00001);
    }
    canvasWeights = new TCanvas("canvasWeights", "canvasWeights");
    canvasWeights->SetCanvasSize(1350, 850);
    canvasWeights->SetWindowSize(1400, 900);
    canvasWeights->Divide(2, 1);
    canvasWeights->cd(1);
    gPad->SetRightMargin(0.15);
    vec_weights.at(0)->Draw("COLZ");
    canvasWeights->cd(2);
    gPad->SetRightMargin(0.15);
    vec_weights.at(1)->Draw("COLZ");
  }
  watch->Stop();
  watch->Print();

  // Alter some save names based on presence of efficiency correction and mixed
  // events
  if(!hasMixedEvents && hasEffCorr){
    plotSaveNames[2] = "rawSignalNoR";
    plotSaveNames[4] = "correctedSignalNoR";
  }
  else if(!hasMixedEvents && hasEffCorr){
    plotSaveNames[2] = "rawSignalNoR";
    plotSaveNames[4] = "rawSignalScaledNoR";
  }
  else if(hasMixedEvents && !hasEffCorr){
    plotSaveNames[4]= "rawSignalScaled";
  }

  TFile* outFile = 0x0;
  if(saveName != ""){
    saveName.Prepend(paths::results_loc + paths::spectra_dir + "corrected/diElec/");
    outFile = TFile::Open(saveName+".root", "RECREATE");
      if(!outFile){
        printError("File could not be opened for writing!");
      }
      else{
        outFile->cd();
        // Save corrected spectra in 2D
        spectrum->Write("spectrum2D_corrected");
        // Save raw 2D inputs
        for(Int_t i = 0; i < vec_inputSpectra2D.size(); ++i){
          vec_inputSpectra2D.at(i)->Write();
        }
        // Save projection plots
        for(Int_t i = 0; i < vec_xAxisProjections.size(); ++i){
          // Skip R factor if not present (only contains NULL pointer)
          if(vec_xAxisProjections.at(i) == 0x0){
            continue;
          }
          vec_xAxisProjections.at(i)->Write(plotSaveNames[i] + "_mass");
          vec_yAxisProjections.at(i)->Write(plotSaveNames[i] + "_pairPt");
        }
        // Save specific canvases
        if(hasMixedEvents){
          canvXaxisRfacComp->Write("canvRfacImpact_mass");
          canvYaxisRfacComp->Write("canvRfacImpact_mass");
        }
        if(hasEffCorr){
          canvasEff2D->Write("canvEff2D");
          canvasEff1D->Write("canvEff1D");
        }
        printInfo("Plots saved to -> " + saveName + ".root");
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

// Get correct bins found using 3D plots
std::vector<Int_t> getBinRanges(Float_t xRange[], Float_t yRange[], Float_t zRange[], const TH3D* inputSpectrum){

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

// Get correct bins found using 2D plots.
std::vector<Int_t> getBinRanges(Float_t xRange[], Float_t yRange[], const TH2D* inputSpectrum){

  std::vector<Int_t> binRanges;

  //##########################################################
  // ------------  Get range restrictions
  binRanges.push_back(inputSpectrum->GetXaxis()->FindBin(xRange[0]+0.00001));
  binRanges.push_back(inputSpectrum->GetXaxis()->FindBin(xRange[1]-0.00001));
  binRanges.push_back(inputSpectrum->GetYaxis()->FindBin(yRange[0]+0.00001));
  binRanges.push_back(inputSpectrum->GetYaxis()->FindBin(yRange[1]-0.00001));

  std::cout << "X axis range: " << xRange[0] << ", " << xRange[1];
  std::cout << ", Bins: " << binRanges[0] << ", " << binRanges[1] << std::endl;
  std::cout << "Y axis range: " << yRange[0] << ", " << yRange[1];
  std::cout << ", Bins: " << binRanges[2] << ", " << binRanges[3] << std::endl;

  return binRanges;
}

TH1D* getRfactorProjections(const std::vector<TH2D*> vec_inputs, const std::vector<Int_t> binRanges, Int_t axis = 0){

  //TODO: remove hard coded numbers
  TH1D* tempPN = 0x0;
  TH1D* tempNN = 0x0;
  TH1D* tempPP = 0x0;
  TH1D* rFacProjection = 0x0;
  if(axis == 0){
    tempPN = (TH1D*)vec_inputs.at(3)->ProjectionX("tempPN", binRanges[2], binRanges[3], "e");
    tempPP = (TH1D*)vec_inputs.at(4)->ProjectionX("tempPP", binRanges[2], binRanges[3], "e");
    tempNN = (TH1D*)vec_inputs.at(5)->ProjectionX("tempNN", binRanges[2], binRanges[3], "e");
  }
  else if(axis == 1){
    tempPN = (TH1D*)vec_inputs.at(3)->ProjectionY("tempPN", binRanges[0], binRanges[1], "e");
    tempPP = (TH1D*)vec_inputs.at(4)->ProjectionY("tempPP", binRanges[0], binRanges[1], "e");
    tempNN = (TH1D*)vec_inputs.at(5)->ProjectionY("tempNN", binRanges[0], binRanges[1], "e");
  }
  if(!tempPN || !tempNN || !tempPP){
    printError("Projections for R factor failed.");
    return 0x0;
  }else{
    rFacProjection = calcDiElecRfactor(tempPN, tempPP, tempNN, (axis == 0 ? 5 : 10));
    return rFacProjection;
  }

}

TH1D* getEffProjections(const std::vector<TH2D*> vec_inputs, const std::vector<Int_t>& binRanges, Int_t axis = 0, Int_t whichProj = 0){

  Int_t proj1, proj2;
  if(whichProj == 0){
    proj1 = 1;
    proj2 = 2;
  }else{
    proj1 = 3;
    proj2 = 4;
  }
  TH1D* effGen = 0x0;
  TH1D* effAcc = 0x0;
  if(axis == 0){
    effGen = (TH1D*)vec_inputs.at(proj1)->ProjectionX("effGen", binRanges[2], binRanges[3], "e");
    effAcc = (TH1D*)vec_inputs.at(proj2)->ProjectionX("effAcc", binRanges[2], binRanges[3], "e");
  }else if(axis == 1){
    effGen = (TH1D*)vec_inputs.at(proj1)->ProjectionY("effGen", binRanges[0], binRanges[1], "e");
    effAcc = (TH1D*)vec_inputs.at(proj2)->ProjectionY("effAcc", binRanges[0], binRanges[1], "e");
  }

  TH1D* effDielec = (TH1D*)effGen->Clone("effDielec");
  effDielec->Divide(effAcc, effGen, 1, 1, "B");

  delete effGen;
  delete effAcc;

  return effDielec;
}
