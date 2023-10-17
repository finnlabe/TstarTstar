#include <iostream>
#include <memory>
#include <string>

// UHH2 stuff
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopPtReweight.h"

// TstarTstar
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/NeuralNetworkModules.h"
#include "UHH2/TstarTstar/include/TstarTstarSignalRegionHists.h"
#include "UHH2/TstarTstar/include/TstarTstarScaleFactors.h"
#include "UHH2/TstarTstar/include/TstarTstarDDTHists.h"
#include "UHH2/TstarTstar/include/ElecTriggerSF.h"

// other
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {

class TstarTstarDNNModule: public AnalysisModule {
public:

  explicit TstarTstarDNNModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  // ###### Modules ######
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<NeuralNetworkIncluder> DNN_Includer;             // main module here, running the DNN inference
  std::unique_ptr<uhh2::AnalysisModule> TstarTstarSpinSwitcher;    // used in case we want to estimate spin 3/2 from spin 1/2 samples
  std::unique_ptr<AnalysisModule> sf_ele_trigger;                  // needed to re-apply electron trigger SFs in case they changed since sel
  std::unique_ptr<AnalysisModule> ttgenprod;
  std::unique_ptr<AnalysisModule> TopPtReweighting;                // re-running this here just to get the event weight :)

  // ##### Histograms #####

  // GEN histograms
  std::unique_ptr<Hists> h_gencheck;
  
  // some "regular" hist collections
  std::unique_ptr<Hists> h_crosscheck, h_hists_SR, h_hists_VR, h_hists_btagCR, h_hists_ttbarCR;
  std::unique_ptr<Hists> h_crosscheck_ele, h_hists_SR_ele, h_hists_VR_ele, h_hists_btagCR_ele, h_hists_ttbarCR_ele;
  std::unique_ptr<Hists> h_crosscheck_mu, h_hists_SR_mu, h_hists_VR_mu, h_hists_btagCR_mu, h_hists_ttbarCR_mu;
  std::unique_ptr<Hists> h_hists_VR_noElecTrigSFs, h_hists_ttbarCR_noElecTrigSFs; // this one is needed to check the effect of the SFs
  std::unique_ptr<Hists> h_AfterDNNcut_06, h_NotDNNcut_06; // used to compare DDT results to simple cut

  // temp check hists
  std::unique_ptr<Hists> h_hists_VR_ele_highHT, h_hists_VR_ele_highST, h_hists_VR_ele_METo300, h_hists_VR_ele_METu300, h_hists_VR_ele_ele300, h_hists_VR_ele_ele500, h_hists_VR_ele_eleu300;

  // DNN output plots
  std::unique_ptr<Hists> h_DNN, h_DNN_DDT, h_DNN_btagCR;

  // now the main "region histograms" which implement only ST, but with all systematic variations
  std::unique_ptr<Hists> h_SignalRegion_total,      h_SignalRegion_mu,        h_SignalRegion_ele;
  std::unique_ptr<Hists> h_ValidationRegion_total,  h_ValidationRegion_mu,    h_ValidationRegion_ele;
  std::unique_ptr<Hists> h_ControlRegion_total,     h_ControlRegion_mu,       h_ControlRegion_ele;
  std::unique_ptr<Hists> h_ttbarControlRegion_total, h_ttbarControlRegion_mu, h_ttbarControlRegion_ele;

  // in case this is data and datadriven BG estimation running, these histograms will be filled for all systematic variations
  std::unique_ptr<Hists> h_SignalRegion_datadriven_FuncUp_total,      h_SignalRegion_datadriven_FuncDown_total;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_FuncUp_mu,         h_SignalRegion_datadriven_FuncDown_mu;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_FuncUp_ele,        h_SignalRegion_datadriven_FuncDown_ele;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_BtagUp_total,      h_SignalRegion_datadriven_BtagDown_total;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_BtagUp_mu,         h_SignalRegion_datadriven_BtagDown_mu;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_BtagUp_ele,        h_SignalRegion_datadriven_BtagDown_ele;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_channelUp_total,   h_SignalRegion_datadriven_channelDown_total;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_channelUp_mu,      h_SignalRegion_datadriven_channelDown_mu;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_channelUp_ele,     h_SignalRegion_datadriven_channelDown_ele;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_yearUp_total,      h_SignalRegion_datadriven_yearDown_total;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_yearUp_mu,         h_SignalRegion_datadriven_yearDown_mu;
  std::unique_ptr<Hists> h_SignalRegion_datadriven_yearUp_ele,        h_SignalRegion_datadriven_yearDown_ele;

  // and the same for the validation region
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_FuncUp_total,     h_ValidationRegion_datadriven_FuncDown_total;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_FuncUp_mu,        h_ValidationRegion_datadriven_FuncDown_mu;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_FuncUp_ele,       h_ValidationRegion_datadriven_FuncDown_ele;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_BtagUp_total,     h_ValidationRegion_datadriven_BtagDown_total;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_BtagUp_mu,        h_ValidationRegion_datadriven_BtagDown_mu;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_BtagUp_ele,       h_ValidationRegion_datadriven_BtagDown_ele;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_channelUp_total,  h_ValidationRegion_datadriven_channelDown_total;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_channelUp_mu,     h_ValidationRegion_datadriven_channelDown_mu;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_channelUp_ele,    h_ValidationRegion_datadriven_channelDown_ele;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_yearUp_total,     h_ValidationRegion_datadriven_yearDown_total;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_yearUp_mu,        h_ValidationRegion_datadriven_yearDown_mu;
  std::unique_ptr<Hists> h_ValidationRegion_datadriven_yearUp_ele,       h_ValidationRegion_datadriven_yearDown_ele;

  // some histograms to test various DDT working points
  std::unique_ptr<uhh2::TstarTstarDDTHists> h_DDTtestHists;  

  // ###### Handles ######
  uhh2::Event::Handle<double> h_evt_weight;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;

  uhh2::Event::Handle<bool> h_trigger_decision;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<bool> h_flag_muonevent;
  uhh2::Event::Handle<bool> h_is_btagevent;

  uhh2::Event::Handle<double> h_DNN_output;
  uhh2::Event::Handle<double> h_ST_AK4;
  uhh2::Event::Handle<double> h_ST_HOTVR;
  uhh2::Event::Handle<double> h_DDT_score;

  uhh2::Event::Handle<TString> h_region;

  uhh2::Event::Handle<float> h_weight_sfelec_triggerNominal;


  // ###### other parameters ######
  bool debug = false;
  bool is_MC;
  bool is_datadriven_BG_run = false;
  TString year;

  // background estimation functions
  TGraphAsymmErrors* bgest_purity;
  TF1 * backgroundEstimationFunctionNominal_SR,   * backgroundEstimationFunctionNominal_VR;
  
  TF1 * backgroundEstimationFunctionFuncUp_SR,    * backgroundEstimationFunctionFuncUp_VR;
  TF1 * backgroundEstimationFunctionFuncDown_SR,  * backgroundEstimationFunctionFuncDown_VR;

  TF1 * backgroundEstimationFunctionBtagUp_SR,    * backgroundEstimationFunctionBtagUp_VR;
  TF1 * backgroundEstimationFunctionBtagDown_SR,  * backgroundEstimationFunctionBtagDown_VR;

  TF1 * backgroundEstimationFunctionChannelUp_SR,    * backgroundEstimationFunctionChannelUp_VR;
  TF1 * backgroundEstimationFunctionChannelDown_SR,  * backgroundEstimationFunctionChannelDown_VR;

  TF1 * backgroundEstimationFunctionYearUp_SR,    * backgroundEstimationFunctionYearUp_VR;
  TF1 * backgroundEstimationFunctionYearDown_SR,  * backgroundEstimationFunctionYearDown_VR;

  // DDT function(s)
  std::vector<TF1*> DDTFunctions;
  TF1* BestDDTFunction;

};


TstarTstarDNNModule::TstarTstarDNNModule(Context & ctx){

  // 0. environment setup

  // setting debug from xml file
  if(ctx.get("debug", "<not set>") == "true") debug = true;

  // ###### 0. Setting variables ######
  // MC or real data
  is_MC = ctx.get("dataset_type") == "MC";

  // year of samples
  year = ctx.get("year", "<not set>");
  if(year == "<not set>"){
    if(ctx.get("dataset_version").find("2016") != std::string::npos) year = "2016";
    else if(ctx.get("dataset_version").find("2017") != std::string::npos) year = "2017";
    else if(ctx.get("dataset_version").find("2018") != std::string::npos) year = "2018";
    else if(ctx.get("dataset_version").find("UL16preVFP") != std::string::npos) year = "UL16preVFP";
    else if(ctx.get("dataset_version").find("UL16postVFP") != std::string::npos) year = "UL16postVFP";
    else if(ctx.get("dataset_version").find("UL17") != std::string::npos) year = "UL17";
    else if(ctx.get("dataset_version").find("UL18") != std::string::npos) year = "UL18";
    else throw "No year found in dataset name!";
  }
  if(debug) cout << "Year is " << year << "." << endl;

  // get which running mode to use for data
  if(!is_MC) is_datadriven_BG_run = ctx.get("use_data_for", "regular") == "background_extrapolation"; 

  // 1. set up modules 
  reco_primlep.reset(new PrimaryLepton(ctx));
  DNN_Includer.reset(new NeuralNetworkIncluder(ctx, false));
  TstarTstarSpinSwitcher.reset(new TstarTstarSpinScale(ctx, "/nfs/dust/cms/user/flabe/TstarTstar/Jupyter/reweight_factors/"));
  sf_ele_trigger.reset( new uhh2::ElecTriggerSF(ctx, "central", "eta_ptbins", year) );

  if(is_MC) ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));

  TopPtReweighting.reset( new TopPtReweight(ctx, 0.0615, -0.0005, "ttbargen", "weight_ttbar", false) ); // just rerunning to get the weight, not applying!


  // 2. set up histograms:

  // GEN
  h_gencheck.reset(new TstarTstarAllGenHists(ctx, "gencheck")); // was previously named topcheck

  // full histograms
  h_crosscheck.reset(new TstarTstarHists(ctx, "crosscheck"));
  h_hists_SR.reset(new TstarTstarHists(ctx, "hists_SR"));
  h_hists_VR.reset(new TstarTstarHists(ctx, "hists_VR"));
  h_hists_btagCR.reset(new TstarTstarHists(ctx, "hists_btagCR"));
  h_hists_ttbarCR.reset(new TstarTstarHists(ctx, "hists_ttbarCR"));

  h_crosscheck_ele.reset(new TstarTstarHists(ctx, "crosscheck_ele"));
  h_hists_SR_ele.reset(new TstarTstarHists(ctx, "hists_SR_ele"));
  h_hists_VR_ele.reset(new TstarTstarHists(ctx, "hists_VR_ele"));
  h_hists_btagCR_ele.reset(new TstarTstarHists(ctx, "hists_btagCR_ele"));
  h_hists_ttbarCR_ele.reset(new TstarTstarHists(ctx, "hists_ttbarCR_ele"));

  h_crosscheck_mu.reset(new TstarTstarHists(ctx, "crosscheck_mu"));
  h_hists_SR_mu.reset(new TstarTstarHists(ctx, "hists_SR_mu"));
  h_hists_VR_mu.reset(new TstarTstarHists(ctx, "hists_VR_mu"));
  h_hists_btagCR_mu.reset(new TstarTstarHists(ctx, "hists_btagCR_mu"));
  h_hists_ttbarCR_mu.reset(new TstarTstarHists(ctx, "hists_ttbarCR_mu"));

  h_hists_VR_ele_highHT.reset(new TstarTstarHists(ctx, "hists_VR_ele_highHT"));
  h_hists_VR_ele_highST.reset(new TstarTstarHists(ctx, "hists_VR_ele_highST"));
  h_hists_VR_ele_METo300.reset(new TstarTstarHists(ctx, "hists_VR_ele_METo300"));
  h_hists_VR_ele_METu300.reset(new TstarTstarHists(ctx, "hists_VR_ele_METu300"));
  h_hists_VR_ele_ele300.reset(new TstarTstarHists(ctx, "hists_VR_ele_ele300"));
  h_hists_VR_ele_eleu300.reset(new TstarTstarHists(ctx, "hists_VR_ele_eleu300"));
  h_hists_VR_ele_ele500.reset(new TstarTstarHists(ctx, "hists_VR_ele_ele500"));

  h_hists_VR_noElecTrigSFs.reset(new TstarTstarHists(ctx, "hists_VR_noElecTrigSFs"));
  h_hists_ttbarCR_noElecTrigSFs.reset(new TstarTstarHists(ctx, "hists_ttbarCR_noElecTrigSFs"));
  h_AfterDNNcut_06.reset(new TstarTstarHists(ctx, "AfterDNNcut_06"));
  h_NotDNNcut_06.reset(new TstarTstarHists(ctx, "NotDNNcut_06"));

  // DNN output plots
  h_DNN.reset(new TstarTstarDNNHists(ctx, "DNN"));
  h_DNN_DDT.reset(new TstarTstarDNNHists(ctx, "DNN_DDT"));
  h_DNN_btagCR.reset(new TstarTstarDNNHists(ctx, "DNN_btagCR"));

  // region plots
  // the TstarTstarSignalRegionHists has all systematic variations included
  h_SignalRegion_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_total"));
  h_SignalRegion_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_mu"));
  h_SignalRegion_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_ele"));

  h_ValidationRegion_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_total"));
  h_ValidationRegion_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_mu"));
  h_ValidationRegion_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_ele"));

  h_ControlRegion_total.reset(new TstarTstarSignalRegionHists(ctx, "ControlRegion_total"));
  h_ControlRegion_mu.reset(new TstarTstarSignalRegionHists(ctx, "ControlRegion_mu"));
  h_ControlRegion_ele.reset(new TstarTstarSignalRegionHists(ctx, "ControlRegion_ele"));

  h_ttbarControlRegion_total.reset(new TstarTstarSignalRegionHists(ctx, "ttbarControlRegion_total"));
  h_ttbarControlRegion_mu.reset(new TstarTstarSignalRegionHists(ctx, "ttbarControlRegion_mu"));
  h_ttbarControlRegion_ele.reset(new TstarTstarSignalRegionHists(ctx, "ttbarControlRegion_ele"));

  h_SignalRegion_datadriven_FuncUp_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncUp_total"));
  h_SignalRegion_datadriven_FuncDown_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncDown_total"));
  h_ValidationRegion_datadriven_FuncUp_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncUp_total")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_FuncDown_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncDown_total")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_FuncUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncUp_mu"));
  h_SignalRegion_datadriven_FuncDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncDown_mu"));
  h_ValidationRegion_datadriven_FuncUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncUp_mu")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_FuncDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncDown_mu")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_FuncUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncUp_ele"));
  h_SignalRegion_datadriven_FuncDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_FuncDown_ele"));
  h_ValidationRegion_datadriven_FuncUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncUp_ele")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_FuncDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_FuncDown_ele")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_BtagUp_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagUp_total"));
  h_SignalRegion_datadriven_BtagDown_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagDown_total"));
  h_ValidationRegion_datadriven_BtagUp_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagUp_total")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_BtagDown_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagDown_total")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_BtagUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagUp_mu"));
  h_SignalRegion_datadriven_BtagDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagDown_mu"));
  h_ValidationRegion_datadriven_BtagUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagUp_mu")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_BtagDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagDown_mu")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_BtagUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagUp_ele"));
  h_SignalRegion_datadriven_BtagDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_BtagDown_ele"));
  h_ValidationRegion_datadriven_BtagUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagUp_ele")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_BtagDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_BtagDown_ele")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_channelUp_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_channelUp_total"));
  h_SignalRegion_datadriven_channelDown_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_channelDown_total"));
  h_ValidationRegion_datadriven_channelUp_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_channelUp_total")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_channelDown_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_channelDown_total")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_channelUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_channelUp_mu"));
  h_SignalRegion_datadriven_channelDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_channelDown_mu"));
  h_ValidationRegion_datadriven_channelUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_channelUp_mu")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_channelDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_channelDown_mu")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_channelUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_channelUp_ele"));
  h_SignalRegion_datadriven_channelDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_channelDown_ele"));
  h_ValidationRegion_datadriven_channelUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_channelUp_ele")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_channelDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_channelDown_ele")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_yearUp_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_yearUp_total"));
  h_SignalRegion_datadriven_yearDown_total.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_yearDown_total"));
  h_ValidationRegion_datadriven_yearUp_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_yearUp_total")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_yearDown_total.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_yearDown_total")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_yearUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_yearUp_mu"));
  h_SignalRegion_datadriven_yearDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_yearDown_mu"));
  h_ValidationRegion_datadriven_yearUp_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_yearUp_mu")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_yearDown_mu.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_yearDown_mu")); // using SR hists here, although it is not

  h_SignalRegion_datadriven_yearUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_yearUp_ele"));
  h_SignalRegion_datadriven_yearDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "SignalRegion_datadriven_yearDown_ele"));
  h_ValidationRegion_datadriven_yearUp_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_yearUp_ele")); // using SR hists here, although it is not
  h_ValidationRegion_datadriven_yearDown_ele.reset(new TstarTstarSignalRegionHists(ctx, "ValidationRegion_datadriven_yearDown_ele")); // using SR hists here, although it is not

  // DDT check histograms, at various cut values
  std::vector<TString> DDT_points_to_check = {"0p3"}; // also used below to load the files
  h_DDTtestHists.reset( new TstarTstarDDTHists(ctx, "DDTHists", DDT_points_to_check) ) ; 

  // 3. init handles
  ctx.undeclare_all_event_output(); // we'll not store the tree, as its not needed  (this is the last UHH2 step) -> saving space!
  h_evt_weight = ctx.get_handle<double>("evt_weight");

  h_trigger_decision = ctx.get_handle<bool>("trigger_decision");
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");
  h_flag_muonevent = ctx.get_handle<bool>("is_muevt");
  h_is_btagevent = ctx.get_handle<bool>("is_btagevent");

  h_DNN_output = ctx.get_handle<double>("DNN_output");
  h_ST_AK4 = ctx.get_handle<double>("ST_AK4");
  h_ST_HOTVR = ctx.get_handle<double>("ST_HOTVR");
  h_DDT_score = ctx.declare_event_output<double>("DDT_score");
  
  h_region = ctx.declare_event_output<TString>("region");

  h_weight_sfelec_triggerNominal = ctx.get_handle<float>("weight_sfelec_trigger");


  // 4. load files

  // datadriven estimation functions only if needed
  if(is_datadriven_BG_run) {
    TString path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/bgest/";

    TString filename_nominal_SR = "alphaFunction_HOTVR__SR_total.root";
    TString filename_nominal_VR = "alphaFunction_HOTVR___VR_total.root"; // TODO hard-coded UL17 dummy

    TString filename_btagUp_SR = "alphaFunction_HOTVR__SR_total_btagging_totalUp.root";
    TString filename_btagDown_SR = "alphaFunction_HOTVR__SR_total_btagging_totalDown.root";
    TString filename_btagUp_VR = "alphaFunction_HOTVR__VR_total_btagging_totalUp.root";
    TString filename_btagDown_VR = "alphaFunction_HOTVR__VR_total_btagging_totalDown.root";

    // using ele and mu channel as up and down variations
    // could also think of just doing these split entirely -> to check what the effect is!
    TString filename_channelUp_SR = "alphaFunction_HOTVR__SR_mu.root";
    TString filename_channelDown_SR = "alphaFunction_HOTVR__SR_ele.root";
    TString filename_channelUp_VR = "alphaFunction_HOTVR__VR_mu.root";
    TString filename_channelDown_VR = "alphaFunction_HOTVR__VR_ele.root";

    // using year extremes as up and down variations
    // these were selected manually and need to be replaced once re-running the analysis!
    // UL16preVFP and UL17 are found to be the extremes
    TString filename_yearUp_SR = "alphaFunction_HOTVR_UL16preVFP_SR_total.root";
    TString filename_yearDown_SR = "alphaFunction_HOTVR_UL17_SR_total.root";
    TString filename_yearUp_VR = "alphaFunction_HOTVR_UL16preVFP_VR_total.root";
    TString filename_yearDown_VR = "alphaFunction_HOTVR_UL17_VR_total.root";


    // main file for signal region
    TFile *file_nominal_SR = new TFile(path+filename_nominal_SR);
    backgroundEstimationFunctionNominal_SR = (TF1*)file_nominal_SR->Get("fit_mean");
    backgroundEstimationFunctionFuncUp_SR = (TF1*)file_nominal_SR->Get("fit1");
    backgroundEstimationFunctionFuncDown_SR = (TF1*)file_nominal_SR->Get("fit2");

    // main file for validation region
    TFile *file_nominal_VR = new TFile(path+filename_nominal_VR);
    backgroundEstimationFunctionNominal_VR = (TF1*)file_nominal_VR->Get("fit_mean");
    backgroundEstimationFunctionFuncUp_VR = (TF1*)file_nominal_VR->Get("fit1");
    backgroundEstimationFunctionFuncDown_VR = (TF1*)file_nominal_VR->Get("fit2");


    // now getting the btag variations
    TFile *file_btagUp_SR = new TFile(path+filename_btagUp_SR);
    backgroundEstimationFunctionBtagUp_SR = (TF1*)file_btagUp_SR->Get("fit_mean");
    TFile *file_btagDown_SR = new TFile(path+filename_btagDown_SR);
    backgroundEstimationFunctionBtagDown_SR = (TF1*)file_btagDown_SR->Get("fit_mean");
    TFile *file_btagUp_VR = new TFile(path+filename_btagUp_VR);
    backgroundEstimationFunctionBtagUp_VR = (TF1*)file_btagUp_VR->Get("fit_mean");
    TFile *file_btagDown_VR = new TFile(path+filename_btagDown_VR);
    backgroundEstimationFunctionBtagDown_VR = (TF1*)file_btagDown_VR->Get("fit_mean");


    // channel variations
    TFile *file_channelUp_SR = new TFile(path+filename_channelUp_SR);
    backgroundEstimationFunctionChannelUp_SR = (TF1*)file_channelUp_SR->Get("fit_mean");
    TFile *file_channelDown_SR = new TFile(path+filename_channelDown_SR);
    backgroundEstimationFunctionChannelDown_SR = (TF1*)file_channelDown_SR->Get("fit_mean");
    TFile *file_channelUp_VR = new TFile(path+filename_channelUp_VR);
    backgroundEstimationFunctionChannelUp_VR = (TF1*)file_channelUp_VR->Get("fit_mean");
    TFile *file_channelDown_VR = new TFile(path+filename_channelDown_VR);
    backgroundEstimationFunctionChannelDown_VR = (TF1*)file_channelDown_VR->Get("fit_mean");


    // year variations
    TFile *file_yearUp_SR = new TFile(path+filename_yearUp_SR);
    backgroundEstimationFunctionYearUp_SR = (TF1*)file_yearUp_SR->Get("fit_mean");
    TFile *file_yearDown_SR = new TFile(path+filename_yearDown_SR);
    backgroundEstimationFunctionYearDown_SR = (TF1*)file_yearDown_SR->Get("fit_mean");
    TFile *file_yearUp_VR = new TFile(path+filename_yearUp_VR);
    backgroundEstimationFunctionYearUp_VR = (TF1*)file_yearUp_VR->Get("fit_mean");
    TFile *file_yearDown_VR = new TFile(path+filename_yearDown_VR);
    backgroundEstimationFunctionYearDown_VR = (TF1*)file_yearDown_VR->Get("fit_mean");

    // finally, the purity
    TString background_estimation_purity_filepath = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/bgest/purity_HOTVR__total.root"; // TODO hard-coded dummy
    TFile *f = new TFile(background_estimation_purity_filepath);
    if(!f) throw std::runtime_error("ERROR: cant open background estimation purity at " + ctx.get("background_estimation_purity_file"));
    bgest_purity = (TGraphAsymmErrors*)f->Get("purity");
    if(!bgest_purity) throw std::runtime_error("ERROR: File at " + ctx.get("background_estimation_purity_file") + " does not contain a purity graph!");
  }

  // loading the DDT fit function(s)
  TString path = "/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/";
  TString bestFunction = "0p3";
  TString filename_base = "DDTfunc_";
  
  // getting the best function
  TFile *f = new TFile(path+filename_base+bestFunction+".root");
  BestDDTFunction = (TF1*)f->Get("mean");

  // now lets load all the others
  for (auto point : DDT_points_to_check) {
    TFile *f = new TFile(path+filename_base+point+".root");
    auto DDTFunction = (TF1*)f->Get("mean");
    DDTFunctions.push_back(DDTFunction);
  }

}


bool TstarTstarDNNModule::process(Event & event) {

  // at the moment, throw out all events that do not pass the trigger here
  // later, this info can be used to calculate trigger efficiencies in the different regions, if needed
  if (!event.get(h_trigger_decision)) return false;

  // debug message
  if(debug){cout << endl << "TstarTstarDNNModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  // fixing poitential empty handle
  event.set(h_region, "not set");

  // reapply weights
  event.weight = event.get(h_evt_weight);
  if(debug) cout << "weights applied." << endl;

  // reapply electron trigger SFs (as they changed since running preselection)
  if(!event.get(h_flag_muonevent)) {
    if( event.get(h_weight_sfelec_triggerNominal) != 0) event.weight /= event.get(h_weight_sfelec_triggerNominal);
    sf_ele_trigger->process(event);
  }

  if(is_MC) ttgenprod->process(event);
  TopPtReweighting->process(event); // will set the weight 

  // apply spin reweighting
  if(!(TstarTstarSpinSwitcher->process(event))) return false; // this will reweight, but only if set so in config file
  event.set(h_evt_weight, event.weight); // we'll use this later, so need to update

  // set primary lepton
  reco_primlep->process(event);

  // a cross-check whether the btagging flag is correct
  // can be removed if its not failing in the next running iteration
  BTag bJetID = BTag(BTag::algo::DEEPJET, BTag::wp::WP_MEDIUM);
  bool pass_btagcut = false;
  for (const auto & jet: *event.jets){
    if(bJetID(jet, event)) pass_btagcut = true;
  }
  assert(pass_btagcut == event.get(h_is_btagevent));

  // do additional ST cut
  if(event.get(h_ST_HOTVR) < 600) return false;

  // filling crosscheck histograms after all initial steps are done
  if(event.get(h_is_btagevent)) {
    h_crosscheck->fill(event);

    if(event.get(h_flag_muonevent)) {
      h_crosscheck_mu->fill(event);
    } else {
      h_crosscheck_ele->fill(event);
    }

  }


  // ################
  // ### DNN Part ###
  // ################

  // including trained DNN model
  // this will set the h_DNN_output handle
  if(debug) cout << "Include DNN model" << endl;
  DNN_Includer->process(event);
  if(debug) cout << "Done DNN include" << endl;

  // hists
  if(event.get(h_is_btagevent)) {
    if(debug) cout << "Start filling hists" << endl;

    h_DNN->fill(event);
    h_gencheck->fill(event);

    if(event.get(h_DNN_output) > 0.6) h_AfterDNNcut_06->fill(event);
    else h_NotDNNcut_06->fill(event);

  }


  // ################
  // ### DDT Part ###
  // ################

  if(debug) cout << "Start DDT part" << endl;

  // calculating and applying the DDT shift for the best function
  double DDTshift = BestDDTFunction->Eval(event.get(h_ST_HOTVR));
  double DDTscore = event.get(h_DNN_output) - ( 1 - DDTshift );
  event.set(h_DDT_score, DDTscore);
  
  // lets store all the other tagger outputs in a vector
  std::vector<double> newTaggerResults;
  for (auto function : DDTFunctions) {
    double value = event.get(h_DNN_output) - ( 1 - function->Eval(event.get(h_ST_HOTVR)) );
    newTaggerResults.push_back(value);
  }

  if(event.get(h_is_btagevent)) h_DNN_DDT->fill(event);


  // ##################################
  // ### background estimation part ###
  // ##################################

  if(debug) cout << "Start BG estimation part" << endl;

  // initializing some values (set to one in case background estimation is not done)
  double transfer_weight_nominal_SR = 1;
  double transfer_weight_nominal_VR = 1;
  double transfer_weight_funcUp_SR = 1;
  double transfer_weight_funcUp_VR = 1;
  double transfer_weight_funcDown_SR = 1;
  double transfer_weight_funcDown_VR = 1;
  double transfer_weight_btagUp_SR = 1;
  double transfer_weight_btagUp_VR = 1;
  double transfer_weight_btagDown_SR = 1;
  double transfer_weight_btagDown_VR = 1;
  double transfer_weight_channelUp_SR = 1;
  double transfer_weight_channelUp_VR = 1;
  double transfer_weight_channelDown_SR = 1;
  double transfer_weight_channelDown_VR = 1;
  double transfer_weight_yearUp_SR = 1;
  double transfer_weight_yearUp_VR = 1;
  double transfer_weight_yearDown_SR = 1;
  double transfer_weight_yearDown_VR = 1;

  // storing the event weight to be able to reset later
  double weight_for_resetting = event.weight; 
  double purity_value = 1;
  if(is_datadriven_BG_run && !pass_btagcut) {

    if(debug) cout << "Doing datadriven BG estimation" << endl;

    transfer_weight_nominal_SR = backgroundEstimationFunctionNominal_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_nominal_VR = backgroundEstimationFunctionNominal_VR->Eval(event.get(h_ST_HOTVR));

    transfer_weight_funcUp_SR = backgroundEstimationFunctionFuncUp_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_funcUp_VR = backgroundEstimationFunctionFuncUp_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_funcDown_SR = backgroundEstimationFunctionFuncDown_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_funcDown_VR = backgroundEstimationFunctionFuncDown_VR->Eval(event.get(h_ST_HOTVR));

    transfer_weight_btagUp_SR = backgroundEstimationFunctionBtagUp_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_btagUp_VR = backgroundEstimationFunctionBtagUp_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_btagDown_SR = backgroundEstimationFunctionBtagDown_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_btagDown_VR = backgroundEstimationFunctionBtagDown_VR->Eval(event.get(h_ST_HOTVR));

    transfer_weight_channelUp_SR = backgroundEstimationFunctionChannelUp_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_channelUp_VR = backgroundEstimationFunctionChannelUp_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_channelDown_SR = backgroundEstimationFunctionChannelDown_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_channelDown_VR = backgroundEstimationFunctionChannelDown_VR->Eval(event.get(h_ST_HOTVR));

    transfer_weight_yearUp_SR = backgroundEstimationFunctionYearUp_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_yearUp_VR = backgroundEstimationFunctionYearUp_VR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_yearDown_SR = backgroundEstimationFunctionYearDown_SR->Eval(event.get(h_ST_HOTVR));
    transfer_weight_yearDown_VR = backgroundEstimationFunctionYearDown_VR->Eval(event.get(h_ST_HOTVR));

    if(debug) cout << "Gotten all values" << endl;
    
    // fetching the purity (which needs to be account for)
    purity_value = bgest_purity->Eval(event.get(h_ST_HOTVR));
    if(debug) cout << "purity: " << purity_value << endl;

    // fill all the up and down-variations of the signal region and validation region plots
    // the nominal case is handled below to put it into the same place as all other main signal region histograms
    // this will make things easier for combine later

    // function choice
    event.weight *= transfer_weight_funcUp_SR * purity_value;
    h_SignalRegion_datadriven_FuncUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_FuncUp_mu->fill(event);
    else h_SignalRegion_datadriven_FuncUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_funcDown_SR * purity_value;
    h_SignalRegion_datadriven_FuncDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_FuncDown_mu->fill(event);
    else h_SignalRegion_datadriven_FuncDown_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_funcUp_VR * purity_value;
    h_ValidationRegion_datadriven_FuncUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_FuncUp_mu->fill(event);
    else h_ValidationRegion_datadriven_FuncUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_funcDown_VR * purity_value;
    h_ValidationRegion_datadriven_FuncDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_FuncDown_mu->fill(event);
    else h_ValidationRegion_datadriven_FuncDown_ele->fill(event);
    event.weight = weight_for_resetting;

    // btag variation
    event.weight *= transfer_weight_btagUp_SR * purity_value;
    h_SignalRegion_datadriven_BtagUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_BtagUp_mu->fill(event);
    else h_SignalRegion_datadriven_BtagUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_btagDown_SR * purity_value;
    h_SignalRegion_datadriven_BtagDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_BtagDown_mu->fill(event);
    else h_SignalRegion_datadriven_BtagDown_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_btagUp_VR * purity_value;
    h_ValidationRegion_datadriven_BtagUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_BtagUp_mu->fill(event);
    else h_ValidationRegion_datadriven_BtagUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_btagDown_VR * purity_value;
    h_ValidationRegion_datadriven_BtagDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_BtagDown_mu->fill(event);
    else h_ValidationRegion_datadriven_BtagDown_ele->fill(event);
    event.weight = weight_for_resetting;

    // channel variation
    event.weight *= transfer_weight_channelUp_SR * purity_value;
    h_SignalRegion_datadriven_channelUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_channelUp_mu->fill(event);
    else h_SignalRegion_datadriven_channelUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_channelDown_SR * purity_value;
    h_SignalRegion_datadriven_channelDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_channelDown_mu->fill(event);
    else h_SignalRegion_datadriven_channelDown_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_channelUp_VR * purity_value;
    h_ValidationRegion_datadriven_channelUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_channelUp_mu->fill(event);
    else h_ValidationRegion_datadriven_channelUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_channelDown_VR * purity_value;
    h_ValidationRegion_datadriven_channelDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_channelDown_mu->fill(event);
    else h_ValidationRegion_datadriven_channelDown_ele->fill(event);
    event.weight = weight_for_resetting;

    // year variation
    event.weight *= transfer_weight_yearUp_SR * purity_value;
    h_SignalRegion_datadriven_yearUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_yearUp_mu->fill(event);
    else h_SignalRegion_datadriven_yearUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_yearDown_SR * purity_value;
    h_SignalRegion_datadriven_yearDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_SignalRegion_datadriven_yearDown_mu->fill(event);
    else h_SignalRegion_datadriven_yearDown_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_yearUp_VR * purity_value;
    h_ValidationRegion_datadriven_yearUp_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_yearUp_mu->fill(event);
    else h_ValidationRegion_datadriven_yearUp_ele->fill(event);
    event.weight = weight_for_resetting;

    event.weight *= transfer_weight_yearDown_VR * purity_value;
    h_ValidationRegion_datadriven_yearDown_total->fill(event);
    if(event.get(h_flag_muonevent)) h_ValidationRegion_datadriven_yearDown_mu->fill(event);
    else h_ValidationRegion_datadriven_yearDown_ele->fill(event);
    event.weight = weight_for_resetting;

  } // else fill them normally here -> to double check, at least for VR!

  
  // ##############################
  // ### histogram filling part ###
  // ##############################

  if(debug) cout << "Start hist filling part" << endl;

  // note the inversion of pass_btagcut in case we are datadriven, as here we *only* want to have the events not passing the btagging requirement
  if( (!is_datadriven_BG_run && pass_btagcut) || (is_datadriven_BG_run && !pass_btagcut) ) {

    // these bools decide where we'll fill
    // just simplicity to not duplicate these conditions
    bool fillSR = false;
    bool fillVR = false;

    if(is_datadriven_BG_run) {
      // datadriven events are filled into *both* region - but the weight factor will be different!
      fillSR = true;
      fillVR = true;
      event.set(h_region, "datadrivenEst");
      if(debug) cout << "Region: datadriven est." << endl;
    } else if (!is_datadriven_BG_run) {
      if(event.get(h_DDT_score) > 0) {
        fillSR = true;
        fillVR = false;
        event.set(h_region, "SR");
        if(debug) cout << "Region: Signal Region" << endl;
      } else {
        fillSR = false;
        fillVR = true;
        event.set(h_region, "VR");
        if(debug) cout << "Region: Valdation Region" << endl;
      }
    }
   
    if(debug) cout << "Region filling" << endl;
    // now the actual filling
    if(fillSR) {
      if (is_datadriven_BG_run) event.weight *= transfer_weight_nominal_SR * purity_value; // here we change the weight if we are datadriven BG
      if(is_MC /* blinding */ || is_datadriven_BG_run) {
        h_hists_SR->fill(event);
        h_SignalRegion_total->fill(event);
        if(event.get(h_flag_muonevent)) {
          h_hists_SR_mu->fill(event);
          h_SignalRegion_mu->fill(event);
        }
        else {
          h_hists_SR_ele->fill(event);
          h_SignalRegion_ele->fill(event);
        }
      }
      if (is_datadriven_BG_run) event.weight = weight_for_resetting; // resetting the weight
    }
    if(fillVR) {
      if (is_datadriven_BG_run) event.weight *= transfer_weight_nominal_VR * purity_value; // here we change the weight if we are datadriven BG
      h_hists_VR->fill(event);
      h_ValidationRegion_total->fill(event);
      if(event.get(h_flag_muonevent)) {
        h_hists_VR_mu->fill(event);
        h_ValidationRegion_mu->fill(event);
      }
      else {

        h_hists_VR_ele->fill(event);
        h_ValidationRegion_ele->fill(event);

        // quickly calculate HT and ST to check region
        double ht = 0.;
        for(const auto & topjet : *event.topjets) {
          ht += topjet.pt();
        }

        if (ht > 3000) h_hists_VR_ele_highHT->fill(event);
        if (event.get(h_ST_HOTVR) > 3000) h_hists_VR_ele_highST->fill(event);
        if (event.met->pt() > 300) h_hists_VR_ele_METo300->fill(event);
        else h_hists_VR_ele_METu300->fill(event);
        if (event.electrons->at(0).pt() > 300) h_hists_VR_ele_ele300->fill(event);
        else h_hists_VR_ele_eleu300->fill(event);
        if (event.electrons->at(0).pt() > 500) h_hists_VR_ele_ele500->fill(event);

        // special case here: plot same, but without electron trigger SF
        double reset_weight_etrigger = event.weight;
        if( event.get(h_weight_sfelec_triggerNominal) != 0) event.weight /= event.get(h_weight_sfelec_triggerNominal);
        h_hists_VR_noElecTrigSFs->fill(event);
        event.weight = reset_weight_etrigger;
          
      }

      // a subsest of this region will be used as a ttbar CR
      // three different variants of this are defined here

      if(debug) cout << "Filling ttbar CRs" << endl;

      // HOTVR top tag definition
      TopJetId topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
      bool passHOTVRtoptag = false;
      for (const auto & jet: *event.topjets){
        if(topjetID(jet, event)) passHOTVRtoptag = true;
      }

      if(passHOTVRtoptag) {

        // tighter b-tag selection
        BTag bJetID_medium = BTag(BTag::algo::DEEPJET, BTag::wp::WP_MEDIUM);
        int N_btag_medium = 0;
        for (const auto & jet: *event.jets){
          if(bJetID_medium(jet, event)) N_btag_medium++;
        }

        if(N_btag_medium >= 2){
          h_hists_ttbarCR->fill(event);
          h_ttbarControlRegion_total->fill(event);

          if(event.get(h_flag_muonevent)) {
            h_hists_ttbarCR_mu->fill(event);
            h_ttbarControlRegion_mu->fill(event);
          }
          else {
            
            h_hists_ttbarCR_ele->fill(event);
            h_ttbarControlRegion_ele->fill(event);

            // special case here: plot same, but without electron trigger SF
            double reset_weight_etrigger = event.weight;
            if( event.get(h_weight_sfelec_triggerNominal) != 0) event.weight /= event.get(h_weight_sfelec_triggerNominal);
            h_hists_ttbarCR_noElecTrigSFs->fill(event);
            event.weight = reset_weight_etrigger;

          }

        }

      }

      if (is_datadriven_BG_run) event.weight = weight_for_resetting; // resetting the weight

    }

    // lets fill some other histogram class that will give us the output for each point to check
    if(debug) cout << "Filling DDT test hists" << endl;
    h_DDTtestHists->fill(event, newTaggerResults);

  } else { // here we enter the no b-tag CR

    h_hists_btagCR->fill(event);

    h_ControlRegion_total->fill(event);
    if(event.get(h_flag_muonevent)) {
      h_hists_btagCR_mu->fill(event);
      h_ControlRegion_mu->fill(event);
    }
    else {
      h_hists_btagCR_ele->fill(event);
      h_ControlRegion_ele->fill(event);
    }

    event.set(h_region,"CR");
    if(debug) cout << "Region: no b-tag region" << endl;
  }


  if(debug){cout << "Done ##################################" << endl;}
  return true;

}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarDNNModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarDNNModule)

}
