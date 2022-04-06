#include "UHH2/TstarTstar/include/TstarTstarSignalRegionHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

TstarTstarSignalRegionHists::TstarTstarSignalRegionHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  // ST handle
  h_ST = ctx.get_handle<double>("ST");

  // weight handles
  h_weight_puNominal = ctx.get_handle<double>("weight_pu");
  h_weight_puUp = ctx.get_handle<double>("weight_pu_up");
  h_weight_puDown = ctx.get_handle<double>("weight_pu_down");

  h_prefiringWeightNominal = ctx.get_handle<double>("prefiringWeight");
  h_prefiringWeightUp = ctx.get_handle<double>("prefiringWeightUp");
  h_prefiringWeightDown = ctx.get_handle<double>("prefiringWeightDown");

  h_weight_btagdiscNominal = ctx.get_handle<double>("weight_btagdisc_central");

  h_weight_btagdisc_jesUp = ctx.get_handle<double>("weight_btagdisc_jes_up");
  h_weight_btagdisc_jesDown = ctx.get_handle<double>("weight_btagdisc_jes_down");

  h_weight_btagdisc_lfUp = ctx.get_handle<double>("weight_btagdisc_lf_up");
  h_weight_btagdisc_lfDown = ctx.get_handle<double>("weight_btagdisc_lf_down");

  h_weight_btagdisc_hfUp = ctx.get_handle<double>("weight_btagdisc_hf_up");
  h_weight_btagdisc_hfDown = ctx.get_handle<double>("weight_btagdisc_hf_down");

  h_weight_btagdisc_hfstats1Up = ctx.get_handle<double>("weight_btagdisc_hfstats1_up");
  h_weight_btagdisc_hfstats1Down = ctx.get_handle<double>("weight_btagdisc_hfstats1_down");

  h_weight_btagdisc_hfstats2Up = ctx.get_handle<double>("weight_btagdisc_hfstats2_up");
  h_weight_btagdisc_hfstats2Down = ctx.get_handle<double>("weight_btagdisc_hfstats2_down");

  h_weight_btagdisc_lfstats1Up = ctx.get_handle<double>("weight_btagdisc_lfstats1_up");
  h_weight_btagdisc_lfstats1Down = ctx.get_handle<double>("weight_btagdisc_lfstats1_down");

  h_weight_btagdisc_lfstats2Up = ctx.get_handle<double>("weight_btagdisc_lfstats2_up");
  h_weight_btagdisc_lfstats2Down = ctx.get_handle<double>("weight_btagdisc_lfstats2_down");

  h_weight_btagdisc_cferr1Up = ctx.get_handle<double>("weight_btagdisc_cferr1_up");
  h_weight_btagdisc_cferr1Down = ctx.get_handle<double>("weight_btagdisc_cferr1_down");

  h_weight_btagdisc_cferr2Up = ctx.get_handle<double>("weight_btagdisc_cferr2_up");
  h_weight_btagdisc_cferr2Down = ctx.get_handle<double>("weight_btagdisc_cferr2_down");

  h_weight_sfelec_idNominal = ctx.get_handle<double>("weight_sfelec_id");
  h_weight_sfelec_idUp = ctx.get_handle<double>("weight_sfelec_id_up");
  h_weight_sfelec_idDown = ctx.get_handle<double>("weight_sfelec_id_down");

  h_weight_sfelec_triggerNominal = ctx.get_handle<double>("weight_sfelec_trigger");
  h_weight_sfelec_triggerUp = ctx.get_handle<double>("weight_sfelec_trigger_up");
  h_weight_sfelec_triggerDown = ctx.get_handle<double>("weight_sfelec_trigger_down");

  h_weight_sfmu_idNominal = ctx.get_handle<double>("weight_sfmu_id");
  h_weight_sfmu_idUp = ctx.get_handle<double>("weight_sfmu_id_up");
  h_weight_sfmu_idDown = ctx.get_handle<double>("weight_sfmu_id_down");

  h_weight_sfmu_isolationNominal = ctx.get_handle<double>("weight_sfmu_isolation");
  h_weight_sfmu_isolationUp = ctx.get_handle<double>("weight_sfmu_isolation_up");
  h_weight_sfmu_isolationDown = ctx.get_handle<double>("weight_sfmu_isolation_down");

  h_weight_sfmu_triggerNominal = ctx.get_handle<double>("weight_sfmu_trigger");
  h_weight_sfmu_triggerUp = ctx.get_handle<double>("weight_sfmu_trigger_up");
  h_weight_sfmu_triggerDown = ctx.get_handle<double>("weight_sfmu_trigger_down");

  // binning definition
  const int nbins = 34;
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
    2600, 2700, 2800, 2900, 3000, 3250, 4000, 6000};

  // nominal histogram
  book<TH1F>("pt_ST_nominal", "S_{T} [GeV]", nbins-1, bins);

  // variation histograms
  book<TH1F>("pt_ST_puUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_puDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_prefiringUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_prefiringDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_jesUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_jesDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_lfUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_lfDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_hfUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_hfDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_hfstats1Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_hfstats1Down", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_hfstats2Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_hfstats2Down", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_lfstats1Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_lfstats1Down", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_lfstats2Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_lfstats2Down", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_cferr1Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_cferr1Down", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_cferr2Up", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_btagdisc_cferr2Down", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfelec_idUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfelec_idDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfelec_triggerUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfelec_triggerDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_idUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_idDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_isolationUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_isolationDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_triggerUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_sfmu_triggerDown", "S_{T} [GeV]", nbins-1, bins);

}


void TstarTstarSignalRegionHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  bool debug = false;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  if(debug) cout << "Starting Tstar Hists." << endl;

  try { // this will of course only work once ST has been defined

    double st = event.get(h_ST);
    if(st > 6000) st = 5999.9; // handling overflow

    // fill nominal
    hist("pt_ST_nominal")->Fill(st, weight);

    // example dummy code
    // double event_weight_up = event.get(h_some_variation_up);
    // double event_weight_down = event.get(h_some_variation_down);
    // hist("pt_ST_some_variation_up")->Fill(st, event_weight_up);
    // hist("pt_ST_some_variation_down")->Fill(st, event_weight_down);

    // pu
    hist("pt_ST_puUp")->Fill(st, event.get(h_weight_puUp)*weight/event.get(h_weight_puNominal));
    hist("pt_ST_puDown")->Fill(st, event.get(h_weight_puDown)*weight/event.get(h_weight_puNominal));

    // prefiring
    hist("pt_ST_prefiringUp")->Fill(st, event.get(h_prefiringWeightUp)*weight/event.get(h_prefiringWeightNominal));
    hist("pt_ST_prefiringDown")->Fill(st, event.get(h_prefiringWeightDown)*weight/event.get(h_prefiringWeightNominal));

    // put b-tagging stuff here

    // sfelec_id
    hist("pt_ST_sfelec_idUp")->Fill(st, event.get(h_weight_sfelec_idUp)*weight/event.get(h_weight_sfelec_idNominal));
    hist("pt_ST_sfelec_idDown")->Fill(st, event.get(h_weight_sfelec_idDown)*weight/event.get(h_weight_sfelec_idNominal));

    // sfelec_trigger
    hist("pt_ST_sfelec_triggerUp")->Fill(st, event.get(h_weight_sfelec_triggerUp)*weight/event.get(h_weight_sfelec_triggerNominal));
    hist("pt_ST_sfelec_triggerDown")->Fill(st, event.get(h_weight_sfelec_triggerDown)*weight/event.get(h_weight_sfelec_triggerNominal));

    // sfmu_id
    hist("pt_ST_sfmu_idUp")->Fill(st, event.get(h_weight_sfmu_idUp)*weight/event.get(h_weight_sfmu_idNominal));
    hist("pt_ST_sfmu_idDown")->Fill(st, event.get(h_weight_sfmu_idDown)*weight/event.get(h_weight_sfmu_idNominal));

    // sfmu_isolation
    hist("pt_ST_sfmu_isolationUp")->Fill(st, event.get(h_weight_sfmu_isolationUp)*weight/event.get(h_weight_sfmu_isolationNominal));
    hist("pt_ST_sfmu_isolationDown")->Fill(st, event.get(h_weight_sfmu_isolationDown)*weight/event.get(h_weight_sfmu_isolationNominal));

    // sfmu_isolation
    hist("pt_ST_sfmu_triggerUp")->Fill(st, event.get(h_weight_sfmu_triggerUp)*weight/event.get(h_weight_sfmu_triggerNominal));
    hist("pt_ST_sfmu_triggerDown")->Fill(st, event.get(h_weight_sfmu_triggerDown)*weight/event.get(h_weight_sfmu_triggerNominal));

  } catch(...) {}

  }

TstarTstarSignalRegionHists::~TstarTstarSignalRegionHists(){}
