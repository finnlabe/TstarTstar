#include "UHH2/TstarTstar/include/TstarTstarVariationHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarVariationHists::TstarTstarVariationHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  // get sensitive variable
  h_ST = ctx.get_handle<double>("ST");

  // get weights
  h_weight_sfmu_id = ctx.get_handle<float>("weight_sfmu_id");
  h_weight_sfmu_id_up = ctx.get_handle<float>("weight_sfmu_id_up");
  h_weight_sfmu_id_down = ctx.get_handle<float>("weight_sfmu_id_down");

  h_weight_sfmu_isolation = ctx.get_handle<float>("weight_sfmu_isolation");
  h_weight_sfmu_isolation_up = ctx.get_handle<float>("weight_sfmu_isolation_up");
  h_weight_sfmu_isolation_down = ctx.get_handle<float>("weight_sfmu_isolation_down");

  h_weight_sfelec_id = ctx.get_handle<float>("weight_sfelec_id");
  h_weight_sfelec_id_up = ctx.get_handle<float>("weight_sfelec_id_up");
  h_weight_sfelec_id_down = ctx.get_handle<float>("weight_sfelec_id_down");

  // binning
  const int nbins = 36;
  double bins[nbins] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500,
    2600, 2700, 2800, 2900, 3000, 3200, 3400, 3600, 4000, 6000};

  // nominal
  book<TH1F>("pt_ST_rebinned", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_muonIDUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_muonIDDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_muonIsoUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_muonIsoDown", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_eleIDUp", "S_{T} [GeV]", nbins-1, bins);
  book<TH1F>("pt_ST_rebinned_eleIDDown", "S_{T} [GeV]", nbins-1, bins);

  book<TH1F>("pt_mu", "p^{#mu}_{T} [GeV]",50, 0, 1000);
  book<TH1F>("pt_mu_muonIDUp", "p^{#mu}_{T} [GeV]", 50, 0, 1000);
  book<TH1F>("pt_mu_muonIDDown", "p^{#mu}_{T} [GeV]", 50, 0, 1000);
  book<TH1F>("pt_mu_muonIsoUp", "p^{#mu}_{T} [GeV]", 50, 0, 1000);
  book<TH1F>("pt_mu_muonIsoDown", "p^{#mu}_{T} [GeV]" ,50, 0, 1000);
  book<TH1F>("pt_mu_eleIDUp", "p^{#mu}_{T} [GeV]", 50, 0, 1000);
  book<TH1F>("pt_mu_eleIDDown", "p^{#mu}_{T} [GeV]", 50, 0, 1000);

}


void TstarTstarVariationHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  double st = event.get(h_ST);
  double pt_mu = 9999999;
  if(event.muons->size() > 0) pt_mu = event.muons->at(0).pt();
  if(st > 6000) st = 5999.9; // shouln't be much

  // nominal
  hist("pt_ST_rebinned")->Fill(st, weight);
  hist("pt_mu")->Fill(pt_mu, weight);

  // muon ID
  double muID_nominal = event.get(h_weight_sfmu_id);
  double muID_up = event.get(h_weight_sfmu_id_up);
  double muID_down = event.get(h_weight_sfmu_id_down);
  hist("pt_ST_rebinned_muonIDUp")->Fill(st, weight*muID_up/muID_nominal);
  hist("pt_ST_rebinned_muonIDDown")->Fill(st, weight*muID_down/muID_nominal);
  hist("pt_mu_muonIDUp")->Fill(pt_mu, weight*muID_up/muID_nominal);
  hist("pt_mu_muonIDDown")->Fill(pt_mu, weight*muID_down/muID_nominal);

  // muon isolation
  double muIso_nominal = event.get(h_weight_sfmu_isolation);
  double muIso_up = event.get(h_weight_sfmu_isolation_up);
  double muIso_down = event.get(h_weight_sfmu_isolation_down);
  hist("pt_ST_rebinned_muonIsoUp")->Fill(st, weight*muIso_up/muIso_nominal);
  hist("pt_ST_rebinned_muonIsoDown")->Fill(st, weight*muIso_down/muIso_nominal);
  hist("pt_mu_muonIsoUp")->Fill(pt_mu, weight*muIso_up/muIso_nominal);
  hist("pt_mu_muonIsoDown")->Fill(pt_mu, weight*muIso_down/muIso_nominal);

  // electron ID
  double eleID_nominal = event.get(h_weight_sfelec_id);
  double eleID_up = event.get(h_weight_sfelec_id_up);
  double eleID_down = event.get(h_weight_sfelec_id_down);
  hist("pt_ST_rebinned_eleIDUp")->Fill(st, weight*eleID_up/eleID_nominal);
  hist("pt_ST_rebinned_eleIDDown")->Fill(st, weight*eleID_down/eleID_nominal);
  hist("pt_mu_eleIDUp")->Fill(pt_mu, weight*eleID_up/eleID_nominal);
  hist("pt_mu_eleIDDown")->Fill(pt_mu, weight*eleID_down/eleID_nominal);

  // check if I want to apply top tagging SF here. If so, I would like to know whether event is merged, semi-merged or not!
  // probably its best to check it the HOTVR scaling thingy can output and store that in a handle

}

TstarTstarVariationHists::~TstarTstarVariationHists(){}
