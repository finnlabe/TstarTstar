#include "UHH2/TstarTstar/include/TstarTstarAfterDNNHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

float inv_mass_2(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }


TstarTstarAfterDNNHists::TstarTstarAfterDNNHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_DNN_output = ctx.get_handle<double>("DNN_output");

  // book all histograms here
  book<TH1F>("DNN_output", "DNN output", 20, 0, 1);

  book<TH1F>("tau32", "HOTVR-1 #tau_{32}", 20, 0, 1);
  book<TH1F>("jetmass", "m_{HOTVR-1}", 25, 0, 500);

}


void TstarTstarAfterDNNHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  double DNNoutput = event.get(h_DNN_output);

  hist("DNN_output")->Fill(DNNoutput, weight);

  hist("tau32")->Fill(event.topjets->at(0).tau3_groomed()/event.topjets->at(0).tau2_groomed(), weight);
  hist("jetmass")->Fill(inv_mass_2(event.topjets->at(0).v4()), weight);

}

TstarTstarAfterDNNHists::~TstarTstarAfterDNNHists(){}
