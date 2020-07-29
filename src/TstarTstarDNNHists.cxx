#include "UHH2/TstarTstar/include/TstarTstarDNNHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarDNNHists::TstarTstarDNNHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_DNN_output = ctx.get_handle<double>("DNN_output");
  h_DNN_Inputs = ctx.get_handle<std::vector<double>>("DNN_Inputs");

  // book all histograms here
  book<TH1F>("DNN_output", "DNN output", 20, 0, 1);
  DNN_2D_ST = book<TH2D>("2D_DNN_ST", "DNN output against ST", 20, 0, 3000, 20, 0, 1);
  DNN_2D_1 = book<TH2D>("2D_DNN_1", "DNN output against lepton pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_2 = book<TH2D>("2D_DNN_2", "DNN output against lepton eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_3 = book<TH2D>("2D_DNN_3", "DNN output against lepton phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_4 = book<TH2D>("2D_DNN_4", "DNN output against lepton energy", 20, -1, 1, 20, 0, 1);
  DNN_2D_5 = book<TH2D>("2D_DNN_5", "DNN output against lepton relIso", 20, -1, 1, 20, 0, 1);
  DNN_2D_6 = book<TH2D>("2D_DNN_6", "DNN output against lepton dR closest jet", 20, -1, 1, 20, 0, 1);
  DNN_2D_7 = book<TH2D>("2D_DNN_7", "DNN output against lepton rel pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_8 = book<TH2D>("2D_DNN_8", "DNN output against AK4-1 pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_9 = book<TH2D>("2D_DNN_9", "DNN output against AK4-1 eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_10 = book<TH2D>("2D_DNN_10", "DNN output against AK4-1 phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_11 = book<TH2D>("2D_DNN_11", "DNN output against AK4-1 inv mass", 20, -1, 1, 20, 0, 1);
  DNN_2D_12 = book<TH2D>("2D_DNN_12", "DNN output against AK4-1 btag", 20, -1, 1, 20, 0, 1);
  DNN_2D_13 = book<TH2D>("2D_DNN_13", "DNN output against AK4-2 pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_14 = book<TH2D>("2D_DNN_14", "DNN output against AK4-2 eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_15 = book<TH2D>("2D_DNN_15", "DNN output against AK4-2 phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_16 = book<TH2D>("2D_DNN_16", "DNN output against AK4-2 inv mass", 20, -1, 1, 20, 0, 1);
  DNN_2D_17 = book<TH2D>("2D_DNN_17", "DNN output against AK4-2 btag", 20, -1, 1, 20, 0, 1);
  DNN_2D_18 = book<TH2D>("2D_DNN_18", "DNN output against AK4-3 pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_19 = book<TH2D>("2D_DNN_19", "DNN output against AK4-3 eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_20 = book<TH2D>("2D_DNN_20", "DNN output against AK4-3 phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_21 = book<TH2D>("2D_DNN_21", "DNN output against AK4-3 inv mass", 20, -1, 1, 20, 0, 1);
  DNN_2D_22 = book<TH2D>("2D_DNN_22", "DNN output against AK4-3 btag", 20, -1, 1, 20, 0, 1);
  DNN_2D_23 = book<TH2D>("2D_DNN_23", "DNN output against AK4-4 pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_24 = book<TH2D>("2D_DNN_24", "DNN output against AK4-4 eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_25 = book<TH2D>("2D_DNN_25", "DNN output against AK4-4 phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_26 = book<TH2D>("2D_DNN_26", "DNN output against AK4-4 inv mass", 20, -1, 1, 20, 0, 1);
  DNN_2D_27 = book<TH2D>("2D_DNN_27", "DNN output against AK4-4 btag", 20, -1, 1, 20, 0, 1);
  DNN_2D_28 = book<TH2D>("2D_DNN_28", "DNN output against HOTVR-1 pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_29 = book<TH2D>("2D_DNN_29", "DNN output against HOTVR-1 eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_30 = book<TH2D>("2D_DNN_30", "DNN output against HOTVR-1 phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_31 = book<TH2D>("2D_DNN_31", "DNN output against HOTVR-1 inv mass", 20, -1, 1, 20, 0, 1);
  DNN_2D_32 = book<TH2D>("2D_DNN_32", "DNN output against HOTVR-1 tau1", 20, -1, 1, 20, 0, 1);
  DNN_2D_33 = book<TH2D>("2D_DNN_33", "DNN output against HOTVR-1 tau2", 20, -1, 1, 20, 0, 1);
  DNN_2D_34 = book<TH2D>("2D_DNN_34", "DNN output against HOTVR-1 tau3", 20, -1, 1, 20, 0, 1);
  DNN_2D_35 = book<TH2D>("2D_DNN_35", "DNN output against HOTVR-1 subjets", 20, -1, 1, 20, 0, 1);
  DNN_2D_36 = book<TH2D>("2D_DNN_36", "DNN output against HOTVR-1 pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_37 = book<TH2D>("2D_DNN_37", "DNN output against HOTVR-1 eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_38 = book<TH2D>("2D_DNN_38", "DNN output against HOTVR-1 phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_39 = book<TH2D>("2D_DNN_39", "DNN output against HOTVR-1 inv mass", 20, -1, 1, 20, 0, 1);
  DNN_2D_40 = book<TH2D>("2D_DNN_40", "DNN output against HOTVR-1 tau1", 20, -1, 1, 20, 0, 1);
  DNN_2D_41 = book<TH2D>("2D_DNN_41", "DNN output against HOTVR-1 tau2", 20, -1, 1, 20, 0, 1);
  DNN_2D_42 = book<TH2D>("2D_DNN_42", "DNN output against HOTVR-1 tau2", 20, -1, 1, 20, 0, 1);
  DNN_2D_43 = book<TH2D>("2D_DNN_43", "DNN output against HOTVR-1 sunjets", 20, -1, 1, 20, 0, 1);
  DNN_2D_44 = book<TH2D>("2D_DNN_44", "DNN output against Neutrino pt", 20, -1, 1, 20, 0, 1);
  DNN_2D_45 = book<TH2D>("2D_DNN_45", "DNN output against Neutrino eta", 20, -1, 1, 20, 0, 1);
  DNN_2D_46 = book<TH2D>("2D_DNN_46", "DNN output against Neutrino phi", 20, -1, 1, 20, 0, 1);
  DNN_2D_47 = book<TH2D>("2D_DNN_47", "DNN output against AK4 count", 20, -1, 1, 20, 0, 1);
  DNN_2D_48 = book<TH2D>("2D_DNN_48", "DNN output against HOTVR count", 20, -1, 1, 20, 0, 1);
  DNN_2D_49 = book<TH2D>("2D_DNN_49", "DNN output against Sphericity px px", 20, -1, 1, 20, 0, 1);
  DNN_2D_50 = book<TH2D>("2D_DNN_50", "DNN output against Sphericity px py", 20, -1, 1, 20, 0, 1);
  DNN_2D_51 = book<TH2D>("2D_DNN_51", "DNN output against Sphericity px pz", 20, -1, 1, 20, 0, 1);
  DNN_2D_52 = book<TH2D>("2D_DNN_52", "DNN output against Sphericity py py", 20, -1, 1, 20, 0, 1);
  DNN_2D_53 = book<TH2D>("2D_DNN_53", "DNN output against Sphericity py pz", 20, -1, 1, 20, 0, 1);
  DNN_2D_54 = book<TH2D>("2D_DNN_54", "DNN output against Sphericity pz pz", 20, -1, 1, 20, 0, 1);

}


void TstarTstarDNNHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  double DNNoutput = event.get(h_DNN_output);
  std::vector<double> inputs = event.get(h_DNN_Inputs);

  hist("DNN_output")->Fill(DNNoutput, weight);
  double st_jets = 0;
  for(const auto & jet : *event.topjets) st_jets += jet.pt();
  DNN_2D_ST->Fill(st_jets, DNNoutput, weight);

  DNN_2D_1->Fill(inputs.at(0), DNNoutput, weight);
  DNN_2D_2->Fill(inputs.at(1), DNNoutput, weight);
  DNN_2D_3->Fill(inputs.at(2), DNNoutput, weight);
  DNN_2D_4->Fill(inputs.at(3), DNNoutput, weight);
  DNN_2D_5->Fill(inputs.at(4), DNNoutput, weight);
  DNN_2D_6->Fill(inputs.at(5), DNNoutput, weight);
  DNN_2D_7->Fill(inputs.at(6), DNNoutput, weight);
  DNN_2D_8->Fill(inputs.at(7), DNNoutput, weight);
  DNN_2D_9->Fill(inputs.at(8), DNNoutput, weight);
  DNN_2D_10->Fill(inputs.at(9), DNNoutput, weight);
  DNN_2D_11->Fill(inputs.at(10), DNNoutput, weight);
  DNN_2D_12->Fill(inputs.at(11), DNNoutput, weight);
  DNN_2D_13->Fill(inputs.at(12), DNNoutput, weight);
  DNN_2D_14->Fill(inputs.at(13), DNNoutput, weight);
  DNN_2D_15->Fill(inputs.at(14), DNNoutput, weight);
  DNN_2D_16->Fill(inputs.at(15), DNNoutput, weight);
  DNN_2D_17->Fill(inputs.at(16), DNNoutput, weight);
  DNN_2D_18->Fill(inputs.at(17), DNNoutput, weight);
  DNN_2D_19->Fill(inputs.at(18), DNNoutput, weight);
  DNN_2D_20->Fill(inputs.at(19), DNNoutput, weight);
  DNN_2D_21->Fill(inputs.at(20), DNNoutput, weight);
  DNN_2D_22->Fill(inputs.at(21), DNNoutput, weight);
  DNN_2D_23->Fill(inputs.at(22), DNNoutput, weight);
  DNN_2D_24->Fill(inputs.at(23), DNNoutput, weight);
  DNN_2D_25->Fill(inputs.at(24), DNNoutput, weight);
  DNN_2D_26->Fill(inputs.at(25), DNNoutput, weight);
  DNN_2D_27->Fill(inputs.at(26), DNNoutput, weight);
  DNN_2D_28->Fill(inputs.at(27), DNNoutput, weight);
  DNN_2D_29->Fill(inputs.at(28), DNNoutput, weight);
  DNN_2D_30->Fill(inputs.at(29), DNNoutput, weight);
  DNN_2D_31->Fill(inputs.at(30), DNNoutput, weight);
  DNN_2D_32->Fill(inputs.at(31), DNNoutput, weight);
  DNN_2D_33->Fill(inputs.at(32), DNNoutput, weight);
  DNN_2D_34->Fill(inputs.at(33), DNNoutput, weight);
  DNN_2D_35->Fill(inputs.at(34), DNNoutput, weight);
  DNN_2D_36->Fill(inputs.at(35), DNNoutput, weight);
  DNN_2D_37->Fill(inputs.at(36), DNNoutput, weight);
  DNN_2D_38->Fill(inputs.at(37), DNNoutput, weight);
  DNN_2D_39->Fill(inputs.at(38), DNNoutput, weight);
  DNN_2D_40->Fill(inputs.at(39), DNNoutput, weight);
  DNN_2D_41->Fill(inputs.at(40), DNNoutput, weight);
  DNN_2D_42->Fill(inputs.at(41), DNNoutput, weight);
  DNN_2D_43->Fill(inputs.at(42), DNNoutput, weight);
  DNN_2D_44->Fill(inputs.at(43), DNNoutput, weight);
  DNN_2D_45->Fill(inputs.at(44), DNNoutput, weight);
  DNN_2D_46->Fill(inputs.at(45), DNNoutput, weight);
  DNN_2D_47->Fill(inputs.at(46), DNNoutput, weight);
  DNN_2D_48->Fill(inputs.at(47), DNNoutput, weight);
  DNN_2D_49->Fill(inputs.at(48), DNNoutput, weight);
  DNN_2D_50->Fill(inputs.at(49), DNNoutput, weight);
  DNN_2D_51->Fill(inputs.at(50), DNNoutput, weight);
  DNN_2D_52->Fill(inputs.at(51), DNNoutput, weight);
  DNN_2D_53->Fill(inputs.at(52), DNNoutput, weight);
  DNN_2D_54->Fill(inputs.at(53), DNNoutput, weight);

}

TstarTstarDNNHists::~TstarTstarDNNHists(){}
