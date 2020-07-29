#include "UHH2/TstarTstar/include/TstarTstarDNNInputHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <string>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarDNNInputHists::TstarTstarDNNInputHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  h_DNN_Inputs = ctx.get_handle<std::vector<double>>("DNN_Inputs");

  // book all histograms here
  DNN_2D_1 = book<TH1D>("DNN_Input_1", "lepton pt", 20, 0, 2000);
  DNN_2D_2 = book<TH1D>("DNN_Input_2", "lepton eta", 20, -4, 4);
  DNN_2D_3 = book<TH1D>("DNN_Input_3", "lepton phi", 20, -4, 4);
  DNN_2D_4 = book<TH1D>("DNN_Input_4", "lepton energy", 20, 0, 3000);
  DNN_2D_5 = book<TH1D>("DNN_Input_5", "lepton relIso", 20, 0, 1);
  DNN_2D_6 = book<TH1D>("DNN_Input_6", "lepton dR closest jet", 20, 0, 6);
  DNN_2D_7 = book<TH1D>("DNN_Input_7", "lepton rel pt", 20, 0, 5);
  DNN_2D_8 = book<TH1D>("DNN_Input_8", "AK4-1 pt", 20, 0, 2000);
  DNN_2D_9 = book<TH1D>("DNN_Input_9", "AK4-1 eta", 20, -4, 4);
  DNN_2D_10 = book<TH1D>("DNN_Input_10", "AK4-1 phi", 20, -4, 4);
  DNN_2D_11 = book<TH1D>("DNN_Input_11", "AK4-1 inv mass", 20, 0, 800);
  DNN_2D_12 = book<TH1D>("DNN_Input_12", "AK4-1 btag", 20, 0, 1);
  DNN_2D_13 = book<TH1D>("DNN_Input_13", "AK4-2 pt", 20, 0, 2000);
  DNN_2D_14 = book<TH1D>("DNN_Input_14", "AK4-2 eta", 20, -4, 4);
  DNN_2D_15 = book<TH1D>("DNN_Input_15", "AK4-2 phi", 20, -4, 4);
  DNN_2D_16 = book<TH1D>("DNN_Input_16", "AK4-2 inv mass", 20, 0, 800);
  DNN_2D_17 = book<TH1D>("DNN_Input_17", "AK4-2 btag", 20, 0, 1);
  DNN_2D_18 = book<TH1D>("DNN_Input_18", "AK4-3 pt", 20, 0, 2000);
  DNN_2D_19 = book<TH1D>("DNN_Input_19", "AK4-3 eta", 20, -4, 4);
  DNN_2D_20 = book<TH1D>("DNN_Input_20", "AK4-3 phi", 20, -4, 4);
  DNN_2D_21 = book<TH1D>("DNN_Input_21", "AK4-3 inv mass", 20, 0, 800);
  DNN_2D_22 = book<TH1D>("DNN_Input_22", "AK4-3 btag", 20, 0, 1);
  DNN_2D_23 = book<TH1D>("DNN_Input_23", "AK4-4 pt", 20, 0, 2000);
  DNN_2D_24 = book<TH1D>("DNN_Input_24", "AK4-4 eta", 20, -4, 4);
  DNN_2D_25 = book<TH1D>("DNN_Input_25", "AK4-4 phi", 20, -4, 4);
  DNN_2D_26 = book<TH1D>("DNN_Input_26", "AK4-4 inv mass", 20, 0, 800);
  DNN_2D_27 = book<TH1D>("DNN_Input_27", "AK4-4 btag", 20, 0, 1);
  DNN_2D_28 = book<TH1D>("DNN_Input_28", "HOTVR-1 pt", 20, 0, 2000);
  DNN_2D_29 = book<TH1D>("DNN_Input_29", "HOTVR-1 eta", 20, -4, 4);
  DNN_2D_30 = book<TH1D>("DNN_Input_30", "HOTVR-1 phi", 20, -4, 4);
  DNN_2D_31 = book<TH1D>("DNN_Input_31", "HOTVR-1 inv mass", 20, 0, 1200);
  DNN_2D_32 = book<TH1D>("DNN_Input_32", "HOTVR-1 tau1", 20, 0, 1);
  DNN_2D_33 = book<TH1D>("DNN_Input_33", "HOTVR-1 tau2", 20, 0, 1);
  DNN_2D_34 = book<TH1D>("DNN_Input_34", "HOTVR-1 tau3", 20, 0, 1);
  DNN_2D_35 = book<TH1D>("DNN_Input_35", "HOTVR-1 subjets", 10, 0, 10);
  DNN_2D_36 = book<TH1D>("DNN_Input_36", "HOTVR-2 pt", 20, 0, 2000);
  DNN_2D_37 = book<TH1D>("DNN_Input_37", "HOTVR-2 eta", 20, -4, 4);
  DNN_2D_38 = book<TH1D>("DNN_Input_38", "HOTVR-2 phi", 20, -4, 4);
  DNN_2D_39 = book<TH1D>("DNN_Input_39", "HOTVR-2 inv mass", 20, 0, 1200);
  DNN_2D_40 = book<TH1D>("DNN_Input_40", "HOTVR-2 tau1", 20, 0, 1);
  DNN_2D_41 = book<TH1D>("DNN_Input_41", "HOTVR-2 tau2", 20, 0, 1);
  DNN_2D_42 = book<TH1D>("DNN_Input_42", "HOTVR-2 tau2", 20, 0, 1);
  DNN_2D_43 = book<TH1D>("DNN_Input_43", "HOTVR-2 subjets", 10, 0, 10);
  DNN_2D_44 = book<TH1D>("DNN_Input_44", "Neutrino pt", 20, 0, 2000);
  DNN_2D_45 = book<TH1D>("DNN_Input_45", "Neutrino eta", 20, -4, 4);
  DNN_2D_46 = book<TH1D>("DNN_Input_46", "Neutrino phi", 20, -4, 4);
  DNN_2D_47 = book<TH1D>("DNN_Input_47", "AK4 count", 10, 0, 10);
  DNN_2D_48 = book<TH1D>("DNN_Input_48", "HOTVR count", 10, 0, 10);
  DNN_2D_49 = book<TH1D>("DNN_Input_49", "Sphericity px px", 20, -1, 1);
  DNN_2D_50 = book<TH1D>("DNN_Input_50", "Sphericity px py", 20, -1, 1);
  DNN_2D_51 = book<TH1D>("DNN_Input_51", "Sphericity px pz", 20, -1, 1);
  DNN_2D_52 = book<TH1D>("DNN_Input_52", "Sphericity py py", 20, -1, 1);
  DNN_2D_53 = book<TH1D>("DNN_Input_53", "Sphericity py pz", 20, -1, 1);
  DNN_2D_54 = book<TH1D>("DNN_Input_54", "Sphericity pz pz", 20, -1, 1);

}


void TstarTstarDNNInputHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  std::vector<double> inputs = event.get(h_DNN_Inputs);

  DNN_2D_1->Fill(inputs.at(0), weight);
  DNN_2D_2->Fill(inputs.at(1), weight);
  DNN_2D_3->Fill(inputs.at(2), weight);
  DNN_2D_4->Fill(inputs.at(3), weight);
  DNN_2D_5->Fill(inputs.at(4), weight);
  DNN_2D_6->Fill(inputs.at(5), weight);
  DNN_2D_7->Fill(inputs.at(6), weight);
  DNN_2D_8->Fill(inputs.at(7), weight);
  DNN_2D_9->Fill(inputs.at(8), weight);
  DNN_2D_10->Fill(inputs.at(9), weight);
  DNN_2D_11->Fill(inputs.at(10), weight);
  DNN_2D_12->Fill(inputs.at(11), weight);
  DNN_2D_13->Fill(inputs.at(12), weight);
  DNN_2D_14->Fill(inputs.at(13), weight);
  DNN_2D_15->Fill(inputs.at(14), weight);
  DNN_2D_16->Fill(inputs.at(15), weight);
  DNN_2D_17->Fill(inputs.at(16), weight);
  DNN_2D_18->Fill(inputs.at(17), weight);
  DNN_2D_19->Fill(inputs.at(18), weight);
  DNN_2D_20->Fill(inputs.at(19), weight);
  DNN_2D_21->Fill(inputs.at(20), weight);
  DNN_2D_22->Fill(inputs.at(21), weight);
  DNN_2D_23->Fill(inputs.at(22), weight);
  DNN_2D_24->Fill(inputs.at(23), weight);
  DNN_2D_25->Fill(inputs.at(24), weight);
  DNN_2D_26->Fill(inputs.at(25), weight);
  DNN_2D_27->Fill(inputs.at(26), weight);
  DNN_2D_28->Fill(inputs.at(27), weight);
  DNN_2D_29->Fill(inputs.at(28), weight);
  DNN_2D_30->Fill(inputs.at(29), weight);
  DNN_2D_31->Fill(inputs.at(30), weight);
  DNN_2D_32->Fill(inputs.at(31), weight);
  DNN_2D_33->Fill(inputs.at(32), weight);
  DNN_2D_34->Fill(inputs.at(33), weight);
  DNN_2D_35->Fill(inputs.at(34), weight);
  DNN_2D_36->Fill(inputs.at(35), weight);
  DNN_2D_37->Fill(inputs.at(36), weight);
  DNN_2D_38->Fill(inputs.at(37), weight);
  DNN_2D_39->Fill(inputs.at(38), weight);
  DNN_2D_40->Fill(inputs.at(39), weight);
  DNN_2D_41->Fill(inputs.at(40), weight);
  DNN_2D_42->Fill(inputs.at(41), weight);
  DNN_2D_43->Fill(inputs.at(42), weight);
  DNN_2D_44->Fill(inputs.at(43), weight);
  DNN_2D_45->Fill(inputs.at(44), weight);
  DNN_2D_46->Fill(inputs.at(45), weight);
  DNN_2D_47->Fill(inputs.at(46), weight);
  DNN_2D_48->Fill(inputs.at(47), weight);
  DNN_2D_49->Fill(inputs.at(48), weight);
  DNN_2D_50->Fill(inputs.at(49), weight);
  DNN_2D_51->Fill(inputs.at(50), weight);
  DNN_2D_52->Fill(inputs.at(51), weight);
  DNN_2D_53->Fill(inputs.at(52), weight);
  DNN_2D_54->Fill(inputs.at(52), weight);

}

TstarTstarDNNInputHists::~TstarTstarDNNInputHists(){}
