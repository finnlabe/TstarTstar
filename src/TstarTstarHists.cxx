#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

TstarTstarHists::TstarTstarHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1F>("N_jets", "N_{AK4 jets}", 20, 0, 20);  
  book<TH1F>("N_AK8jets", "N_{AK8 jets}", 20, 0, 20);  
  book<TH1F>("N_PU", "N_{PU}", 100, 0, 100);  
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 50, -5.2, 5.2);
  book<TH1F>("pt_jet1", "p_{T}^{jet 1}", 100, 10, 2000);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 50, -5.2, 5.2);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 50, -5.2, 5.2);
  book<TH1F>("eta_jet4", "#eta^{jet 4}", 50, -5.2, 5.2);
  book<TH1F>("pt_jet2", "p_{T}^{jet 2}", 100, 10, 2000);
  book<TH1F>("pt_jet3", "p_{T}^{jet 3}", 100, 10, 2000);
  book<TH1F>("pt_jet4", "p_{T}^{jet 4}", 100, 10, 2000);

  book<TH1F>("EMcharged_jet1", "EMcharged_jet1", 100,0.0,1.0);
  book<TH1F>("EMneutral_jet1", "EMneutral_jet1", 100,0.0,1.0);
  book<TH1F>("HADcharged_jet1", "HADcharged_jet1", 100,0.0,1.0);
  book<TH1F>("HADneutral_jet1", "HADneutral_jet1", 100,0.0,1.0);

  book<TH1F>("eta_ak8jet1", "#eta^{ak8jet 1}", 50, -5.2, 5.2);
  book<TH1F>("pt_ak8jet1", "p_{T}^{ak8jet 1}", 100, 10, 2000);
  book<TH1F>("eta_ak8jet2", "#eta^{ak8jet 2}", 50, -5.2, 5.2);
  book<TH1F>("pt_ak8jet2", "p_{T}^{ak8jet 2}", 100, 10, 2000);
  book<TH1F>("eta_ak8jet3", "#eta^{ak8jet 3}", 50, -5.2, 5.2);
  book<TH1F>("pt_ak8jet3", "p_{T}^{ak8jet 3}", 100, 10, 2000);

  book<TH2D>("pt_mu_pt_ak4jet1", ";p_{T}^{#mu}; p_{T}^{AK4 jet 1}", 60, 10, 1100, 70, 10, 2500);
  book<TH2D>("pt_ele_pt_ak4jet1", ";p_{T}^{ele}; p_{T}^{AK4 jet 1}", 60, 10, 1100, 70, 10, 2500);
  book<TH2D>("pt_mu_pt_ak8jet1", ";p_{T}^{#mu}; p_{T}^{AK8 jet 1}", 60, 10, 1100, 70, 10, 2500);
  book<TH2D>("pt_ele_pt_ak8jet1", ";p_{T}^{ele}; p_{T}^{AK8 jet 1}", 60, 10, 1100, 70, 10, 2500);

  // book<TH2D>("EMcharged_vs_eta_jet1","EMcharged vs #eta; #eta; EMcharged",100,-6,6,100,0.0,1.0);   
  // book<TH2D>("EMneutral_vs_eta_jet1","EMneutral vs #eta; #eta; EMneutral",100,-6,6,100,0.0,1.0);   
  // book<TH2D>("HADcharged_vs_eta_jet1","HADcharged vs #eta; #eta; HADcharged",100,-6,6,100,0.0,1.0);   
  // book<TH2D>("HADneutral_vs_eta_jet1","HADneutral vs #eta; #eta; HADneutral",100,-6,6,100,0.0,1.0);   
  // book<TH2D>("EMcharged_vs_PU_jet1","EMcharged vs PU; PU; EMcharged",100,0,100,100,0.0,1.0);   
  // book<TH2D>("EMneutral_vs_PU_jet1","EMneutral vs PU; PU; EMneutral",100,0,100,100,0.0,1.0);   
  // book<TH2D>("HADcharged_vs_PU_jet1","HADcharged vs PU; PU; HADcharged",100,0,100,100,0.0,1.0);   
  // book<TH2D>("HADneutral_vs_PU_jet1","HADneutral vs PU; PU; HADneutral",100,0,100,100,0.0,1.0);   


  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 50, 10, 1000);
  book<TH1F>("eta_mu", "#eta^{#mu}", 50, -2.5, 2.5);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  book<TH1F>("N_ele", "N^{e}", 10, 0, 10);
  book<TH1F>("pt_ele", "p_{T}^{ele} [GeV/c]", 50, 10, 1000);
  book<TH1F>("eta_ele", "#eta^{ele}", 50, -2.5, 2.5);

  book<TH1F>("N_photon", "N^{#gamma}", 10, 0, 10);
  book<TH1F>("pt_photon", "p_{T}^{#gamma} [GeV/c]", 50, 10, 2000);
  book<TH1F>("eta_photon", "#eta^{#gamma}", 50, -5.2, 5.2);

  book<TH1F>("pt_photon_1", "p_{T}^{leading #gamma} [GeV/c]", 50, 10, 2000);
  book<TH1F>("pt_photon_2", "p_{T}^{second #gamma} [GeV/c]", 50, 10, 2000);


  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);

  //MET and HT
  book<TH1F>("pt_MET", "missing E_{T} [GeV]", 100, 10, 2000);
  book<TH1F>("pt_ST_jets", "S_{T}^{jets}=#sum_{i}|AK4 p^{i}_{T}| [GeV/c]", 100, 10, 4000);

  //M_Tstar
  h_M_Tstar_gluon_ = ctx.get_handle< float >("M_Tstar_gluon");
  h_M_Tstar_gamma_ = ctx.get_handle< float >("M_Tstar_gamma");
  book<TH1F>("M_Tstar_gluon", "M_{T^{*}_{g}} [GeV]", 30, 0, 3000);
  book<TH1F>("M_Tstar_gamma", "M_{T^{*}_{#gamma}} [GeV]", 30, 0, 3000);

  //Delta_R
  h_DeltaR_toplep_ak8jet1_ = ctx.get_handle< float >("DeltaR_toplep_ak8jet1");
  h_DeltaR_tophad_ak8jet1_ = ctx.get_handle< float >("DeltaR_tophad_ak8jet1");
  h_DeltaR_toplep_ak8jet2_ = ctx.get_handle< float >("DeltaR_toplep_ak8jet2");
  h_DeltaR_tophad_ak8jet2_ = ctx.get_handle< float >("DeltaR_tophad_ak8jet2");
  book<TH1F>("DeltaR_toplep_ak8jet1", "DeltaR_toplep_ak8jet1", 40, 0, 4);
  book<TH1F>("DeltaR_tophad_ak8jet1", "DeltaR_tophad_ak8jet1", 40, 0, 4);
  book<TH1F>("DeltaR_toplep_ak8jet2", "DeltaR_toplep_ak8jet2", 40, 0, 4);
  book<TH1F>("DeltaR_tophad_ak8jet2", "DeltaR_tophad_ak8jet2", 40, 0, 4);
}


void TstarTstarHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);


  if(!event.isRealData)  hist("N_PU")->Fill(event.genInfo->pileup_TrueNumInteractions(), weight);

  if(Njets>=1){
    hist("eta_jet1")->Fill(jets->at(0).eta(), weight);
    hist("pt_jet1")->Fill(jets->at(0).pt(), weight);
    hist("EMcharged_jet1")->Fill(jets->at(0).chargedEmEnergyFraction(), weight);
    hist("EMneutral_jet1")->Fill(jets->at(0).neutralEmEnergyFraction(), weight);
    hist("HADcharged_jet1")->Fill(jets->at(0).chargedHadronEnergyFraction(), weight);
    hist("HADneutral_jet1")->Fill(jets->at(0).neutralHadronEnergyFraction(), weight);
    
    // ((TH2D*)hist("EMcharged_vs_eta_jet1"))->Fill(jets->at(0).eta(),jets->at(0).chargedEmEnergyFraction(), weight);
    // ((TH2D*)hist("EMneutral_vs_eta_jet1"))->Fill(jets->at(0).eta(),jets->at(0).neutralEmEnergyFraction(), weight);
    // ((TH2D*)hist("HADcharged_vs_eta_jet1"))->Fill(jets->at(0).eta(),jets->at(0).chargedHadronEnergyFraction(), weight);
    // ((TH2D*)hist("HADneutral_vs_eta_jet1"))->Fill(jets->at(0).eta(),jets->at(0).neutralHadronEnergyFraction(), weight);
    // if(!event.isRealData){
    //   ((TH2D*)hist("EMcharged_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(),jets->at(0).chargedEmEnergyFraction(), weight);
    //   ((TH2D*)hist("EMneutral_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(),jets->at(0).neutralEmEnergyFraction(), weight);
    //   ((TH2D*)hist("HADcharged_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(),jets->at(0).chargedHadronEnergyFraction(), weight);
    //   ((TH2D*)hist("HADneutral_vs_PU_jet1"))->Fill(event.genInfo->pileup_TrueNumInteractions(),jets->at(0).neutralHadronEnergyFraction(), weight);
    // }
  }
  if(Njets>=2){
    hist("eta_jet2")->Fill(jets->at(1).eta(), weight);
    hist("pt_jet2")->Fill(jets->at(1).pt(), weight);
  }
  if(Njets>=3){
    hist("eta_jet3")->Fill(jets->at(2).eta(), weight);
    hist("pt_jet3")->Fill(jets->at(2).pt(), weight);
  }
  if(Njets>=4){
    hist("eta_jet4")->Fill(jets->at(3).eta(), weight);
    hist("pt_jet4")->Fill(jets->at(3).pt(), weight);
  }

  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      hist("pt_mu")->Fill(thismu.pt(), weight);
      hist("eta_mu")->Fill(thismu.eta(), weight);
      hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }
  
  int Nele = event.electrons->size();
  hist("N_ele")->Fill(Nele, weight);
  for (const Electron & thisele : *event.electrons){
      hist("pt_ele")->Fill(thisele.pt(), weight);
      hist("eta_ele")->Fill(thisele.eta(), weight);
  }

  int Ngamma = event.photons->size();
  hist("N_photon")->Fill(Ngamma, weight);
  int countgamma = 1;
  for (const Photon & thisgamma : *event.photons){
    hist("pt_photon")->Fill(thisgamma.pt(), weight);
    if(countgamma == 1){hist("pt_photon_1")->Fill(thisgamma.pt(), weight);}
    if(countgamma == 2){hist("pt_photon_2")->Fill(thisgamma.pt(), weight);}
    hist("eta_photon")->Fill(thisgamma.eta(), weight);
    countgamma++;
  }

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);

  hist("N_AK8jets")->Fill(event.topjets->size(), weight);
  if(event.topjets->size()>0){
    hist("eta_ak8jet1")->Fill(event.topjets->at(0).eta(), weight);
    hist("pt_ak8jet1")->Fill(event.topjets->at(0).pt(), weight);
  }
  if(event.topjets->size()>1){
    hist("eta_ak8jet2")->Fill(event.topjets->at(1).eta(), weight);
    hist("pt_ak8jet2")->Fill(event.topjets->at(1).pt(), weight);
  }
  if(event.topjets->size()>2){
    hist("eta_ak8jet3")->Fill(event.topjets->at(2).eta(), weight);
    hist("pt_ak8jet3")->Fill(event.topjets->at(2).pt(), weight);
  }
  if(Nele>0 && Njets>0) ((TH2D*)hist("pt_ele_pt_ak4jet1"))->Fill(event.electrons->at(0).pt(),jets->at(0).pt(), weight);
  if(Nmuons>0 && Njets>0) ((TH2D*)hist("pt_mu_pt_ak4jet1"))->Fill(event.muons->at(0).pt(),jets->at(0).pt(), weight);
  if(Nele>0 && event.topjets->size()>0) ((TH2D*)hist("pt_ele_pt_ak8jet1"))->Fill(event.electrons->at(0).pt(),event.topjets->at(0).pt(), weight);
  if(Nmuons>0 && event.topjets->size()>0) ((TH2D*)hist("pt_mu_pt_ak8jet1"))->Fill(event.muons->at(0).pt(),event.topjets->at(0).pt(), weight);

  hist("pt_MET")->Fill(event.met->pt(), weight);

  double st_jets = 0.;
  for(const auto & jet : *jets) st_jets += jet.pt();
  hist("pt_ST_jets")->Fill(st_jets, weight);

  hist("M_Tstar_gluon")->Fill(event.get(h_M_Tstar_gluon_), weight);
  hist("M_Tstar_gamma")->Fill(event.get(h_M_Tstar_gamma_), weight);

  float val_ = event.get(h_DeltaR_toplep_ak8jet1_);
  float binwidth_ = hist("DeltaR_toplep_ak8jet1")->GetXaxis()->GetBinWidth( hist("DeltaR_toplep_ak8jet1")->GetXaxis()->FindBin(val_) );
  hist("DeltaR_toplep_ak8jet1")->Fill(val_, weight / binwidth_);

  val_ = event.get(h_DeltaR_tophad_ak8jet1_);
  binwidth_ = hist("DeltaR_tophad_ak8jet1")->GetXaxis()->GetBinWidth( hist("DeltaR_tophad_ak8jet1")->GetXaxis()->FindBin(val_) );
  hist("DeltaR_tophad_ak8jet1")->Fill(val_, weight / binwidth_);

  val_ = event.get(h_DeltaR_toplep_ak8jet2_);
  binwidth_ = hist("DeltaR_toplep_ak8jet2")->GetXaxis()->GetBinWidth( hist("DeltaR_toplep_ak8jet2")->GetXaxis()->FindBin(val_) );
  hist("DeltaR_toplep_ak8jet2")->Fill(val_, weight / binwidth_);

  val_ = event.get(h_DeltaR_tophad_ak8jet2_);
  binwidth_ = hist("DeltaR_tophad_ak8jet2")->GetXaxis()->GetBinWidth( hist("DeltaR_tophad_ak8jet2")->GetXaxis()->FindBin(val_) );
  hist("DeltaR_tophad_ak8jet2")->Fill(val_, weight / binwidth_);

}

TstarTstarHists::~TstarTstarHists(){}
