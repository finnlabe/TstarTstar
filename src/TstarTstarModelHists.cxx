#include "UHH2/TstarTstar/include/TstarTstarModelHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }
}


TstarTstarModelHists::TstarTstarModelHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  book<TH1F>("weight", "weight", 40, 0, 2);

  // #################
  // ## Basic hists ##
  // #################

  book<TH1F>("N_Tstar", "N_{T*} gen", 5, 0, 5);
  book<TH1F>("Tstar_spin", "Spin^{T*} gen", 30, -3, 3);
  book<TH1F>("Tstar_mass", "M^{T*} gen", 100, 0, 2000);
  book<TH1F>("Tstar_pt", "p_{T}^{T*} gen", 50, 0, 2000);
  book<TH1F>("Tstar_eta", "#eta^{T*} gen", 50, -5.2, 5.2);
  book<TH1F>("Tstar_phi", "#phi^{T*} gen", 30, -3.14, 3.14);

  book<TH1F>("N_top", "N_{t} gen", 5, 0, 5);
  book<TH1F>("top_spin", "Spin^{t} gen", 30, -3, 3);
  book<TH1F>("top_mass", "M^{t} gen", 100, 0, 2000);
  book<TH1F>("top_pt", "p_{T}^{t} gen", 50, 0, 2000);
  book<TH1F>("top_pt_initial", "p_{T}^{t (initial)} gen", 50, 0, 2000);
  book<TH1F>("top_pt_fromTstar", "p_{T}^{t (from T*)} gen", 50, 0, 2000);
  book<TH1F>("top_eta", "#eta^{t} gen", 52, -5.2, 5.2);
  book<TH1F>("top_phi", "#phi^{t} gen", 30, -3.14, 3.14);

  book<TH1F>("N_gluon", "N_{g} gen", 5, 0, 5);
  book<TH1F>("gluon_pt", "p_{T}^{g} gen", 50, 0, 2000);
  book<TH1F>("gluon_eta", "#eta^{g} gen", 52, -5.2, 5.2);
  book<TH1F>("gluon_phi", "#phi^{g} gen", 30, -3.14, 3.14);

  book<TH1F>("N_electron", "N_{e} gen", 5, 0, 5);
  book<TH1F>("electron_pt", "p_{T}^{e} gen", 50, 0, 2000);
  book<TH1F>("electron_eta", "#eta^{e} gen", 52, -5.2, 5.2);
  book<TH1F>("electron_phi", "#phi^{e} gen", 30, -3.14, 3.14);

  book<TH1F>("N_muon", "N_{#mu} gen", 5, 0, 5);
  book<TH1F>("muon_pt", "p_{T}^{#mu} gen", 50, 0, 2000);
  book<TH1F>("muon_eta", "#eta^{#mu} gen", 52, -5.2, 5.2);
  book<TH1F>("muon_phi", "#phi^{#mu} gen", 30, -3.14, 3.14);

  book<TH1F>("N_b", "N_{b} gen", 5, 0, 5);
  book<TH1F>("b_pt", "p_{T}^{b} gen", 50, 0, 2000);
  book<TH1F>("b_eta", "#eta^{b} gen", 52, -5.2, 5.2);
  book<TH1F>("b_phi", "#phi^{b} gen", 30, -3.14, 3.14);

  book<TH1F>("N_q", "N_{q} gen", 5, 0, 5);
  book<TH1F>("q_pt", "p_{T}^{q} gen", 50, 0, 2000);
  book<TH1F>("q_eta", "#eta^{q} gen", 52, -5.2, 5.2);
  book<TH1F>("q_phi", "#phi^{q} gen", 30, -3.14, 3.14);

  // #################
  // ## Other hists ##
  // #################

  book<TH1F>("dR_tops", "#DeltaR_{T* T* or t T*} gen", 30, 0, 6);
  book<TH1F>("deta_tops", "#Delta#eta_{T* T* or t T*} gen", 30, 0, 6);
  book<TH1F>("dphi_tops", "#Delta#phi_{T* T* or t T*} gen", 30, 0, 6);

  book<TH1F>("invmass_top", "m_{t t} gen", 100, 0, 2000);
  book<TH1F>("invmass_gluon", "m_{g g} gen", 100, 0, 2000);

  book<TH1F>("dR_gluons", "#DeltaR_{g g} gen", 30, 0, 6);
  book<TH1F>("deta_gluons", "#Delta#eta_{g g} gen", 30, 0, 6);
  book<TH1F>("dphi_gluons", "#Delta#phi_{g g} gen", 30, 0, 6);

  book<TH1F>("dR_gluon_top", "#DeltaR_{g t} gen", 30, 0, 6);
  book<TH1F>("deta_gluon_top", "#Delta#eta_{g t} gen", 30, 0, 6);
  book<TH1F>("dphi_gluon_top", "#Delta#phi_{g t} gen", 30, 0, 6);

  book<TH1F>("dR_qs", "#DeltaR_{q q} gen", 30, 0, 6);
  book<TH1F>("deta_qs", "#Delta#eta_{q q} gen", 30, 0, 6);
  book<TH1F>("dphi_qs", "#Delta#phi_{q q} gen", 30, 0, 6);

  book<TH1F>("recoTstar_mass", "m_{t g} gen", 100, 0, 2000);
  book<TH1F>("recoTstar_pt", "p_{T}^{t g} gen", 50, 0, 2000);

  is_mc = ctx.get("dataset_type") == "MC";

}


void TstarTstarModelHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  hist("weight")->Fill(weight);

  // finding Tstars and determining which model is used
  std::vector<GenParticle> Tstars;
  bool isSpin12 = false;
  for(const GenParticle & gp : *event.genparticles){
    if((gp.pdgId() == 25001) || (gp.pdgId() == -25001)) {
      isSpin12 = true;
      Tstars.push_back(gp);
    }
    else if((gp.pdgId() == 9000005 && (gp.status()==23 || gp.status()==22)) || (gp.pdgId() == -9000005 && (gp.status()==23 || gp.status()==22))) {
      if(isSpin12) std::cout << "Error, models mixed?" << std::endl;
      Tstars.push_back(gp);
    }
  }

  // N_Tstar to determine whether pair production or single production
  int N_Tstar = Tstars.size();

  // filling, depending on used model
  if(isSpin12) {
    std::vector<GenParticle> tops, gluons, electrons, muons, bs, qs;
    for(const GenParticle & gp : *event.genparticles){
      if((gp.pdgId() == 6) || (gp.pdgId() == -6)) tops.push_back(gp);
      else if((gp.pdgId() == 21) && (gp.status() > 0)) gluons.push_back(gp);
      else if((gp.pdgId() == 11) || (gp.pdgId() == -11)) electrons.push_back(gp);
      else if((gp.pdgId() == 13) || (gp.pdgId() == -13)) muons.push_back(gp);
      else if((gp.pdgId() == 5) || (gp.pdgId() == -5)) bs.push_back(gp);
      else if(((abs(gp.pdgId()) == 1) || (abs(gp.pdgId()) == 2) || (abs(gp.pdgId()) == 3) || (abs(gp.pdgId()) == 4)) && (gp.status() > 0)) qs.push_back(gp);
    }
  }

  hist("N_Tstar")->Fill(N_Tstar, weight);
  for (const auto & Tstar : Tstars){
    hist("Tstar_spin")->Fill(Tstar.spin(), weight/N_Tstar);
    hist("Tstar_mass")->Fill(inv_mass(Tstar.v4()), weight/N_Tstar);
    hist("Tstar_pt")->Fill(Tstar.pt(), weight/N_Tstar);
    hist("Tstar_eta")->Fill(Tstar.eta(), weight/N_Tstar);
    hist("Tstar_phi")->Fill(Tstar.phi(), weight/N_Tstar);
  }

  int N_top = tops.size();
  hist("N_top")->Fill(N_top, weight);
  for (const auto & top : tops){
    hist("top_spin")->Fill(top.spin(), weight/N_top);
    hist("top_mass")->Fill(inv_mass(top.v4()), weight/N_top);
    hist("top_pt")->Fill(top.pt(), weight/N_top);
    hist("top_eta")->Fill(top.eta(), weight/N_top);
    hist("top_phi")->Fill(top.phi(), weight/N_top);
  }

  // angular plots
  if(N_Tstar == 2){
    hist("dR_tops")->Fill(deltaR(Tstars.at(0), Tstars.at(1)), weight);
    hist("deta_tops")->Fill(abs(Tstars.at(0).eta() - Tstars.at(1).eta()), weight);
    hist("dphi_tops")->Fill(abs(Tstars.at(0).phi() - Tstars.at(1).phi()), weight);
  }
  else {
    hist("dR_tops")->Fill(deltaR(Tstars.at(0), tops.at(0)), weight);
    hist("deta_tops")->Fill(abs(Tstars.at(0).eta() - tops.at(0).eta()), weight);
    hist("dphi_tops")->Fill(abs(Tstars.at(0).phi() - tops.at(0).phi()), weight);
  }

  int N_gluon = gluons.size();
  hist("N_gluon")->Fill(N_gluon, weight);
  int count = 0; // crosscheck, should in the end be the same as N_Tstar
  for (const auto & gluon : gluons){
    hist("gluon_pt")->Fill(gluon.pt(), weight/N_gluon);
    hist("gluon_eta")->Fill(gluon.eta(), weight/N_gluon);
    hist("gluon_phi")->Fill(gluon.phi(), weight/N_gluon);
    for (const auto & top : tops){
      if(top.mother1() == gluon.mother1()) {
        hist("dR_gluon_top")->Fill(deltaR(gluon, top), weight/N_Tstar);
        hist("deta_gluon_top")->Fill(abs(gluon.eta() - top.eta()), weight/N_Tstar);
        hist("dphi_gluon_top")->Fill(abs(gluon.phi() - top.phi()), weight/N_Tstar);
        hist("recoTstar_pt")->Fill((gluon.v4()+top.v4()).pt(), weight/N_Tstar);
        hist("recoTstar_mass")->Fill(inv_mass(gluon.v4()+top.v4()), weight/N_Tstar);
        hist("top_pt_fromTstar")->Fill(top.pt(), weight/N_top);
        count++;
      }
      else if(N_Tstar == 1) hist("top_pt_initial")->Fill(top.pt(), weight/N_top); // TODO this is an ugly hack. Fix it Finn!
    }
  }
  if(count > N_Tstar) std::cout << "Found too many top-gluon pairs." << std::endl;

  if(N_gluon == 2){
    hist("dR_gluons")->Fill(deltaR(gluons.at(0), gluons.at(1)), weight);
    hist("deta_gluons")->Fill(abs(gluons.at(0).eta() - gluons.at(1).eta()), weight);
    hist("dphi_gluons")->Fill(abs(gluons.at(0).phi() - gluons.at(1).phi()), weight);
  }

  // invariant masses
  // currently only if exactly two tops / gluons are present (which should be the case in the hard process)
  if(N_top == 2) hist("invmass_top")->Fill(inv_mass(tops.at(0).v4() + tops.at(1).v4()), weight);
  if(N_gluon == 2) hist("invmass_gluon")->Fill(inv_mass(gluons.at(0).v4() + gluons.at(1).v4()), weight);

  int N_electron = electrons.size();
  hist("N_electron")->Fill(N_electron, weight);
  for (const auto & electron : electrons){
    hist("electron_pt")->Fill(electron.pt(), weight/N_electron);
    hist("electron_eta")->Fill(electron.eta(), weight/N_electron);
    hist("electron_phi")->Fill(electron.phi(), weight/N_electron);
  }

  int N_muon = muons.size();
  hist("N_muon")->Fill(N_muon, weight);
  for (const auto & muon : muons){
    hist("muon_pt")->Fill(muon.pt(), weight/N_muon);
    hist("muon_eta")->Fill(muon.eta(), weight/N_muon);
    hist("muon_phi")->Fill(muon.phi(), weight/N_muon);
  }

  int N_b = bs.size();
  hist("N_b")->Fill(N_b, weight);
  for (const auto & b : bs){
    hist("b_pt")->Fill(b.pt(), weight/N_b);
    hist("b_eta")->Fill(b.eta(), weight/N_b);
    hist("b_phi")->Fill(b.phi(), weight/N_b);
  }

  int N_q = qs.size();
  hist("N_q")->Fill(N_q, weight);
  for (const auto & q : qs){
    hist("q_pt")->Fill(q.pt(), weight/N_q);
    hist("q_eta")->Fill(q.eta(), weight/N_q);
    hist("q_phi")->Fill(q.phi(), weight/N_q);
  }

  if(N_q == 2){
    hist("dR_qs")->Fill(deltaR(qs.at(0), qs.at(1)), weight);
    hist("deta_qs")->Fill(abs(qs.at(0).eta() - qs.at(1).eta()), weight);
    hist("dphi_qs")->Fill(abs(qs.at(0).phi() - qs.at(1).phi()), weight);
  }

}

TstarTstarModelHists::~TstarTstarModelHists(){}
