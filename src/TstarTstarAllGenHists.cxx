#include "UHH2/TstarTstar/include/TstarTstarAllGenHists.h"
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
#include "Math/Boost.h"
#include "Math/VectorUtil.h"

using namespace std;
using namespace uhh2;
//using namespace uhh2examples;

namespace {
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }
}


TstarTstarAllGenHists::TstarTstarAllGenHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  book<TH1F>("dR_top_closest_gluon", "dR_top_closest_gluon", 24, 0, 6);
  book<TH1F>("dphi_top_closest_gluon", "dphi_top_closest_gluon", 24, 0, 8);
  book<TH1F>("deta_top_closest_gluon", "deta_top_closest_gluon", 24, 0, 8);

  book<TH1F>("matches_found", "matches_found", 10, 0, 10);
  book<TH1F>("pt_mother", "pt_mother", 20, 0, 2000);

  book<TH1F>("dR_top_matched_gluon", "dR_top_matched_gluon", 20, 0, 10);
  book<TH1F>("dphi_top_matched_gluon", "dphi_top_matched_gluon", 24, 0, 8);
  book<TH1F>("deta_top_matched_gluon", "deta_top_matched_gluon", 24, 0, 8);

  book<TH1F>("dR_top_matched_gluon_lowpt", "dR_top_matched_gluon_lowpt", 20, 0, 10);
  book<TH1F>("dphi_top_matched_gluon_lowpt", "dphi_top_matched_gluon_lowpt", 24, 0, 8);
  book<TH1F>("deta_top_matched_gluon_lowpt", "deta_top_matched_gluon_lowpt", 24, 0, 8);

  book<TH1F>("dR_top_matched_gluon_medpt", "dR_top_matched_gluon_medpt", 20, 0, 10);
  book<TH1F>("dphi_top_matched_gluon_medpt", "dphi_top_matched_gluon_medpt", 24, 0, 8);
  book<TH1F>("deta_top_matched_gluon_medpt", "deta_top_matched_gluon_medpt", 24, -0, 8);

  book<TH1F>("dR_top_matched_gluon_highpt", "dR_top_matched_gluon_highpt", 20, 0, 10);
  book<TH1F>("dphi_top_matched_gluon_highpt", "dphi_top_matched_gluon_highpt", 24, 0, 8);
  book<TH1F>("deta_top_matched_gluon_highpt", "deta_top_matched_gluon_highpt", 24, 0, 8);

  book<TH1F>("costheta", "costheta", 24, -1, 1);
  book<TH1F>("costheta_lowpt", "costheta_lowpt", 24, -1, 1);
  book<TH1F>("costheta_medpt", "costheta_medpt", 24, -1, 1);
  book<TH1F>("costheta_highpt", "costheta_highpt", 24, -1, 1);

  is_mc = ctx.get("dataset_type") == "MC";

}


void TstarTstarAllGenHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Find signal, top, and antitop
  GenParticle top, antitop, tstar, antitstar;
  std::vector<GenParticle> gluons;
  // Find tops
  bool found_top = false, found_antitop = false, found_tstar = false, found_antitstar = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 6){
      top = gp;
      if(found_top){
        std::cout << "Found multiple tops" << endl;
        return;
      }
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      if(found_antitop){
        std::cout << "Found multiple antitops" << endl;
        return;
      }
      found_antitop = true;
    }
    else if(gp.pdgId() == 9000005 && (gp.status()==23 || gp.status()==22)){
      tstar = gp;
      found_tstar = true;
    }
    else if(gp.pdgId() == -9000005 && (gp.status()==23 || gp.status()==22)){
      antitstar = gp;
      found_antitstar = true;
    }
    else if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      gluons.push_back(gp);
    }
  }

  if(!found_top || !found_antitop) return;

  //////////////////////////
  // The interesting part //
  //////////////////////////

  // matching
  std::vector<std::vector<GenParticle>> sets;
  for (const auto gluon : gluons) {
    std::vector<GenParticle> set;
    if (gluon.mother1() == top.mother1()) {
      set.push_back(event.genparticles->at(gluon.mother1()));
      set.push_back(top);
      set.push_back(gluon);
    }
    else if (event.genparticles->at(gluon.mother1()) == tstar) {
      set.push_back(tstar);
      set.push_back(top);
      set.push_back(gluon);
    }
    else if (gluon.mother1() == antitop.mother1()) {
      set.push_back(event.genparticles->at(gluon.mother1()));
      set.push_back(antitop);
      set.push_back(gluon);
    }
    else if (event.genparticles->at(gluon.mother1()) == antitstar) {
      set.push_back(antitstar);
      set.push_back(antitop);
      set.push_back(gluon);
    }
    if(set.size() == 3) sets.push_back(set);
    else if (set.size() != 0) { std::cout << "Set not size 3" << endl; return; }
  }

  int set_count = sets.size();
  hist("matches_found")->Fill(set_count, weight);

  for (const auto & set : sets) {
    // transform to rest frame of interaction
    double mother_pt = set.at(0).pt();
    LorentzVector mother_v4 = set.at(0).v4();
    if(mother_pt == 0) mother_pt = set.at(1).pt() + set.at(2).pt();
    LorentzVector top_v4 = set.at(1).v4();
    LorentzVector gluon_v4 = set.at(2).v4();

    auto boostvec_top = top_v4.BoostToCM(gluon_v4);
    auto boostvec_gluon = top_v4.BoostToCM(top_v4);
    ROOT::Math::Boost boost_top = ROOT::Math::Boost(boostvec_top);
    ROOT::Math::Boost boost_gluon = ROOT::Math::Boost(boostvec_gluon);
    top_v4 = boost_top(top_v4);
    gluon_v4 = boost_gluon(gluon_v4);

    double dR = deltaR(top_v4, gluon_v4);
    double deta = abs(top_v4.eta() - gluon_v4.eta());
    double dphi = abs(top_v4.phi() - gluon_v4.phi());
    double costheta = ROOT::Math::VectorUtil::CosTheta(mother_v4,gluon_v4);

    // fill total Hists
    hist("pt_mother")->Fill(mother_pt, weight/set_count);
    hist("dR_top_matched_gluon")->Fill(dR, weight/set_count);
    hist("deta_top_matched_gluon")->Fill(deta, weight/set_count);
    hist("dphi_top_matched_gluon")->Fill(dphi, weight/set_count);
    hist("costheta")->Fill(costheta, weight/set_count);

    // differentiate in mother pt
    if(mother_pt < 400){
      hist("dR_top_matched_gluon_lowpt")->Fill(dR, weight/set_count);
      hist("deta_top_matched_gluon_lowpt")->Fill(deta, weight/set_count);
      hist("dphi_top_matched_gluon_lowpt")->Fill(dphi, weight/set_count);
      hist("costheta_lowpt")->Fill(costheta, weight/set_count);
    }
    else if(mother_pt < 800){
      hist("dR_top_matched_gluon_medpt")->Fill(dR, weight/set_count);
      hist("deta_top_matched_gluon_medpt")->Fill(deta, weight/set_count);
      hist("dphi_top_matched_gluon_medpt")->Fill(dphi, weight/set_count);
      hist("costheta_medpt")->Fill(costheta, weight/set_count);
    }
    else if(mother_pt >= 800){
      hist("dR_top_matched_gluon_highpt")->Fill(dR, weight/set_count);
      hist("deta_top_matched_gluon_highpt")->Fill(deta, weight/set_count);
      hist("dphi_top_matched_gluon_highpt")->Fill(dphi, weight/set_count);
      hist("costheta_highpt")->Fill(costheta, weight/set_count);
    }
  }

  if(gluons.size() == 0) return;

  // find closest gluon
  GenParticle top_closest_gluon;
  double dR_top = 999;
  for(const auto gluon : gluons){
    double dR_step = deltaR(gluon, top);
    if(dR_step < dR_top) {
      dR_top = dR_step;
      top_closest_gluon = gluon;
    }
  }
  GenParticle antitop_closest_gluon;
  double dR_antitop = 999;
  for(const auto gluon : gluons){
    double dR_step = deltaR(gluon, antitop);
    if(dR_step < dR_antitop) {
      dR_antitop = dR_step;
      antitop_closest_gluon = gluon;
    }
  }
  hist("dR_top_closest_gluon")->Fill((dR_top + dR_antitop)/2, weight);
  double dphi_top = abs(top.phi() - top_closest_gluon.phi());
  double dphi_antitop = abs(antitop.phi() - antitop_closest_gluon.phi());
  hist("dphi_top_closest_gluon")->Fill((dphi_top + dphi_antitop)/2, weight);
  double deta_top = abs(top.eta() - top_closest_gluon.eta());
  double deta_antitop = abs(antitop.eta() - antitop_closest_gluon.eta());
  hist("deta_top_closest_gluon")->Fill((deta_top + deta_antitop)/2, weight);


}

TstarTstarAllGenHists::~TstarTstarAllGenHists(){}
