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
  book<TH1F>("dphi_top_closest_gluon", "dphi_top_closest_gluon", 24, -4, 4);
  book<TH1F>("deta_top_closest_gluon", "deta_top_closest_gluon", 24, -4, 4);

  book<TH1F>("dR_top_matched_gluon", "dR_top_matched_gluon", 24, 0, 6);
  book<TH1F>("dphi_top_matched_gluon", "dphi_top_matched_gluon", 24, -4, 4);
  book<TH1F>("deta_top_matched_gluon", "deta_top_matched_gluon", 24, -4, 4);

  is_mc = ctx.get("dataset_type") == "MC";

}


void TstarTstarAllGenHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  // Find signal, top, and antitop
  GenParticle top, antitop;
  std::vector<GenParticle> gluons;
  // Find tops
  bool found_top = false, found_antitop = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 6){
      top = gp;
      found_top = true;
    }
    else if(gp.pdgId() == -6){
      antitop = gp;
      found_antitop = true;
    }
    else if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      gluons.push_back(gp);
    }
  }

  if(!found_top || !found_antitop) return;
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

  // matching
  bool found_topgluon = false;
  bool found_antitopgluon = false;
  GenParticle topgluon, antitopgluon;
  for (const auto gluon : gluons) {
    if((event.genparticles->at(gluon.mother1()) == top) || (gluon.mother1()==top.mother1())) {
      topgluon = gluon;
      found_topgluon = true;
    }
    else if((event.genparticles->at(gluon.mother1()) == antitop) || (gluon.mother1()==antitop.mother1())){
      found_antitopgluon = true;
      antitopgluon = gluon;
    }
  }

  if(found_topgluon && found_antitopgluon) {
    double dR_top_gluon = deltaR(top, topgluon);
    double dR_antitop_gluon = deltaR(antitop, antitopgluon);
    hist("dR_top_matched_gluon")->Fill((dR_top_gluon + dR_antitop_gluon)/2, weight);
    double dphi_top_gluon = abs(top.phi() - topgluon.phi());
    double dphi_antitop_gluon = abs(antitop.phi() - antitopgluon.phi());
    hist("dphi_top_matched_gluon")->Fill((dphi_top_gluon + dphi_antitop_gluon)/2, weight);
    double deta_top_gluon = abs(top.eta() - topgluon.eta());
    double deta_antitop_gluon = abs(antitop.eta() - antitopgluon.eta());
    hist("deta_top_matched_gluon")->Fill((deta_top_gluon + deta_antitop_gluon)/2, weight);
  }
  else if(found_topgluon && !found_antitopgluon) {
    double dR_top_gluon = deltaR(top, topgluon);
    hist("dR_top_matched_gluon")->Fill(dR_top_gluon, weight);
    double dphi_top_gluon = abs(top.phi() - topgluon.phi());
    hist("dphi_top_matched_gluon")->Fill(dphi_top_gluon, weight);
    double deta_top_gluon = abs(top.eta() - topgluon.eta());
    hist("deta_top_matched_gluon")->Fill(deta_top_gluon, weight);
  }
  else if(!found_topgluon && found_antitopgluon) {
    double dR_antitop_gluon = deltaR(antitop, antitopgluon);
    hist("dR_top_matched_gluon")->Fill(dR_antitop_gluon, weight);
    double dphi_antitop_gluon = abs(antitop.phi() - antitopgluon.phi());
    hist("dphi_top_matched_gluon")->Fill(dphi_antitop_gluon, weight);
    double deta_antitop_gluon = abs(antitop.eta() - antitopgluon.eta());
    hist("deta_top_matched_gluon")->Fill(deta_antitop_gluon, weight);
  }
}

TstarTstarAllGenHists::~TstarTstarAllGenHists(){}
