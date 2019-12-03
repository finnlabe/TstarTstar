#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>

using namespace std;
using namespace uhh2;

namespace {
    
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

}

TstarTstarGenMatcher::TstarTstarGenMatcher(uhh2::Context& ctx){

  h_tstartstar_hyps = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TTbarReconstruction");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;

}


bool TstarTstarGenMatcher::process(uhh2::Event& event){
  
  bool debug = false;

  // Define a correct ReconstructionTstarHypothesis
  // Will contain the best ttbar hypothesis (as those are already found by ttbar reco class), and this is only for matchable ones!
  // then manually insert correct gluon ak8 jet(s) and correct photon (for tgtgamma)
  ReconstructionTstarHypothesis correctTstarHyp = ReconstructionTstarHypothesis();

  // ####### Start of TstarTstar Reco stuff ######
  // TstarTstar GENReco
  TTbarGen ttbargen = event.get(h_ttbargen);

  // Define Reco stuff!
  LorentzVector reco_toplep_v4, reco_tophad_v4;
  Electron reco_ele;
  Muon reco_muon;
  bool isEle = false;
  bool isMu = false;

  if(is_tgtgamma){
    TopJet reco_gluon;
    Photon reco_pho;
  }
  if(is_tgtg){
    TopJet reco_gluon1;
    TopJet reco_gluon2;
  }

  // evaluate ttbar hypothesises
  for(const auto & hyp : hyps){
    
    ReconstructionHypothesis ttbarhyp = hyp.ttbar_hyp();

    // evaluate toplep
    double deltaR_toplep = 1e9;
    
    if(debug && (ttbarhyp.toplep_jets().size() > 1)){cout << "More than one jet used for leptonic top..." << endl;};
    for(const auto & jet : ttbarhyp.toplep_jets()){
      deltaR_toplep = deltaR(jet, ttbargen.BLep());
    }
    deltaR_toplep += deltaR(ttbarhyp.lepton(), ttbargen.ChargedLepton());
    
    // evaluate tophad
    double deltaR_tophad = 0.;
    // match to b.
    for(const auto & jet : ttbarhyp.tophad_jets()){
      deltaR_tophad += min( {deltaR(jet, ttbargen.BHad()), deltaR(jet, ttbargen.Q1()), deltaR(jet, ttbargen.Q2())} );
    }
    
    hyps_eval.push_back(deltaR_toplep + deltaR_tophad);
  }
  
  // find ttbar hypothesis with minimum deviation from GEN
  double deltaR_min = 1e9;
  int correctHyp = -1;
  for(unsigned int i = 0; i < hyps_eval.size(); i++){
    if(hyps_eval.at(i) < deltaR_min){
      deltaR_min = hyps_eval.at(i);
      correctHyp = i;
    }
  }

  correctTstarHyp.set_ttbar(hyps.at(correctHyp)); // Finished finding correct ttbar hypothesis!
  
  // match photon and gluon
  GenParticle gluon1, gluon2, photon1, photon2;
  bool found_gluon1 = false,  found_photon1 = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      if(!found_gluon1){
	gluon1 = gp;
	found_gluon1 = true;
      }
      else{
	gluon2 = gp;
      }
    }
    else if(gp.pdgId() == 22 && gp.status()==23){//only photons from Tstar decay
      if(!found_photon1){
	photon1 = gp;
	found_photon1 = true;
      }
      else{
	photon2 = gp;
      }
    }
  }
  
  // Find signatures matching GEN gluons(s) and photon (for tgtgamma)
  if(is_tgtgamma){
    
    double deltaR_gamma = 1e9;
    for(const auto & photon : *event.photons){
      double deltaR_tmp = deltaR(photon, photon1);
      if(deltaR_tmp < deltaR_gamma){
	deltaR_gamma = deltaR_tmp;
	
	// find out which top this belongs to.
	// for this, get parent of GEN photon, check whether it is 

      }
    }
    
    double deltaR_gluon1 = 1e9;
    for(const auto & jet : *event.topjets){
      double deltaR_tmp = deltaR(jet, gluon1);
      if(deltaR_tmp < deltaR_gluon1){
	deltaR_gluon1 = deltaR_tmp;
	//TODO
      }
    }
    
    
  }
  if(is_tgtg){
    
    double deltaR_gluon1 = 1e9;
    for(const auto & jet : *event.topjets){
      double deltaR_tmp = deltaR(jet, gluon1);
      if(deltaR_tmp < deltaR_gluon1){
	deltaR_gluon1 = deltaR_tmp;
	//TODO
      }
    }
    
    double deltaR_gluon2 = 1e9;
    for(const auto & jet : *event.topjets){
      double deltaR_tmp = deltaR(jet, gluon2);
      if(deltaR_tmp < deltaR_gluon2){
	deltaR_gluon2 = deltaR_tmp;
	//TODO
      }
    }
  }

  // TODO Plot!

  return true;
}
