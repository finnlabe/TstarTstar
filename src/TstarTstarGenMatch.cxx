#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include <UHH2/common/include/TTbarReconstruction.h>
#include <UHH2/common/include/TTbarGen.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"


using namespace std;
using namespace uhh2;

namespace {
    
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

}

TstarTstarGenMatcher::TstarTstarGenMatcher(uhh2::Context& ctx){
  
  is_tgtg = false; is_tgtgamma = false;
  if(ctx.get("channel") == "tgtg") is_tgtg = true;
  if(ctx.get("channel") == "tgtgamma") is_tgtgamma = true;

  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());

}


bool TstarTstarGenMatcher::process(uhh2::Event& event){
  
  bool debug = false;
  
  if(debug){cout << "Hello World from TstarTstarGenMatcher!" << endl;} 

  TTbarGen ttbargen = event.get(h_ttbargen);
  if(!ttbargen.IsSemiLeptonicDecay()){return false;}
  const bool pass_ttbarsemilep = TTbarSemiLepMatchable_selection->passes(event); // just to make sure this has been used.

  // ####### Start of TstarTstar Reco stuff ######
  // TstarTstar GENReco
  ReconstructionTstarHypothesis correctTstarHyp;
  std::vector<ReconstructionTstarHypothesis> tstartstar_all_hyps = event.get(h_tstartstar_hyp_vector); // get all ttbar hyps

  std::vector<double> tstartstsar_hyp_eval; // vector to save deltaRs in.

  // Finding GEN particles (other than ttbar which is in ttbargen)
  if(debug){cout << "Start finding GEN gluon(s) and photon." << endl;}
  GenParticle gluon1, gluon2, photon1, photon2, top1, top2;
  bool found_gluon1 = false,  found_photon1 = false, found_gluon2 = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 21 && gp.status()==23){ //only gluons from Tstar decay
      if(!found_gluon1){
	gluon1 = gp;
	found_gluon1 = true;
      }
      else{
	gluon2 = gp;
	found_gluon2 = true;
      }
    }
    else if(gp.pdgId() == 22 && gp.status()==23){ //only photons from Tstar decay
      if(!found_photon1){
	photon1 = gp;
	found_photon1 = true;
      }
      else{
	photon2 = gp;
      }
    }
  }
  if(is_tgtgamma && (!found_photon1 || !found_gluon1)){
    if(debug){cout << "Error: tgtgamma, but no photon or gluon found" << endl;}
    return false;
  }
  if(is_tgtg && (!found_gluon1 || !found_gluon2)){
    if(debug){cout << "Error: tgtg, but no photon or gluon found" << endl;}
    return false;
  } // Done finding GEN
  
  // Big LOOP over all TstarHypothesissesss
  double mindRsum = 1e6; int best_match_hyp = -1; // For finding best hyp
  for(unsigned int i=0; i<tstartstar_all_hyps.size(); i++){
    double deltaR_current;
 
    // Getting deltaR for ttbar
    std::pair<bool,double> pass_check_reco_ttbar_pair = TTbarSemiLepMatchable_selection->check_reco(tstartstar_all_hyps.at(i).ttbar_hyp());
    deltaR_current = pass_check_reco_ttbar_pair.second; // save deltaR value for ttbar.

    // Getting deltaR for gluon(s) / photon(s)
    if(is_tgtg){
      double deltaR_gluons = min( deltaR(gluon1.v4(), tstartstar_all_hyps.at(i).gluon1_v4()) + deltaR(gluon2.v4(), tstartstar_all_hyps.at(i).gluon2_v4())  , deltaR(gluon2.v4(), tstartstar_all_hyps.at(i).gluon1_v4()) + deltaR(gluon1.v4(), tstartstar_all_hyps.at(i).gluon1_v4()));
      deltaR_current += deltaR_gluons;
    }

    // Checking if new best
    if(deltaR_current < mindRsum){
      mindRsum = deltaR_current;
      best_match_hyp = i;
    }
  }
  
  if(debug){cout << "Done finding best ttbarhyp. Writing..." << endl;} 

  event.set(h_tstartstar_hyp, tstartstar_all_hyps.at(best_match_hyp));

  if(debug){cout << "Done writing, return to main!" << endl;}

  return true;
}
