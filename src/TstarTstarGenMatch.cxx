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

  if(is_tgtg){
    const std::string tstartstar_hyps_label("TstarTstar_tgtg");
    h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
  }
  if(is_tgtgamma){   
    const std::string tstartstar_hyps_label("TstarTstar_tgtgamma");
    h_recohyp_tstartstar = ctx.declare_event_output<ReconstructionTstarHypothesis>(tstartstar_hyps_label+"_best");
  }

  h_ttbar_hyps = ctx.get_handle<std::vector<ReconstructionHypothesis>>("TTbarReconstruction");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  TTbarSemiLepMatchable_selection.reset(new TTbarSemiLepMatchableSelection());

}


bool TstarTstarGenMatcher::process(uhh2::Event& event){
  
  bool debug = false;
  
  LorentzVector reco_topphotonic_v4, reco_topgluonic_v4;
  LorentzVector reco_toplep_v4, reco_tophad_v4;

  const bool pass_ttbarsemilep = TTbarSemiLepMatchable_selection->passes(event); // just to make sure this has been used.

  if(debug){cout << "Hello World from TstarTstarGenMatcher!" << endl;} 

  // ####### Start of TstarTstar Reco stuff ######
  // TstarTstar GENReco
  TTbarGen ttbargen = event.get(h_ttbargen);
  ReconstructionTstarHypothesis correctTstarHyp = ReconstructionTstarHypothesis();

  std::vector<ReconstructionHypothesis> ttbar_all_hyps = event.get(h_ttbar_hyps); // get all ttbar hyps
  double mindRsum = 1e6; int best_match_hyp = -1;
  bool pass_check_reco_ttbar = false;
  for(unsigned int i=0; i<ttbar_all_hyps.size(); i++){ // evaluate "correct" ttbar hyp
    std::pair<bool,double> pass_check_reco_ttbar_pair = TTbarSemiLepMatchable_selection->check_reco(ttbar_all_hyps.at(i));
    if(pass_check_reco_ttbar_pair.first && pass_check_reco_ttbar_pair.second<mindRsum){
      mindRsum = pass_check_reco_ttbar_pair.second; 
      best_match_hyp = i;
      pass_check_reco_ttbar = true;
    }
  }
  
  if(!pass_check_reco_ttbar){
    if(debug){cout << "no good enough ttbar hyp found. aborting." << endl;}
    return false;
  }

  if(debug){cout << "Done finding best ttbarhyp. Writing..." << endl;} 

  correctTstarHyp.set_ttbar(ttbar_all_hyps.at(best_match_hyp));

  reco_toplep_v4 = ttbar_all_hyps.at(best_match_hyp).toplep_v4();
  reco_tophad_v4 = ttbar_all_hyps.at(best_match_hyp).tophad_v4();

  // For tgtgamma
  TopJet reco_gluon;
  Photon reco_gamma;
  // For tgtg
  TopJet reco_gluon_lep;
  TopJet reco_gluon_had;
  
  if(debug){cout << "Start finding correct gluon(s) and photon." << endl;}
  // match photon and gluon
  GenParticle gluon1, gluon2, photon1, photon2, top1, top2;
  bool found_gluon1 = false,  found_photon1 = false, found_gluon2 = false;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 21 && gp.status()==23){//only gluons from Tstar decay
      if(!found_gluon1){
	gluon1 = gp;
	found_gluon1 = true;
      }
      else{
	gluon2 = gp;
	found_gluon2 = true;
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
  if(is_tgtgamma && (!found_photon1 || !found_gluon1)){
    if(debug){cout << "Error: tgtgamma, but no photon or gluon found" << endl;}
    return false;
  }
  if(is_tgtg && (!found_gluon1 || !found_gluon2)){
    if(debug){cout << "Error: tgtg, but no photon or gluon found" << endl;}
    return false;
  }

  // Find signatures matching GEN gluons(s) and photon (for tgtgamma)
  if(debug){cout << "Start match GEN to detector for gluon(s) and photon." << endl;}
  if(is_tgtgamma){
    
    if( ttbargen.TopLep().mother1() == photon1.mother1() ){
      reco_topphotonic_v4 = reco_toplep_v4;
      reco_topgluonic_v4 = reco_tophad_v4;
    }
    else if( ttbargen.TopLep().mother1() == gluon1.mother1() ){
      reco_topphotonic_v4 = reco_tophad_v4;
      reco_topgluonic_v4 = reco_toplep_v4;
    }
    else{cout << "Error! t neither has a photon nor gluon sibling!" << endl;}

    double deltaR_gamma = 1e9;
    for(const auto & photon : *event.photons){
      double deltaR_tmp = deltaR(photon, photon1);
      if(deltaR_tmp < deltaR_gamma){
	deltaR_gamma = deltaR_tmp;
	
	reco_gamma = photon;

      }
    }
    
    double deltaR_gluon1 = 1e9;
    for(const auto & jet : *event.topjets){
      double deltaR_tmp = deltaR(jet, gluon1);
      if(deltaR_tmp < deltaR_gluon1){
	deltaR_gluon1 = deltaR_tmp;
	
	reco_gluon = jet;

      }
    }
    
    
  } // end tgtgamma
  if(is_tgtg){
    
    double deltaR_gluon1 = 1e9;
    double deltaR_gluon2 = 1e9;
    for(const auto & jet : *event.topjets){
      double deltaR_tmp_1 = deltaR(jet, gluon1);
      double deltaR_tmp_2 = deltaR(jet, gluon2); 

      if ((deltaR_tmp_1 < deltaR_tmp_2) && (deltaR_tmp_1 < deltaR_gluon1)){

	if((gluon1.mother1()) == (ttbargen.TopLep().mother1())){
	  reco_gluon_lep = jet;
	  deltaR_gluon1 = deltaR_tmp_1;
	}
	else if((gluon1.mother1()) == (ttbargen.TopHad().mother1())){
	  reco_gluon_had = jet;
	  deltaR_gluon1 = deltaR_tmp_1;
	}
	else {
	  cout << "This should not happen. gluon1 has no Tstar as mother!" << endl;
	  return false;
	}

      }
      else if((deltaR_tmp_2 < deltaR_tmp_1) && (deltaR_tmp_2 < deltaR_gluon2)){
	deltaR_gluon2 = deltaR_tmp_2;

	if((gluon2.mother1()) == (ttbargen.TopLep().mother1())){
	  reco_gluon_lep = jet;
	  deltaR_gluon2 = deltaR_tmp_2;
	}
	else if((gluon2.mother1()) == (ttbargen.TopHad().mother1())){
	  reco_gluon_had = jet;
	  deltaR_gluon2 = deltaR_tmp_2;
	}
	else {
	  cout << "This should not happen. gluon2 has no Tstar as mother!" << endl;
	  return false;
	}

      }
    }
  }

  if(debug){cout << "Start writing correct TstarTstar Hyp to event." << endl;}

  if(is_tgtg){
    correctTstarHyp.set_tstarlep_v4(reco_toplep_v4+reco_gluon_lep.v4());
    correctTstarHyp.set_tstarhad_v4(reco_tophad_v4+reco_gluon_had.v4());
  }
  if(is_tgtgamma){
    correctTstarHyp.set_tstar1gamma_v4(reco_topphotonic_v4+reco_gamma.v4());
    correctTstarHyp.set_tstar1gluon_v4(reco_topgluonic_v4+reco_gluon.v4());

    if((photon1.mother1()) == (ttbargen.TopLep().mother1())){
      correctTstarHyp.set_tstarlep_v4(reco_toplep_v4+reco_gamma.v4());
      correctTstarHyp.set_tstarhad_v4(reco_tophad_v4+reco_gluon.v4());
    }
    else if((gluon1.mother1()) == (ttbargen.TopLep().mother1())){
      correctTstarHyp.set_tstarlep_v4(reco_toplep_v4+reco_gluon.v4());
      correctTstarHyp.set_tstarhad_v4(reco_tophad_v4+reco_gamma.v4());
    }
    else {
      cout << "This should not happen. Photon has no Tstar as mother!" << endl;
      return false;
    }
  }
  
  event.set(h_recohyp_tstartstar, correctTstarHyp);

  //cout << "hadronic Top mass is: " << inv_mass(reco_tophad_v4) << endl;
  //if(is_tgtg){cout << "hadronic gluon mass is: " << inv_mass(reco_gluon_had.v4()) << endl << endl;}

  if(debug){cout << "Done writing, return to main!" << endl;}

  return true;
}
