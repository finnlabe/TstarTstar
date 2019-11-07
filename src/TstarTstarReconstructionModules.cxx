#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
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

ttbarChi2Discriminator::ttbarChi2Discriminator(uhh2::Context& ctx){

  //Getting needed Handles
  h_ttbar_hyps_ = ctx.get_handle< std::vector<ReconstructionHypothesis> >("TTbarReconstruction");
  h_is_ttbar_reconstructed_ = ctx.get_handle< bool >("is_ttbar_reconstructed_chi2");
  h_recohyp_ = ctx.declare_event_output<ReconstructionHypothesis>("TTbarReconstruction_best");

  // mtophad_ = 175.;
  // mtophad_ttag_ = 177.;
  // sigmatophad_ = 20.;
  // sigmatophad_ttag_ = 17.;
  // mtoplep_ = 173.;
  // mtoplep_ttag_ = 173.;
  // sigmatoplep_ = 29.;
  // sigmatoplep_ttag_ = 29.;

  mtoplep_ = 175.;
  sigmatoplep_ = 19.;
  mtophad_ = 177.;
  sigmatophad_ = 16.;

  mtoplep_ttag_ = 175.;
  sigmatoplep_ttag_ = 19.;
  mtophad_ttag_ = 173.;
  sigmatophad_ttag_ = 15.;

}

bool ttbarChi2Discriminator::process(uhh2::Event& event){
  
  //bool debug = true;
  bool debug = false;
  if(debug){cout << "Hello World from ttbarChi2Discriminator!" << endl;}  
 
  if(debug){cout << "Reading in the hypothesis vector... ";}
  vector<ReconstructionHypothesis>& candidates = event.get(h_ttbar_hyps_);
  ReconstructionHypothesis bestCand;

  if(debug && (candidates.size() < 1)){cout << "No candidates found!" << endl;}
  if(candidates.size() < 1) {
    bestCand = ReconstructionHypothesis(); //Fixing crash when no BestCandidate is found (handle must not be empty)
    event.set(h_recohyp_, bestCand);
    event.set(h_is_ttbar_reconstructed_, false);
    return false;
  }

  if(debug){cout << "Starting to find best candidate." << endl;}
  float chi2_best = 99999999;
  bestCand = candidates.at(0);
  for(unsigned int i=0; i<candidates.size(); i++){
    
    float chi2_had = 0.;
    float chi2_lep = 0.;
    float mhad = 0.;
    float mlep = 0.;

    /** TODO: What does this do? Do I need it?
    //bool is_toptag_reconstruction = candidates.at(i).is_toptag_reconstruction(); //Was macht das?
    if(is_toptag_reconstruction){
      if(!candidates.at(i).is_puppi_reconstruction()) mhad = candidates.at(i).tophad_topjet_ptr()->softdropmass();
      else{
        LorentzVector SumSubjets(0.,0.,0.,0.);
        for(unsigned int k=0; k<candidates.at(i).tophad_topjet_ptr()->subjets().size(); k++) SumSubjets = SumSubjets + candidates.at(i).tophad_topjet_ptr()->subjets().at(k).v4();
        mhad = inv_mass(SumSubjets);
      }
      mlep = inv_mass(candidates.at(i).top_leptonic_v4());
      chi2_had = pow((mhad - mtophad_ttag_) / sigmatophad_ttag_,2);
      chi2_lep = pow((mlep - mtoplep_ttag_) / sigmatoplep_ttag_,2);
    }
    else{
      mhad = inv_mass(candidates.at(i).top_hadronic_v4());
      mlep = inv_mass(candidates.at(i).top_leptonic_v4());
      chi2_had = pow((mhad - mtophad_) / sigmatophad_,2);
      chi2_lep = pow((mlep - mtoplep_) / sigmatoplep_,2);
    }
    **/

    mhad = inv_mass(candidates.at(i).tophad_v4());
    mlep = inv_mass(candidates.at(i).toplep_v4());
    chi2_had = pow((mhad - mtophad_) / sigmatophad_,2);
    chi2_lep = pow((mlep - mtoplep_) / sigmatoplep_,2);

    float chi2 = chi2_had + chi2_lep;
    
    /** TODO: make it possible to save chi2 in candidate as discriminator (for plotting, later use etc.)
    candidates.at(i).set_discriminators("chi2_hadronic", chi2_had);
    candidates.at(i).set_discriminators("chi2_leptonic", chi2_lep);
    candidates.at(i).set_discriminators("chi2_total", chi2);
    **/

    if(chi2 < chi2_best){
      if(debug){cout << "Best chi2 was " << chi2_best << ". Set new best chi2." << endl;}
      chi2_best = chi2;
      bestCand = candidates.at(i);
    }
  }

  //Writing best candidate into event
  if(debug){cout << "Finished finding best candidate!" << endl;} 
  event.set(h_recohyp_, bestCand);
  event.set(h_is_ttbar_reconstructed_, true);
  if(debug){cout << "Finished writing best candidate!" << endl;}
  return true;
}










// ##### Now for the ugly code:


TstarTstar_Reconstruction::TstarTstar_Reconstruction(uhh2::Context& ctx){
  // A lot of handles. Reconstructed Masses as well as several DeltaR are saved for plotting reasons mainly.
  // TODO later: better save reconstructed Tstar as some kind of object (particle?), to not use so many handles. Then DeltaR can be calculated when needed for plotting.
  h_recohyp_ = ctx.declare_event_output<ReconstructionHypothesis>("TTbarReconstruction_best");
  h_M_Tstar_gluon_ = ctx.get_handle< float >("M_Tstar_gluon");
  h_M_Tstar_gamma_ = ctx.get_handle< float >("M_Tstar_gamma");
  
  h_DeltaR_toplep_ak8jet1_ = ctx.get_handle< float >("DeltaR_toplep_ak8jet1");
  h_DeltaR_tophad_ak8jet1_ = ctx.get_handle< float >("DeltaR_tophad_ak8jet1");
  h_DeltaR_toplep_ak8jet2_ = ctx.get_handle< float >("DeltaR_toplep_ak8jet2");
  h_DeltaR_tophad_ak8jet2_ = ctx.get_handle< float >("DeltaR_tophad_ak8jet2");
}

bool TstarTstar_Reconstruction::process(uhh2::Event& event){

  //bool debug = true;
  bool debug = false;

  if(debug){cout << "Hello World from Tstartstar_Reconstruction!" << endl;}

  // Get Information about jets that have been used in Reconstruction
  if(debug){cout << "Reading in Jets from best Hypothesis & Event" << endl;}
  ReconstructionHypothesis hyp = event.get(h_recohyp_);
  std::vector<Jet> used_jets_ = hyp.tophad_jets();
  used_jets_.insert( used_jets_.end(), hyp.toplep_jets().begin(), hyp.toplep_jets().end() );
  
  std::vector<TopJet> all_topjets_ = *event.topjets;

  if(debug){cout << "Searching for separated AK8Jet... " << endl;}
  
  //Calculate DeltaR between ak8 and closest ttbarjet for finding not used leading AK8 Jet Candidate
  float DeltaR_min = 9999;
  Jet leadingAK8;
  bool leadingExists = false;

  for(uint j = 0; j < all_topjets_.size(); j++){
    for (uint i = 0; i < hyp.tophad_jets().size(); i++){
      float DeltaR_new = deltaR(hyp.tophad_jets().at(i), all_topjets_.at(j)); 
      if(DeltaR_new < DeltaR_min){
	DeltaR_min = DeltaR_new;
      }
    }
    for (uint i = 0; i < hyp.toplep_jets().size(); i++){
      float DeltaR_new = deltaR(hyp.toplep_jets().at(i), all_topjets_.at(j)); 
      if(DeltaR_new < DeltaR_min){
	DeltaR_min = DeltaR_new;
      }
    }
    //TODO: Maybe make cut applicable from outside
    if(/*DeltaR_min > 0.4*/ true){
      leadingAK8 = all_topjets_.at(j);
      leadingExists = true;
      break;
    }
  }

  if(!leadingExists){
    if(debug){cout << "Error! No separated AK8Jet found!" << endl;}
    return false;
  }
  
  std::vector<std::vector<float>> M_Tstar_combinations_;
  std::vector< float > M_Tstargluon_had_;
  std::vector< float > M_Tstargluon_lep_;
  std::vector< float > tempvec_;
  
  if(debug){cout << "Calculating all possible M_Tstar pairs..." << endl;}
  //Calculate all possible M_Tstargluon und M_Tstargamma
  //TODO (I dont like how this is done at the moment. Its fast enough at the moment, but maybe it can still be improved)
  M_Tstargluon_had_.push_back(inv_mass(hyp.tophad_v4() + leadingAK8.v4()));
  M_Tstargluon_lep_.push_back(inv_mass(hyp.toplep_v4() + leadingAK8.v4()));

  for (uint j = 0; j < M_Tstargluon_lep_.size(); j++){
    tempvec_.push_back(inv_mass(hyp.tophad_v4() + event.photons->at(0).v4()));
    tempvec_.push_back(M_Tstargluon_lep_.at(j));
    M_Tstar_combinations_.push_back(tempvec_);
  }
  for (uint j = 0; j < M_Tstargluon_had_.size(); j++){
    tempvec_.push_back(inv_mass(hyp.toplep_v4() + event.photons->at(0).v4()));
    tempvec_.push_back(M_Tstargluon_had_.at(j));
    M_Tstar_combinations_.push_back(tempvec_);
  }

  //Find best TstarTstar pair
  if(debug){cout << "Finding the best M_Tstar pair..." << endl;}
  float chi2;
  float chi2_best = 9999999;

  float M_Tstargluon_best_ = 0;
  float M_Tstargamma_best_ = 0;

  for (uint i = 0; i < M_Tstar_combinations_.size(); i++){
    tempvec_ = M_Tstar_combinations_.at(i);
    chi2 = pow(tempvec_.at(0) - tempvec_.at(1),2);
    
    if(chi2 < chi2_best){
      chi2_best = chi2;
      M_Tstargluon_best_ = tempvec_.at(0);
      M_Tstargamma_best_ = tempvec_.at(1);
    }
  }
  
  if(debug){cout << "Best mass pair is:    ( " << M_Tstargluon_best_ << " | " << M_Tstargamma_best_ << " )" << endl;}

  event.set(h_M_Tstar_gluon_, M_Tstargluon_best_);
  event.set(h_M_Tstar_gamma_, M_Tstargamma_best_);

  if(debug){cout << "Masses written, return to main" << endl;}
  
  //Still need to save what the "primary" photon and what the gluon jet is!

  return true;
}
