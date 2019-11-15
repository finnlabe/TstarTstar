#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
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

  // mtoplep_ = 175.;
  // sigmatoplep_ = 19.;
  // mtophad_ = 177.;
  // sigmatophad_ = 16.;

  // mtoplep_ttag_ = 175.;
  // sigmatoplep_ttag_ = 19.;
  // mtophad_ttag_ = 173.;
  // sigmatophad_ttag_ = 15.;

  mtoplep_ = 166.;
  sigmatoplep_ = 22.;
  mtophad_ = 166.;
  sigmatophad_ = 16.;

  mtoplep_ttag_ = 166.;
  sigmatoplep_ttag_ = 22.;
  mtophad_ttag_ = 166.;
  sigmatophad_ttag_ = 16.;

}

bool ttbarChi2Discriminator::process(uhh2::Event& event){
  
  //  bool debug = true;
  bool debug = false;
  if(debug){cout << "Hello World from ttbarChi2Discriminator!" << endl;}  
 
  if(debug){cout << "Reading in the hypothesis vector... ";}
  const std::vector<ReconstructionHypothesis>& candidates = event.get(h_ttbar_hyps_);
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
    //    if(debug) cout<<"Mhad and Mlep candidates: "<<mhad<<" "<<mlep<<" chi2 = "<<chi2<<endl;
    
    /** TODO: make it possible to save chi2 in candidate as discriminator (for plotting, later use etc.)
    candidates.at(i).set_discriminators("chi2_hadronic", chi2_had);
    candidates.at(i).set_discriminators("chi2_leptonic", chi2_lep);
    candidates.at(i).set_discriminators("chi2_total", chi2);
    **/

    if(chi2 < chi2_best){
      if(debug){cout << "Best chi2 was " << chi2_best << ". Set new best chi2." << endl;}
      chi2_best = chi2;
      bestCand = candidates.at(i);
      bestCand.set_tophad_topjet_ptr(candidates.at(i).tophad_topjet_ptr());//TEST
    }
  }

  //Writing best candidate into event
  if(debug){cout << "Finished finding best candidate with chi2:" <<chi2_best <<" and masses had, lep "<<inv_mass(bestCand.tophad_v4())<<", "<<inv_mass(bestCand.toplep_v4())<< endl;} 
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
  // h_M_Tstar_gluon_ = ctx.get_handle< float >("M_Tstar_gluon");
  // h_M_Tstar_gamma_ = ctx.get_handle< float >("M_Tstar_gamma");
  
  // h_DeltaR_toplep_ak8jet1_ = ctx.get_handle< float >("DeltaR_toplep_ak8jet1");
  // h_DeltaR_tophad_ak8jet1_ = ctx.get_handle< float >("DeltaR_tophad_ak8jet1");
  // h_DeltaR_toplep_ak8jet2_ = ctx.get_handle< float >("DeltaR_toplep_ak8jet2");
  // h_DeltaR_tophad_ak8jet2_ = ctx.get_handle< float >("DeltaR_tophad_ak8jet2");
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

  // event.set(h_M_Tstar_gluon_, M_Tstargluon_best_);
  // event.set(h_M_Tstar_gamma_, M_Tstargamma_best_);

  if(debug){cout << "Masses written, return to main" << endl;}
  
  //Still need to save what the "primary" photon and what the gluon jet is!

  return true;
}

/////// T*T* -> t+gluon + t+gluon

TstarTstar_tgluon_tgluon_Reconstruction::TstarTstar_tgluon_tgluon_Reconstruction(uhh2::Context& ctx){
  h_recohyp_ttbar_ = ctx.get_handle<ReconstructionHypothesis>("TTbarReconstruction_best");
  h_recohyp_tstartstar_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtg_best");
}

bool TstarTstar_tgluon_tgluon_Reconstruction::process(uhh2::Event& event){

  //  bool debug = true;
  bool debug = false;

  if(debug){cout << "Hello World from Tstartstar_tgluon_tgluon_Reconstruction!" << endl;}

  // Get Information about jets that have been used in Reconstruction
  if(debug){cout << "Reading in Jets from best Hypothesis & Event" << endl;}
  ReconstructionHypothesis hyp_ttbar = event.get(h_recohyp_ttbar_);//best ttbar hypothesis
  //  cout<<"hyp_ttbar.tophad_topjet_ptr() = "<<hyp_ttbar.tophad_topjet_ptr()<<endl;
  //  std::vector<ReconstructionTstarHypothesis> recoHyps;
  ReconstructionTstarHypothesis recoHyp_best;
  recoHyp_best.set_ttbar(hyp_ttbar);
  //  event.set(h_recohyps_tstartstar_,recoHyps);//fill empty values
  event.set(h_recohyp_tstartstar_best_,recoHyp_best);//fill empty values

  std::vector<Jet> used_jets_ = hyp_ttbar.tophad_jets();
  used_jets_.insert( used_jets_.end(), hyp_ttbar.toplep_jets().begin(), hyp_ttbar.toplep_jets().end() );
  std::vector<TopJet> all_topjets_ = *event.topjets;
  std::vector<TopJet> notused_topjets_; 
  if(debug){cout << "Searching for separated AK8Jet... " << endl;}
  for(uint j = 0; j < all_topjets_.size(); j++){
    for (uint i = 0; i < used_jets_.size(); i++){
      double dR_ =  deltaR(all_topjets_.at(j),used_jets_.at(i));
      if(dR_<0.4) notused_topjets_.push_back(all_topjets_.at(j));
    }
  }

  if(debug) cout<<"Number of not used in ttbar AK8 jets is "<<notused_topjets_.size()<<endl;
  if(notused_topjets_.size()<2) return false;
  //  std::vector< float > M_Tstargluon_had_;
  //  std::vector< float > M_Tstargluon_lep_;
  double M_Tstargluon_had,M_Tstargluon_lep, diffM;
  std::pair<int,int> gluon_cand_min;//gluon from hadronic (first) and leptonic (second) side
  double diffM_min = 1e7;
  for (uint i1 = 0; i1 < notused_topjets_.size(); i1++){
    Jet gluon_cand1 = notused_topjets_.at(i1);
    M_Tstargluon_had = inv_mass(hyp_ttbar.tophad_v4() + gluon_cand1.v4());
    for (uint i2 = 0; i2 < notused_topjets_.size(); i2++){
      Jet gluon_cand2 = notused_topjets_.at(i2);
      if(i1!=i2){
	M_Tstargluon_lep = inv_mass(hyp_ttbar.toplep_v4() + gluon_cand2.v4());
	diffM = abs(M_Tstargluon_had-M_Tstargluon_lep);
      }
      else{
	M_Tstargluon_lep = inv_mass(hyp_ttbar.toplep_v4());
	diffM = 1e7;
      }
      if(diffM<diffM_min){
	diffM_min = diffM;
	gluon_cand_min.first = i1;
	gluon_cand_min.second = i2;
	ReconstructionTstarHypothesis current_tstartstar;
	current_tstartstar.set_ttbar(hyp_ttbar);

	current_tstartstar.set_tstarhad_v4(hyp_ttbar.tophad_v4() + gluon_cand1.v4());
	current_tstartstar.set_tstarlep_v4(hyp_ttbar.toplep_v4() + gluon_cand2.v4());

	current_tstartstar.set_tstar1gluon_v4(hyp_ttbar.toplep_v4() + gluon_cand2.v4());
	current_tstartstar.set_tstar2gluon_v4(hyp_ttbar.tophad_v4() + gluon_cand1.v4());

	current_tstartstar.add_tstarhad_jet(gluon_cand1);
	current_tstartstar.add_tstarlep_jet(gluon_cand2);
	//	recoHyps.push_back(current_tstartstar);
	recoHyp_best = current_tstartstar;
      }
      if(debug) cout<<"i1,i2:"<<i1<<" "<<i2<<" M_Tstargluon_had, M_Tstargluon_lep "<<M_Tstargluon_had<<", "<<M_Tstargluon_lep<<"; diff = "<<diffM<<endl;
    }
  }

  double M_Tstargluon_had_best_ = inv_mass(hyp_ttbar.tophad_v4() + notused_topjets_.at(gluon_cand_min.first).v4());
  double M_Tstargluon_lep_best_ = inv_mass(hyp_ttbar.toplep_v4() + notused_topjets_.at(gluon_cand_min.second).v4());
  if(debug) cout<<"ttbar: lep "<<inv_mass(hyp_ttbar.toplep_v4())<<" had:"<<inv_mass(hyp_ttbar.tophad_v4())<<endl;
  if(debug) cout<<"--- Jetid for hadronic and leptonic gluon: "<<gluon_cand_min.first<<" "<<gluon_cand_min.second<<" had,lep:"<<M_Tstargluon_had_best_<<", "<<M_Tstargluon_lep_best_<<" ---"<<endl;
  // event.set(h_M_Tstar_lep_, M_Tstargluon_lep_best_);
  // event.set(h_M_Tstar_had_, M_Tstargluon_had_best_);

  //  event.set(h_recohyps_tstartstar_,recoHyps);
  event.set(h_recohyp_tstartstar_best_,recoHyp_best);

  if(debug){cout << "Hypothesis written, return to main" << endl;}
  return true;
}
 
