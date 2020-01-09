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

bool TopTagMassWindow::operator()(const TopJet & topjet, const uhh2::Event &) const {
  //auto subjets = topjet.subjets();
  //if(subjets.size() < 3) return false;
  
  float mjet = topjet.softdropmass();
  
  if(mjet < m_mlower) return false;
  if(mjet > m_mupper) return false;
  
  return true;
}

TopTagMassWindow::TopTagMassWindow(double mlower, double mupper): m_mlower(mlower), m_mupper(mupper) {}

bool TopTagSubbtag::operator()(const TopJet & topjet, const uhh2::Event &) const {
  auto subjets = topjet.subjets();
  //if(subjets.size() < 3) return false;
  
  for(const auto & subjet : subjets){
    if(subjet.btag_DeepCSV() > m_btag){return true;}
  }
  
  return false;
}

TopTagSubbtag::TopTagSubbtag(double btag): m_btag(btag) {}



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
  
  bool debug = false;
  if(debug){cout << "Hello World from ttbarChi2Discriminator!" << endl;}  
 
  if(debug){cout << "Reading in the hypothesis vector... ";}
  std::vector<ReconstructionHypothesis>& candidates = event.get(h_ttbar_hyps_);
  ReconstructionHypothesis bestCand;

  if(candidates.size() < 1) {
    if(debug){cout << "No candidates found!" << endl;}
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

    /** TODO: What does this do? Do I need it? TODO this has to do with toptagging!!!
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
    
    candidates.at(i).set_discriminator("chi2_hadronic", chi2_had);
    candidates.at(i).set_discriminator("chi2_leptonic", chi2_lep);
    candidates.at(i).set_discriminator("chi2_total", chi2);
    candidates.at(i).set_discriminator("CorrectMatch", 0);

    if(chi2 < chi2_best){
      if(debug){cout << "Best chi2 was " << chi2_best << ". Set new best chi2." << endl;}
      chi2_best = chi2;
      bestCand = candidates.at(i);
      //      bestCand.set_tophad_topjet_ptr(candidates.at(i).tophad_topjet_ptr());//TEST
    }
  }

  //Writing best candidate into event
  if(debug){cout << "Finished finding best candidate with chi2:" <<chi2_best <<" and masses had, lep "<<inv_mass(bestCand.tophad_v4())<<", "<<inv_mass(bestCand.toplep_v4())<< endl;} 
  event.set(h_recohyp_, bestCand);
  event.set(h_is_ttbar_reconstructed_, true);
  if(debug){cout << "Finished writing best candidate!" << endl;}
  return true;
}

TstarTstar_tgluon_tgamma_Reconstruction::TstarTstar_tgluon_tgamma_Reconstruction(uhh2::Context& ctx){
  h_recohyp_ttbar_ = ctx.get_handle<ReconstructionHypothesis>("TTbarReconstruction_best");
  h_recohyp_tstartstar_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtgamma_best");
  
}

bool TstarTstar_tgluon_tgamma_Reconstruction::process(uhh2::Event& event){

  bool debug = false;
  
  if(debug){cout << "Hello World from Tstartstar_tgluon_tgamma_Reconstruction!" << endl;}

  // Get Information about jets that have been used in Reconstruction
  if(debug){cout << "Reading in Jets from best Hypothesis & Event" << endl;}
  ReconstructionHypothesis hyp_ttbar = event.get(h_recohyp_ttbar_);
  ReconstructionTstarHypothesis recoHyp_best;
  recoHyp_best.set_ttbar(hyp_ttbar);
  event.set(h_recohyp_tstartstar_best_, recoHyp_best);

  std::vector<Jet> used_jets_ = hyp_ttbar.tophad_jets();
  used_jets_.insert( used_jets_.end(), hyp_ttbar.toplep_jets().begin(), hyp_ttbar.toplep_jets().end() );
  std::vector<TopJet> all_topjets_ = *event.topjets;
  std::vector<bool> overlap_flag;

  for(uint j = 0; j < all_topjets_.size(); j++){
    overlap_flag.push_back(false);
  }

  if(debug){cout << "Searching for separated AK8Jet... " << endl;}
  for(uint j = 0; j < all_topjets_.size(); j++){
    for (uint i = 0; i < used_jets_.size(); i++){
      double dR_ =  deltaR(all_topjets_.at(j), used_jets_.at(i));
      if(dR_<1.2) overlap_flag[j] = true;
    }
  }

  int notused_topjets = 0;
  std::vector<TopJet> notused_topjets_;
  for(uint j = 0; j < all_topjets_.size(); j++){
    if(!overlap_flag[j]){
      notused_topjets_.push_back(all_topjets_.at(j));
      notused_topjets++;
    }
  }

  if(debug) cout<<"Number of not used in ttbar AK8 jets is "<<notused_topjets<<endl;
  if(notused_topjets<1) return false; // Require one AK8 Jet

  double M_tophad_gluon, M_toplep_gluon, M_tophad_photon, M_toplep_photon;
  float diffM_min = 1e7;

  if(debug){cout << "Calculating all possible M_Tstar pairs..." << endl;}

  for (uint i = 0; i < notused_topjets_.size(); i++){
    M_tophad_gluon = inv_mass(hyp_ttbar.tophad_v4() + notused_topjets_.at(i).v4());
    M_toplep_gluon = inv_mass(hyp_ttbar.toplep_v4() + notused_topjets_.at(i).v4());
    
    for (uint j = 0; j < event.photons->size(); j++){
      M_tophad_photon = inv_mass(hyp_ttbar.tophad_v4() + event.photons->at(j).v4());
      M_toplep_photon = inv_mass(hyp_ttbar.toplep_v4() + event.photons->at(j).v4());

      if(abs(M_tophad_gluon - M_toplep_photon) < diffM_min){ // hadronic - gluonic    ||    leptonic - photonic
	diffM_min = abs(M_tophad_gluon - M_toplep_photon);

	ReconstructionTstarHypothesis current_tstartstar;
	current_tstartstar.set_ttbar(hyp_ttbar);

	current_tstartstar.set_tstarhad_v4(hyp_ttbar.tophad_v4() + all_topjets_.at(i).v4());
	current_tstartstar.set_tstarlep_v4(hyp_ttbar.toplep_v4() + event.photons->at(j).v4());

	current_tstartstar.set_tstar1gluon_v4(hyp_ttbar.tophad_v4() + all_topjets_.at(i).v4());
	current_tstartstar.set_tstar1gamma_v4(hyp_ttbar.toplep_v4() + event.photons->at(j).v4());

	current_tstartstar.add_tstarhad_jet(all_topjets_.at(i));
	current_tstartstar.add_tstarlep_photon(event.photons->at(j));
	//	recoHyps.push_back(current_tstartstar);
	recoHyp_best = current_tstartstar;
      }
      if(abs(M_toplep_gluon - M_tophad_photon) < diffM_min){ // hadronic - photonic    ||    leptonic - gluonic
	diffM_min = abs(M_toplep_gluon - M_tophad_photon);

	ReconstructionTstarHypothesis current_tstartstar;
	current_tstartstar.set_ttbar(hyp_ttbar);

	current_tstartstar.set_tstarlep_v4(hyp_ttbar.toplep_v4() + all_topjets_.at(i).v4());
	current_tstartstar.set_tstarhad_v4(hyp_ttbar.tophad_v4() + event.photons->at(j).v4());

	current_tstartstar.set_tstar1gluon_v4(hyp_ttbar.toplep_v4() + all_topjets_.at(i).v4());
	current_tstartstar.set_tstar1gamma_v4(hyp_ttbar.tophad_v4() + event.photons->at(j).v4());

	current_tstartstar.add_tstarlep_jet(all_topjets_.at(i));
	current_tstartstar.add_tstarhad_photon(event.photons->at(j));
	//recoHyps.push_back(current_tstartstar);
	recoHyp_best = current_tstartstar;
      }
    }
  }

  event.set(h_recohyp_tstartstar_best_, recoHyp_best);

  if(debug){cout << "Hypothesis written, return to main" << endl;}
  return true;
}

/////// T*T* -> t+gluon + t+gluon

TstarTstar_tgluon_tgluon_Reconstruction::TstarTstar_tgluon_tgluon_Reconstruction(uhh2::Context& ctx){
  h_recohyp_ttbar_ = ctx.get_handle<ReconstructionHypothesis>("TTbarReconstruction_best");
  h_recohyp_tstartstar_best_ = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_tgtg_best");
  h_flag_toptagevent = ctx.get_handle<int>("flag_toptagevent");

}

bool TstarTstar_tgluon_tgluon_Reconstruction::process(uhh2::Event& event){

  bool debug = false;

  if(debug){cout << "Hello World from Tstartstar_tgluon_tgluon_Reconstruction!" << endl;}

  // Get Information about jets that have been used in Reconstruction
  if(debug){cout << "Reading in Jets from best Hypothesis & Event" << endl;}
  ReconstructionHypothesis hyp_ttbar = event.get(h_recohyp_ttbar_);//best ttbar hypothesis
  ReconstructionTstarHypothesis recoHyp_best;
  recoHyp_best.set_ttbar(hyp_ttbar);
  event.set(h_recohyp_tstartstar_best_,recoHyp_best);//fill empty values

  int notused_topjets = 0;
  std::vector<TopJet> notused_topjets_;
  std::vector<bool> overlap_flag;

  // do an if to see whether it it toptag!
  if(event.get(h_flag_toptagevent) == 1){
    TopJet toptagjet = *hyp_ttbar.tophad_topjet_ptr();

    std::vector<TopJet> all_topjets_ = *event.topjets;
    if(debug){cout << "We start with " << all_topjets_.size() << " topjets." << endl;}
    for(const auto & topjet : all_topjets_){
      if(toptagjet.v4() != topjet.v4()){notused_topjets_.push_back(topjet); notused_topjets++;}
    }

    if(all_topjets_.size()-1 != notused_topjets_.size()){cout << "Topjets remove from collection: " << (all_topjets_.size() - notused_topjets_.size()) << endl;}

  }
  else {
    std::vector<Jet> used_jets_ = hyp_ttbar.tophad_jets();
    used_jets_.insert( used_jets_.end(), hyp_ttbar.toplep_jets().begin(), hyp_ttbar.toplep_jets().end() );

    std::vector<TopJet> all_topjets_ = *event.topjets;

    for(uint j = 0; j < all_topjets_.size(); j++){
      overlap_flag.push_back(false);
    }
    if(debug){cout << "Searching for separated AK8Jet... " << endl;}
    for(uint j = 0; j < all_topjets_.size(); j++){
      for (uint i = 0; i < used_jets_.size(); i++){
	double dR_ =  deltaR(all_topjets_.at(j),used_jets_.at(i));
	if(dR_<1.2) overlap_flag[j] = true;
      }
    }
    for(uint j = 0; j < all_topjets_.size(); j++){
      if(!overlap_flag[j]){
	notused_topjets_.push_back(all_topjets_.at(j));
	notused_topjets++;
      }
    }
  }

  if(debug) cout<<"Number of not used in ttbar AK8 jets is "<<notused_topjets<<endl;
  if(notused_topjets<2) return false; // Require two AK8 Jets

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
	current_tstartstar.set_tstarhad_v4(hyp_ttbar.tophad_v4() + gluon_cand1.v4());
	current_tstartstar.set_tstarlep_v4(hyp_ttbar.toplep_v4() + gluon_cand2.v4());

	current_tstartstar.set_tstar1gluon_v4(hyp_ttbar.toplep_v4() + gluon_cand2.v4());
	current_tstartstar.set_tstar2gluon_v4(hyp_ttbar.tophad_v4() + gluon_cand1.v4());

	current_tstartstar.set_gluon1_v4(gluon_cand1.v4());
	current_tstartstar.set_gluon2_v4(gluon_cand2.v4());

	current_tstartstar.add_tstarhad_jet(gluon_cand1);
	current_tstartstar.add_tstarlep_jet(gluon_cand2);
	current_tstartstar.set_used_ttbarjet_flags(overlap_flag);
	//	recoHyps.push_back(current_tstartstar);
	current_tstartstar.set_ttbar(hyp_ttbar);
	recoHyp_best = current_tstartstar;
      }
      if(debug) cout<<"i1,i2:"<<i1<<" "<<i2<<" M_Tstargluon_had, M_Tstargluon_lep "<<M_Tstargluon_had<<", "<<M_Tstargluon_lep<<"; diff = "<<diffM<<endl;
    }
  }

  double M_Tstargluon_had_best_ = inv_mass(hyp_ttbar.tophad_v4() + notused_topjets_.at(gluon_cand_min.first).v4());
  double M_Tstargluon_lep_best_ = inv_mass(hyp_ttbar.toplep_v4() + notused_topjets_.at(gluon_cand_min.second).v4());
  if(debug) cout<<"ttbar: lep "<<inv_mass(hyp_ttbar.toplep_v4())<<" had:"<<inv_mass(hyp_ttbar.tophad_v4())<<" chi2 = "
		<<hyp_ttbar.discriminator("chi2_total")<<" chi2_lep = "<<hyp_ttbar.discriminator("chi2_leptonic")
		<<" chi2_had = "<<hyp_ttbar.discriminator("chi2_hadronic")<<" CorrMatch = "<<hyp_ttbar.discriminator("CorrectMatch")<<endl;
  if(debug) cout<<"  hyp.toplep_jets().size() = "<<hyp_ttbar.toplep_jets().size()<<" hyp.tophad_jets().size() ="<<hyp_ttbar.tophad_jets().size()<<endl;
  if(debug) cout<<"--- Jetid for hadronic and leptonic gluon: "<<gluon_cand_min.first<<" "<<gluon_cand_min.second<<" had,lep:"<<M_Tstargluon_had_best_<<", "<<M_Tstargluon_lep_best_<<" ---"<<endl;

  event.set(h_recohyp_tstartstar_best_,recoHyp_best);

  if(debug){cout << "Hypothesis written, return to main" << endl;}
  return true;
}
 

