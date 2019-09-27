#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/core/include/Event.h"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>
//using namespace uhh2examples;
using namespace uhh2;
using namespace std;

//From Z'->ttbar
bool uhh2::TwoDCut::passes(const uhh2::Event& event){

  assert(event.muons || event.electrons);
  //  if((event.muons->size()+event.electrons->size()) != 1) throw runtime_error("In TwoDCut: Event does not have exactly one muon and electron.");
  if((event.muons->size()+event.electrons->size()) <1) throw runtime_error("In TwoDCut: Event does not have at least one muon and electron.");

  bool pass = false;
  if(event.muons->size()){
    float drmin, ptrel;
    drmin = event.muons->at(0).get_tag(Muon::twodcut_dRmin);
    ptrel = event.muons->at(0).get_tag(Muon::twodcut_pTrel);
    if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
  }
  if(event.electrons->size()){
    float drmin, ptrel;
    drmin = event.electrons->at(0).get_tag(Electron::twodcut_dRmin);
    ptrel = event.electrons->at(0).get_tag(Electron::twodcut_pTrel);
    if((drmin > min_deltaR_) || (ptrel > min_pTrel_)) pass = true;
  }

  return pass;
}

//From Z'->ttbar
TTbarSemiLepMatchableSelection::TTbarSemiLepMatchableSelection(){}
bool TTbarSemiLepMatchableSelection::passes(const Event & event){
  // For pdgID mapping see https://twiki.cern.ch/twiki/bin/view/Main/PdgId
  if(event.isRealData) return false;
  assert(event.genparticles);

  //check, if one top decays had and one decays lep
  bool found_had = false, found_lep = false;

  //Loop over genparticles
  for(const auto & gp : *event.genparticles){
    //    cout<<" gp.pdgId() = "<<gp.pdgId()<<" status() = "<<gp.status()<<endl;
    //Get tops
    if(fabs(gp.pdgId()) == 6){

      //Get b and W
      auto b = gp.daughter(event.genparticles,1);
      auto W = gp.daughter(event.genparticles,2);
      if(fabs(W->pdgId()) == 5 && fabs(b->pdgId()) == 24){
        b = gp.daughter(event.genparticles,2);
        W = gp.daughter(event.genparticles,1);
      }
      //      cout<<" W->pdgId() = "<<W->pdgId()<<" b->pdgId() = "<<b->pdgId()<<endl;
      if(abs(W->pdgId()) != 24) {
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 24)) {
            W = &genp;
            break;
          }
        }
      }
      if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3 && abs(b->pdgId()) != 1) {//b,s,d
        for(unsigned int j = 0; j < event.genparticles->size(); ++j) {
          const GenParticle & genp = event.genparticles->at(j);
          auto m1 = genp.mother(event.genparticles, 1);
          auto m2 = genp.mother(event.genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == gp.index()) || (m2 && m2->index() == gp.index()));
          if(has_top_mother && (abs(genp.pdgId()) == 5 || abs(genp.pdgId()) == 3 || abs(genp.pdgId()) == 1)) {
            b = &genp;
            break;
          }
        }
      }
      //      if(!((fabs(b->pdgId()) == 5 || fabs(b->pdgId()) == 3 || fabs(b->pdgId()) == 1) && fabs(W->pdgId()) == 24)) cout<<"!!! Top decay products (W,b) are not found"<<endl;
      if(!((fabs(b->pdgId()) == 5 || fabs(b->pdgId()) == 3 || fabs(b->pdgId()) == 1) && fabs(W->pdgId()) == 24)) return false;

      //To identify decay type, check ID of W daughters
      auto Wd1 = W->daughter(event.genparticles,1);
      auto Wd2 = W->daughter(event.genparticles,2);
      //      cout<<"  Wd1 ID = "<<fabs(Wd1->pdgId())<<" Wd2 ID = "<<fabs(Wd2->pdgId())<<endl;

      //try to match the b quarks
      bool matched_b_ak4 = false;

      // Consider AK4 jets first
      for(const auto & jet : *event.jets){
	//	cout<<"   Distance between b and AK4 jets:"<<deltaR(*b,jet)<<endl;
        if(deltaR(*b,jet) <= 0.4) matched_b_ak4 = true;
      }

      bool matched_b_ak8 = false;
      // Now consider AK8 jets
      int idx_matched_topjet = -1;
      int idx = 0;
      for(const auto & jet : *event.topjets){
	//	cout<<"    Distance between b and AK8 jets:"<<deltaR(*b,jet)<<endl;
        if(deltaR(*b,jet) <= 0.8){
          matched_b_ak8 = true;
          idx_matched_topjet = idx;
        }
        idx++;
      }
      //      if(!matched_b_ak4 && !matched_b_ak8) cout<<"!!! b jet is not matched to neither AK4 neither to AK8, eta of b gen = "<<b->eta()<<endl;
      if(!matched_b_ak4 && !matched_b_ak8) return false;

      //Check decaymodes of W
    
      //hadronic
      if(fabs(Wd1->pdgId()) < 7 && fabs(Wd2->pdgId()) < 7){
        if(found_had) return false;
        found_had = true;

        //check if both daughters can be matched by jets
        bool matched_d1_ak4 = false, matched_d2_ak4 = false;
        bool matched_d1_ak8 = false, matched_d2_ak8 = false;
        // Consider AK4 jets first
        for(const auto & jet : *event.jets){
          if(deltaR(*Wd1, jet) <= 0.4) matched_d1_ak4 = true;
          if(deltaR(*Wd2, jet) <= 0.4) matched_d2_ak4 = true;
        }

        // Now consider the one AK8 jet also used for the b-jet
        if(!matched_b_ak8){
          matched_d1_ak8 = false;
          matched_d2_ak8 = false;
        }
        else{
          if(deltaR(*Wd1, event.topjets->at(idx_matched_topjet)) <= 0.8) matched_d1_ak8 = true;
          if(deltaR(*Wd2, event.topjets->at(idx_matched_topjet)) <= 0.8) matched_d2_ak8 = true;
        }

        // if(!(matched_d1 && matched_d2)) return false;
	//        if(!(matched_b_ak4 && matched_d1_ak4 && matched_d2_ak4) && !(matched_b_ak8 && matched_d1_ak8 && matched_d2_ak8))
	  //  cout<<"!!! Hadronic side is not matched, ak4:"<<matched_b_ak4<<" "<<matched_d1_ak4<<" "<<matched_d2_ak4<<" ak8:"<<matched_b_ak8<<" "<<matched_d1_ak8<<" "<<matched_d2_ak8<<endl;
        if(!(matched_b_ak4 && matched_d1_ak4 && matched_d2_ak4) && !(matched_b_ak8 && matched_d1_ak8 && matched_d2_ak8)) return false;
      }

      //leptonic
      else if((abs(Wd1->pdgId()) == 11 || abs(Wd1->pdgId()) == 13) || (abs(Wd2->pdgId()) == 11 || abs(Wd2->pdgId()) == 13)){
        if(found_lep) return false;

        // Escape cases where the W radiates an intermediate photon, that splits into llbar
        if(Wd1->pdgId() == -Wd2->pdgId()){
          // cout << "Entered the escape-part" << endl;
          // Find 2 genparts with 11,12,13,14 that follow each other in the list and don't have the same fabs
          int idx = 0;
          for(const auto & genp : *event.genparticles){
            if(found_lep) break;
            if(abs(genp.pdgId()) >= 11 && abs(genp.pdgId()) <= 14){
              bool is_charged = (abs(genp.pdgId()) == 11 || abs(genp.pdgId()) == 13);
              // cout << "Found a genpart at index " << idx << " with id " << genp.pdgId() << ", is_charged: " << is_charged << endl;

              // if the first one is charged, the second one has to have pdgId of +1 wrt. this genpart
              if(is_charged){
                // cout << "(charged) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) + 1){
                  Wd1 = &genp;
                  Wd2 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
              else{
                // cout << "(neutral) Going to check for next particle in list" << endl;
                if(abs(event.genparticles->at(idx+1).pdgId()) == abs(genp.pdgId()) - 1){
                  Wd2 = &genp;
                  Wd1 = &event.genparticles->at(idx+1);
                  found_lep = true;
                }
              }
            }
            idx++;
          }
	  //          if(!found_lep) cout<<"!!! Not found leptonic side";
          if(!found_lep) return false;
        }

        found_lep = true;

        //Find charged lepton
        auto lep = Wd1;
        auto nu = Wd2;
        if(fabs(Wd2->pdgId()) == 11 || fabs(Wd2->pdgId()) == 13){
          lep = Wd2;
          nu = Wd1;
        }
        if(!(abs(lep->pdgId()) == 11 && abs(nu->pdgId()) == 12) && !(abs(lep->pdgId()) == 13 && abs(nu->pdgId()) == 14)) throw runtime_error("In TTbarSemiLepMatchable: The leptonic W does not decay into a lepton and its neutrino.");

        //check, if lepton can be matched
        bool matched_lep = false;
        if(fabs(lep->pdgId()) == 11){
          for(const auto & ele : *event.electrons){
            if(deltaR(*lep,ele) <= 0.2) matched_lep = true;
          }
        }
        else if(fabs(lep->pdgId()) == 13){
          for(const auto & mu : *event.muons){
            if(deltaR(mu,*lep) <= 0.2) matched_lep = true;
          }
        }
        else throw runtime_error("In TTbarSemiLepMatchable: Lepton from W decay is neither e nor mu.");
        if(!matched_lep) return false;
      }
      //tau-decays
      else return false;
    }
  }

  //  if(!(found_had && found_lep)) cout<<"!!! found_had = "<<found_had<<" found_lep = "<<found_lep<<endl;
  if(!(found_had && found_lep)) return false;
  //  cout<<"### found_had = "<<found_had<<" found_lep = "<<found_lep<<endl;
  return true;
}

////////////////////////////////////////////////////////


