#include "UHH2/TstarTstar/include/TstarTstarTESTHists.h"
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


TstarTstarTESTHists::TstarTstarTESTHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here


  // Generated Masses 
  book<TH1F>("x_g1", "x_g1", 20, 0, 1);
  book<TH1F>("x_g2", "x_g2", 20, 0, 1);

  book<TH1F>("Tstar_spin", "Tstar_Spin", 4, 0, 2);

  book<TH1F>("N_g_GEN", "N_g_GEN", 10, 0, 10);
  book<TH1F>("PosSigGlu", "PosSigGlu", 10, 0, 10);
  
  is_mc = ctx.get("dataset_type") == "MC";
}


void TstarTstarTESTHists::fill(const Event & event){
  if(!is_mc) return;

  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  assert(event.genparticles);

  GenInfo* genInfo = event.genInfo;
  
  hist("x_g1")->Fill(genInfo->pdf_x1(), weight);
  hist("x_g2")->Fill(genInfo->pdf_x2(), weight);
  
  std::vector<GenParticle> gluons;
  int ngluons = 0;
  for(const GenParticle & gp : *event.genparticles){
    if(gp.pdgId() == 9000005 && (gp.status()==23 || gp.status()==22)){
      hist("Tstar_spin")->Fill(gp.spin(), weight);
    }
    else if(gp.pdgId() == -9000005 && gp.status()>0){
      hist("Tstar_spin")->Fill(gp.spin(), weight);
    }
    else if(gp.pdgId() == 21 && gp.status()>0){
      gluons.push_back(gp);
      ngluons++;
    }
  }

  //cout << "Number of Gluons: " << ngluons << endl;
  hist("N_g_GEN")->Fill(ngluons, weight);

  // sort vector by pt.
  std::vector<GenParticle> gluons_sorted;
  while(gluons.size()>0){
    double maxPt = -1;
    int maxPtPos = -1;
    int index = 0;
    for(const auto & gluon : gluons){
      if(gluon.pt() > maxPt){
	maxPt = gluon.pt();
	maxPtPos = index;
      }
      index++;
    }
    //cout << "Written glu at " << maxPtPos << " in new vec. Old Vec Size is: " << gluons.size();
    gluons_sorted.push_back(gluons.at(maxPtPos));
    gluons.erase(gluons.begin() + maxPtPos);
    //cout << " ... done. New size: " << gluons.size() << endl;
  }

  // cout << "At least got here..." << endl;

  int count_Tstargluons = 0;
  int index = 0;
  for(const auto & gluon : gluons_sorted){
    if(gluon.mother1() < event.genparticles->size() && gluon.mother1() >= 0){
      //cout << "Mother: " << gluon.mother1() << " | Vector size: " << event.genparticles->size() << endl;;
      if((event.genparticles->at(gluon.mother1()).pdgId() == 9000005) || (event.genparticles->at(gluon.mother1()).pdgId() == -9000005)){
	count_Tstargluons++;
	hist("PosSigGlu")->Fill(index, weight);
      }
    }
    index++;
  }
  
  if(count_Tstargluons > 2){cout << "Error, more than 2 gluons from Tstar found!" << endl;}

}

TstarTstarTESTHists::~TstarTstarTESTHists(){}


