#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/NtupleObjects.h>
#include <UHH2/core/include/LorentzVector.h>
#include <UHH2/common/include/TTbarGen.h>

#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"
#include "UHH2/common/include/ReconstructionHypothesisDiscriminators.h"

class ttbarChi2Discriminator : uhh2::AnalysisModule{

public:
  explicit ttbarChi2Discriminator(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle< std::vector<ReconstructionHypothesis> > h_ttbar_hyps_;
  uhh2::Event::Handle<ReconstructionHypothesis> h_recohyp_;
  uhh2::Event::Handle<bool> h_is_ttbar_reconstructed_;
  float mtoplep_, mtoplep_ttag_;
  float sigmatoplep_, sigmatoplep_ttag_;
  float mtophad_, mtophad_ttag_;
  float sigmatophad_, sigmatophad_ttag_;
  
};



class TstarTstar_Reconstruction : uhh2::AnalysisModule{

public:
  explicit TstarTstar_Reconstruction(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  uhh2::Event::Handle<ReconstructionHypothesis> h_recohyp_;
  uhh2::Event::Handle<float> h_M_Tstar_gluon_;
  uhh2::Event::Handle<float> h_M_Tstar_gamma_;

  uhh2::Event::Handle<float> h_DeltaR_toplep_ak8jet1_;
  uhh2::Event::Handle<float> h_DeltaR_tophad_ak8jet1_;
  uhh2::Event::Handle<float> h_DeltaR_toplep_ak8jet2_;
  uhh2::Event::Handle<float> h_DeltaR_tophad_ak8jet2_;
};
