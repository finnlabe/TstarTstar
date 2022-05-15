#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/MCWeight.h"

enum ParticleID { H=25, W=24, Z=23, t=6, b=5, c=4, s=3, d=2, u=1, g=21};
inline bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}

// Generic Class for Applying SFs
class ScaleFactorsFromHistos : public uhh2::AnalysisModule {

public:
  void LoadHisto(TFile* file, std::string name, std::string hname);
  double Evaluator(std::string hname, double var);

protected:
  std::unordered_map<std::string, std::unique_ptr<TH1F> > histos;

};

// Apply Theory weights //TODO make multiple inheritance
class NLOCorrections : public ScaleFactorsFromHistos {

public:
  explicit NLOCorrections(uhh2::Context& ctx);
  virtual bool process(uhh2::Event&) override;
  double GetPartonObjectPt(uhh2::Event& event, ParticleID objID);

private:
  bool is_Wjets, is_Zjets, is_DY, is_Znn, is2016;

  // Theory Corrections
  const bool do_EWK = true;
  const bool do_QCD_EWK = false;
  const bool do_QCD_NLO  = true;
  const bool do_QCD_NNLO = false;

};
