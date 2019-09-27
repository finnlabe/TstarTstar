#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace uhh2 {
  class TwoDCut : public Selection {
  public:
    explicit TwoDCut(float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {}
    virtual bool passes(const Event&) override;
    
  private:
    float min_deltaR_, min_pTrel_;
  };
  ////

  class TTbarSemiLepMatchableSelection: public uhh2::Selection{
  public:
    TTbarSemiLepMatchableSelection();
    ~TTbarSemiLepMatchableSelection(){};
    virtual bool passes(const uhh2::Event & event);
  };

  /* class TTbarAllHadMatchableSelection: public uhh2::Selection{ */
  /* public: */
  /*   TTbarAllHadMatchableSelection(); */
  /*   ~TTbarAllHadMatchableSelection(){}; */
  /*   virtual bool passes(const uhh2::Event & event); */
  /* }; */

}
