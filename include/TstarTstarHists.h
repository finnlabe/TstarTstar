#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ReconstructionHypothesis.h"

namespace uhh2{

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class TstarTstarHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    TstarTstarHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~TstarTstarHists();

  private:
    uhh2::Event::Handle<ReconstructionHypothesis> h_recohyp_;
    
    uhh2::Event::Handle< float > h_M_Tstar_gluon_;
    uhh2::Event::Handle< float > h_M_Tstar_gamma_;
    uhh2::Event::Handle< float > h_DeltaR_toplep_ak8jet1_;
    uhh2::Event::Handle< float > h_DeltaR_tophad_ak8jet1_;
    uhh2::Event::Handle< float > h_DeltaR_toplep_ak8jet2_;
    uhh2::Event::Handle< float > h_DeltaR_tophad_ak8jet2_;
};

}
