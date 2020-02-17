#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/PhotonIds.h"
#include <UHH2/common/include/MuonIds.h>
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarTESTHists.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {

/** bla
 *  
 * blub
 *
 */
class TstarTstarTESTGEN: public AnalysisModule {
public:
    
  explicit TstarTstarTESTGEN(Context & ctx);
  virtual bool process(Event & event) override;

private:
    
  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  std::unique_ptr<Hists> h_GenInfo;
  std::unique_ptr<Hists> h_GenHists;

  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


};


TstarTstarTESTGEN::TstarTstarTESTGEN(Context & ctx){
  
  h_GenInfo.reset(new TstarTstarTESTHists(ctx, "GenInfo"));
  h_GenHists.reset(new TstarTstarGenHists(ctx, "GenHists"));  
    
  ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
 
}


bool TstarTstarTESTGEN::process(Event & event) {
   
  ttgenprod->process(event);

  h_GenInfo->fill(event);
  h_GenHists->fill(event);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarPreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarTESTGEN)

}
