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
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/TstarTstar/include/TstarTstarSelections.h"
#include "UHH2/TstarTstar/include/TstarTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarRecoTstarHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenHists.h"
#include "UHH2/TstarTstar/include/TstarTstarGenRecoMatchedHists.h"
#include "UHH2/TstarTstar/include/TstarTstarReconstructionModules.h"
#include "UHH2/TstarTstar/include/ReconstructionTstarHypothesis.h"
#include "UHH2/TstarTstar/include/TstarTstarGenMatch.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#pragma GCC diagnostic push
// turn off the specific warning. Can also use "-Wall"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#pragma GCC diagnostic pop


using namespace std;
using namespace uhh2;

namespace uhh2 {

  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }


/** \brief Module for the T*T*->ttbar gg MC based study
 *
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class TstarTstarMCStudyModule: public AnalysisModule {
public:

  explicit TstarTstarMCStudyModule(Context & ctx);
  virtual bool process(Event & event) override;

private:

  unique_ptr<TTbarSemiLepMatchableSelection> TTbarSemiLepMatchable_selection;
  unique_ptr<Selection> toptagevt_sel;

  // ##### Histograms #####
  // Store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.

  std::unique_ptr<Hists> h_beforeReco, h_beforeReco_ttag, h_beforeReco_nottag, h_afterPrimlep, h_afterHypCreation, h_afterReco_Full, h_afterReco_ttag, h_afterReco_nottag, h_afterGEN, h_afterGEN_onlyttbar, h_afterGEN_onlyttbar_ttag, h_afterGEN_onlyttbar_nottag;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_Full, h_RecoPlots_ttag, h_RecoPlots_nottag;
  std::unique_ptr<TstarTstarRecoTstarHists> h_RecoPlots_GEN, h_RecoPlots_GEN_onlyttbar, h_RecoPlots_GEN_onlyttbar_ttag, h_RecoPlots_GEN_onlyttbar_nottag;

  std::unique_ptr<Hists> h_GEN_Hists, h_GEN_Hists_pre;

  // Bools for Debugging/Options
  bool debug = false;

  bool includeDNNmodel = true;
  // save pt, eta, phi for ttagged jet, leptop, MET, leptopjet, gluon candidates
  bool forDNN = true;
  uhh2::Event::Handle<double> h_DNN_ttaggedjet_pt;
  uhh2::Event::Handle<double> h_DNN_ttaggedjet_eta;
  uhh2::Event::Handle<double> h_DNN_ttaggedjet_phi;
  uhh2::Event::Handle<double> h_DNN_lepton_pt;
  uhh2::Event::Handle<double> h_DNN_lepton_eta;
  uhh2::Event::Handle<double> h_DNN_lepton_phi;
  uhh2::Event::Handle<double> h_DNN_MET_pt;
  uhh2::Event::Handle<double> h_DNN_MET_eta;
  uhh2::Event::Handle<double> h_DNN_MET_phi;
  uhh2::Event::Handle<double> h_DNN_leptopjet_pt;
  uhh2::Event::Handle<double> h_DNN_leptopjet_eta;
  uhh2::Event::Handle<double> h_DNN_leptopjet_phi;
  uhh2::Event::Handle<double> h_DNN_gluon1_pt;
  uhh2::Event::Handle<double> h_DNN_gluon1_eta;
  uhh2::Event::Handle<double> h_DNN_gluon1_phi;
  uhh2::Event::Handle<double> h_DNN_gluon2_pt;
  uhh2::Event::Handle<double> h_DNN_gluon2_eta;
  uhh2::Event::Handle<double> h_DNN_gluon2_phi;

  // DNN spplication stuff
  std::vector<tensorflow::Tensor> lp_tensors_;
  tensorflow::Session* session_3;

  uhh2::Event::Handle<double> h_masspoint;
  uhh2::Event::Handle<double> h_DNN_output;

  // primlep
  uhh2::Event::Handle<FlavorParticle> h_primlep;

  // Modules
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<TstarTstarGenDiscriminator> genmatcher;
  std::unique_ptr<TstarTstarGenDiscriminator> genmatcher_onlyttbar;
  std::unique_ptr<TstarTstar_tgtg_TopTag_Reconstruction> TstarTstarHypCreator;
  std::unique_ptr<TstarTstar_Discrimination> TstarTstarHypSelector;

  std::unique_ptr<uhh2::AnalysisModule> MCWeight;

  // Handles
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  uhh2::Event::Handle<int> h_flag_toptagevent;
  uhh2::Event::Handle<int> h_flag_muonevent;
  uhh2::Event::Handle<std::vector<ReconstructionTstarHypothesis>> h_tstartstar_hyp_vector;
  uhh2::Event::Handle<ReconstructionTstarHypothesis> h_tstartstar_hyp;

  uhh2::Event::Handle<int> jets_thrown_away;


  // bools for channel and stuff. will be read in later
  bool isTrigger;
  bool is_MC;
};




TstarTstarMCStudyModule::TstarTstarMCStudyModule(Context & ctx){

  if(debug) {
    cout << "Hello World from TstarTstarMCStudyModule!" << endl;

    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
   }

  // Dont save everything!
  //ctx.undeclare_all_event_output();

  // 0. Reading in whether MC and if so, which channel
  is_MC = ctx.get("dataset_type") == "MC";


  if(is_HOTVR && debug) cout << "We are using HOTVR jets because they are awesome!!!" << endl;

  // Prepare GEN
  if(is_MC){
    ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    genmatcher.reset(new TstarTstarGenDiscriminator(ctx));
    genmatcher_onlyttbar.reset(new TstarTstarGenDiscriminator(ctx, true));
  }

  //TopTag
  TopJetId topjetID;
  if(is_HOTVR) topjetID = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.56));
  else topjetID = AndId<TopJet>(TopTagMassWindow(), Tau32(0.56));
  toptagevt_sel.reset(new TopTagEventSelection(topjetID));
  h_flag_toptagevent = ctx.declare_event_output<int>("flag_toptagevent");

  // Check which lepton is present and save in
  h_flag_muonevent = ctx.declare_event_output<int>("flag_muonevent");


  // 3. Set up Hists classes:

  h_beforeReco.reset(new TstarTstarHists(ctx, "beforeReco"));
  h_beforeReco_ttag.reset(new TstarTstarHists(ctx, "beforeReco_ttag"));
  h_beforeReco_nottag.reset(new TstarTstarHists(ctx, "beforeReco_nottag"));
  h_afterPrimlep.reset(new TstarTstarHists(ctx, "afterPrimlep"));
  h_afterHypCreation.reset(new TstarTstarHists(ctx, "afterHypCreation"));
  h_afterReco_Full.reset(new TstarTstarHists(ctx, "AfterReco_Full"));
  h_afterReco_ttag.reset(new TstarTstarHists(ctx, "AfterReco_ttag"));
  h_afterReco_nottag.reset(new TstarTstarHists(ctx, "AfterReco_nottag"));
  h_afterGEN.reset(new TstarTstarHists(ctx, "AfterGEN"));
  h_afterGEN_onlyttbar.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar"));
  h_afterGEN_onlyttbar_ttag.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar_ttag"));
  h_afterGEN_onlyttbar_nottag.reset(new TstarTstarHists(ctx, "AfterGEN_onlyttbar_nottag"));

  h_RecoPlots_Full.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_Full"));
  h_RecoPlots_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_ttag"));
  h_RecoPlots_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_nottag"));
  h_RecoPlots_GEN.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN"));
  h_RecoPlots_GEN_onlyttbar.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar"));
  h_RecoPlots_GEN_onlyttbar_ttag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar_ttag"));
  h_RecoPlots_GEN_onlyttbar_nottag.reset(new TstarTstarRecoTstarHists(ctx, "RecoPlots_GEN_onlyttbar_nottag"));

  h_GEN_Hists_pre.reset(new TstarTstarGenHists(ctx, "GEN_Hists_beforeReco"));
  h_GEN_Hists.reset(new TstarTstarGenHists(ctx, "GEN_Hists_AfterReco"));


  //4. Set up ttbar reconstruction
  reco_primlep.reset(new PrimaryLepton(ctx));

  h_tstartstar_hyp_vector = ctx.get_handle<std::vector<ReconstructionTstarHypothesis>>("TstarTstar_Hyp_Vector");
  h_tstartstar_hyp = ctx.get_handle<ReconstructionTstarHypothesis>("TstarTstar_Hyp");

  TstarTstarHypCreator.reset(new TstarTstar_tgtg_TopTag_Reconstruction(ctx, NeutrinoReconstruction, topjetID));
  TstarTstarHypSelector.reset(new TstarTstar_Discrimination(ctx));

  MCWeight.reset(new MCLumiWeight(ctx));

  // 5. Handles for DNN
  if(forDNN){
    h_DNN_ttaggedjet_pt = ctx.declare_event_output<double>("DNN_ttaggedjet_pt");
    h_DNN_ttaggedjet_eta = ctx.declare_event_output<double>("DNN_ttaggedjet_eta");
    h_DNN_ttaggedjet_phi = ctx.declare_event_output<double>("DNN_ttaggedjet_phi");
    h_DNN_lepton_pt = ctx.declare_event_output<double>("DNN_lepton_pt");
    h_DNN_lepton_eta = ctx.declare_event_output<double>("DNN_lepton_eta");
    h_DNN_lepton_phi = ctx.declare_event_output<double>("DNN_lepton_phi");
    h_DNN_MET_pt = ctx.declare_event_output<double>("DNN_MET_pt");
    h_DNN_MET_eta = ctx.declare_event_output<double>("DNN_MET_eta");
    h_DNN_MET_phi = ctx.declare_event_output<double>("DNN_MET_phi");
    h_DNN_leptopjet_pt = ctx.declare_event_output<double>("DNN_leptopjet_pt");
    h_DNN_leptopjet_eta = ctx.declare_event_output<double>("DNN_leptopjet_eta");
    h_DNN_leptopjet_phi = ctx.declare_event_output<double>("DNN_leptopjet_phi");
    h_DNN_gluon1_pt = ctx.declare_event_output<double>("DNN_gluon1_pt");
    h_DNN_gluon1_eta = ctx.declare_event_output<double>("DNN_gluon1_eta");
    h_DNN_gluon1_phi = ctx.declare_event_output<double>("DNN_gluon1_phi");
    h_DNN_gluon2_pt = ctx.declare_event_output<double>("DNN_gluon2_pt");
    h_DNN_gluon2_eta = ctx.declare_event_output<double>("DNN_gluon2_eta");
    h_DNN_gluon2_phi = ctx.declare_event_output<double>("DNN_gluon2_phi");

    h_masspoint = ctx.declare_event_output<double>("h_masspoint");
    h_DNN_output = ctx.declare_event_output<double>("h_DNN_output");

    reco_primlep.reset(new PrimaryLepton(ctx));
    h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  }

  if(includeDNNmodel){
    tensorflow::setLogging("3");
      tensorflow::GraphDef* graphDefFinn = tensorflow::loadGraphDef("/nfs/dust/cms/user/gunnep/HiggsAnalysisUHH_102/CMSSW_10_2_11/src/UHH2/DNNUsageTest/model_Finn.pb");
      // create a new, empty session
      session_3 = tensorflow::createSession(graphDefFinn);

      std::srand(std::time(nullptr)); // Initialize random
  }

  //jets_thrown_away = ctx.declare_event_output<double>("jets_thrown_away");

}


bool TstarTstarMCStudyModule::process(Event & event) {

  if(debug){cout << endl << "TstarTstarMCStudyModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;}

  MCWeight->process(event);
  if(is_MC) ttgenprod->process(event);

  h_beforeReco->fill(event);
  h_GEN_Hists_pre->fill(event);

  // check wether ttag is present
  const bool pass_ttag = toptagevt_sel->passes(event);
  event.set(h_flag_toptagevent, int(pass_ttag));
  if(debug) cout << "Done checking ttag: " << pass_ttag << endl;

  if(pass_ttag){h_beforeReco_ttag->fill(event);}
  else {h_beforeReco_nottag->fill(event);}

  // check lepton channel
  const bool muon_evt = (event.muons->size() == 1);
  event.set(h_flag_muonevent, int(muon_evt));

  // Fix empty handle problem
  // (Dummy values are filled into handles so that they are not empty)
  event.set(h_tstartstar_hyp, ReconstructionTstarHypothesis());
  if(debug){cout << "Finished initialization of Handle Variables" << endl;}

  // ########### TstarTstarReco! ############
  if(debug) cout<<"Starting TstarTstar Reconstruction part"<<endl;

  reco_primlep->process(event);//set "primary lepton"

  bool TstarHypsCreated = false;
  bool bestHypFound = false;

  h_afterPrimlep->fill(event);

  if(debug){cout << "Starting to construct all TstarTstar Hypothesiseseses" << endl;}
  TstarHypsCreated = TstarTstarHypCreator->process(event);
  if(TstarHypsCreated){
    h_afterHypCreation->fill(event);
    bestHypFound = TstarTstarHypSelector->process(event);
    if(bestHypFound){
      h_RecoPlots_Full->fill(event);
      h_afterReco_Full->fill(event);
      h_GEN_Hists->fill(event);

      if(pass_ttag){
	       h_RecoPlots_ttag->fill(event);
	       h_afterReco_ttag->fill(event);
      }
      else{
	       h_RecoPlots_nottag->fill(event);
	       h_afterReco_nottag->fill(event);
      }
    }
  }
  if(!TstarHypsCreated || !bestHypFound){
    // TODO FILL A HIST FOR RECO
  }

  // Filling output for DNN
  // TODO put this in extra file that you pass an hypothesis
  event.set(h_masspoint, 0); // failsafe
  vector<double> DNNInputs;
  if(forDNN && bestHypFound){
    if(debug) cout << "Start filling of DNN stuff" << endl;
    ReconstructionTstarHypothesis best_hyp = event.get(h_tstartstar_hyp);
    ReconstructionHypothesis best_hyp_ttbar = best_hyp.ttbar_hyp();

    // MET
    event.set(h_DNN_MET_eta, best_hyp_ttbar.neutrino_v4().eta());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.neutrino_v4().eta());
    event.set(h_DNN_MET_phi, best_hyp_ttbar.neutrino_v4().phi());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.neutrino_v4().phi());
    event.set(h_DNN_MET_pt, best_hyp_ttbar.neutrino_v4().pt());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.neutrino_v4().pt());
    if(debug) cout << "Done with neutrino." << endl;

    // gluon 1
    event.set(h_DNN_gluon1_eta, best_hyp.gluon1_v4().eta());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp.gluon1_v4().eta());
    event.set(h_DNN_gluon1_phi, best_hyp.gluon1_v4().phi());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp.gluon1_v4().phi());
    event.set(h_DNN_gluon1_pt, best_hyp.gluon1_v4().pt());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp.gluon1_v4().pt());
    if(debug) cout << "Done with gluon1." << endl;

    // gluon 2
    event.set(h_DNN_gluon2_eta, best_hyp.gluon2_v4().eta());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp.gluon2_v4().eta());
    event.set(h_DNN_gluon2_phi, best_hyp.gluon2_v4().phi());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp.gluon2_v4().phi());
    event.set(h_DNN_gluon2_pt, best_hyp.gluon2_v4().pt());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp.gluon2_v4().pt());
    if(debug) cout << "Done with gluon2." << endl;

    // lepton
    event.set(h_DNN_lepton_eta, best_hyp_ttbar.lepton().eta());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.lepton().eta());
    event.set(h_DNN_lepton_phi, best_hyp_ttbar.lepton().phi());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.lepton().phi());
    event.set(h_DNN_lepton_pt, best_hyp_ttbar.lepton().pt());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.lepton().pt());
    if(debug) cout << "Done with lepton." << endl;

      // leptopjet
    event.set(h_DNN_leptopjet_eta, best_hyp_ttbar.blep_v4().eta());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.blep_v4().eta());
    event.set(h_DNN_leptopjet_phi, best_hyp_ttbar.blep_v4().phi());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.blep_v4().phi());
    event.set(h_DNN_leptopjet_pt, best_hyp_ttbar.blep_v4().pt());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.blep_v4().pt());
    if(debug) cout << "Done with leptopjet." << endl;

    // hadtopjet
    event.set(h_DNN_ttaggedjet_eta, best_hyp_ttbar.tophad_v4().eta());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.tophad_v4().eta());
    event.set(h_DNN_ttaggedjet_phi, best_hyp_ttbar.tophad_v4().phi());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.tophad_v4().phi());
    event.set(h_DNN_ttaggedjet_pt, best_hyp_ttbar.tophad_v4().pt());
    if(includeDNNmodel) DNNInputs.push_back(best_hyp_ttbar.tophad_v4().pt());
    if(debug) cout << "Done with ttagjet." << endl;

    // Additional Value for parametrized learning has to be set.
    double masspoint = -1;
    if(is_MC){
      assert(event.genparticles);
      for(const GenParticle & gp : *event.genparticles){
        if((gp.pdgId() == 9000005 || gp.pdgId() == -9000005) && (gp.status()==23 || gp.status()==22)){
          masspoint = inv_mass(gp.v4());
        }
      }
    }
    if(masspoint == -1){ // if no TstarFound (eg Data or Background) set random between 700 and 1600
      masspoint = std::rand() % 1801 + 200; // random between 200 and 2000, inclusively!
    }
    if(includeDNNmodel) DNNInputs.push_back(masspoint);
    event.set(h_masspoint, masspoint);
  }
  // Normalisation of input vector
  vector<double> DNNInputs_mean = {0.002327, 0.008487, 188.7, 0.006092, -0.002269, 488.3, 0.004281, 0.00198, 411.7, 0.0006802, -0.01058, 153.5, -0.0007065, 0.002809, 96.87, 0.004934, -0.007805, 434.3, 1166};
  vector<double> DNNInputs_std = {0.8738, 1.831, 154.1, 1.057, 1.808, 307, 0.9522, 1.813, 270.2, 0.9399, 1.82, 133.4, 0.93, 1.816, 99.26, 0.9954, 1.817, 268.4, 277.8};
  for(uint i = 0; i < DNNInputs.size(); i++){
    DNNInputs.at(i) = (DNNInputs.at(i)-DNNInputs_mean.at(i))/DNNInputs_std.at(i);
  }

  if(is_MC && TstarHypsCreated){
    if(debug){ cout << "Doing GEN matching check" << endl;}
    // ##### GEN Matching
    {
      ReconstructionTstarHypothesis hyp_tmp = event.get(h_tstartstar_hyp);
      if(genmatcher->process(event)){
      	if(debug)cout << "GEN MATCHED!" << endl;
      	h_RecoPlots_GEN->fill(event);
      	h_afterGEN->fill(event);
      	event.set(h_tstartstar_hyp, hyp_tmp);
      }
    }
    {
      ReconstructionTstarHypothesis hyp_tmp = event.get(h_tstartstar_hyp);
      if(genmatcher_onlyttbar->process(event)){
	if(debug)cout << "GEN MATCHED!" << endl;
	h_RecoPlots_GEN_onlyttbar->fill(event);
	h_afterGEN_onlyttbar->fill(event);
	if(pass_ttag){
	  h_RecoPlots_GEN_onlyttbar_ttag->fill(event);
	  h_afterGEN_onlyttbar_ttag->fill(event);
	}
	if(pass_ttag){
	  h_RecoPlots_GEN_onlyttbar_nottag->fill(event);
	  h_afterGEN_onlyttbar_nottag->fill(event);
	}
	event.set(h_tstartstar_hyp, hyp_tmp);
      }
    }
  }


  if(forDNN && bestHypFound && includeDNNmodel){
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::Tensor input(tensorflow::DT_FLOAT, {1, 19});
    for (int i = 0; i < 19; i++){
      input.matrix<float>()(0, i) = DNNInputs.at(i);
    }

    //DENSE FINN
    tensorflow::run(session_3, {{"dense_1_input",input}}, {"dense_4/Tanh"}, &outputs);

    for (int i = 0; i<1 ; i++){
      event.set(h_DNN_output, outputs[0].matrix<float>()(0, i)); // Output to console.
    }
  }

  if(debug){cout << "Done ##################################" << endl;}
  return TstarHypsCreated;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the TstarTstarMCStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TstarTstarMCStudyModule)

}
