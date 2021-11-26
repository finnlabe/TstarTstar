#include <UHH2/TstarTstar/include/ElecTriggerSF.h>

ElecTriggerSF::ElecTriggerSF(uhh2::Context & ctx, std::string var_, TString pe, TString year):var(var_){

  h_ele_weight      = ctx.declare_event_output<float>("weight_sfelec_trigger");
  h_ele_weight_up   = ctx.declare_event_output<float>("weight_sfelec_trigger_up");
  h_ele_weight_down = ctx.declare_event_output<float>("weight_sfelec_trigger_down");

  pteta = pe;
  if(pteta != "pt" && pteta != "eta" && pteta != "eta_ptbins"){
    cout << "Warning: You should select 'pt', 'eta' or 'eta_ptbins' for ElecTriggerSF class" << endl;
    return;
  }

  auto dataset_type = ctx.get("dataset_type");
  isMC = dataset_type == "MC";
  if (!isMC) {
    cout << "Warning: MCElecScaleFactor will not have an effect on "
    <<" this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }

  TString dir = "/nfs/dust/cms/user/paaschal/UHH2_102X_v2/CMSSW_10_2_17/src/UHH2/MTopJet/ScaleFactors/Electrons/";
  TString filename = dir + "ElecTriggerSF"+year+".root";
  TFile *file = new TFile(filename);

  TString histname = "Central";
  TString histUp = "Up";
  TString histDown = "Down";
  if(var == "up")        histname = "Up";
  else if(var == "down") histname = "Down";

  // for the options pt and eta the low and high histograms are exactly the same
  // its just to make it easier to code
  if(pteta == "pt"){
    cout << "ElecTriggerSF: You selected the pt dependent SF: " << histname << endl;
    h_sf_lo = (TH1F*)file->Get(histname+"_pt");
    h_sf_me = (TH1F*)file->Get(histname+"_pt");
    h_sf_hi = (TH1F*)file->Get(histname+"_pt");

    h_sf_lo_up = (TH1F*)file->Get(histUp+"_pt");
    h_sf_me_up = (TH1F*)file->Get(histUp+"_pt");
    h_sf_hi_up = (TH1F*)file->Get(histUp+"_pt");

    h_sf_lo_down = (TH1F*)file->Get(histDown+"_pt");
    h_sf_me_down = (TH1F*)file->Get(histDown+"_pt");
    h_sf_hi_down = (TH1F*)file->Get(histDown+"_pt");
  }
  else if(pteta == "eta"){
    cout << "ElecTriggerSF: You selected the eta dependent SF: " << histname << endl;
    h_sf_lo = (TH1F*)file->Get(histname+"_eta");
    h_sf_me = (TH1F*)file->Get(histname+"_eta");
    h_sf_hi = (TH1F*)file->Get(histname+"_eta");

    h_sf_lo_up = (TH1F*)file->Get(histUp+"_eta");
    h_sf_me_up = (TH1F*)file->Get(histUp+"_eta");
    h_sf_hi_up = (TH1F*)file->Get(histUp+"_eta");

    h_sf_lo_down = (TH1F*)file->Get(histDown+"_eta");
    h_sf_me_down = (TH1F*)file->Get(histDown+"_eta");
    h_sf_hi_down = (TH1F*)file->Get(histDown+"_eta");
  }
  else{
    cout << "ElecTriggerSF: You selected the eta (pt bins) dependent SF: " << histname << endl;
    h_sf_lo = (TH1F*)file->Get(histname+"_eta_lowpt");
    h_sf_me = (TH1F*)file->Get(histname+"_eta_midpt");
    h_sf_hi = (TH1F*)file->Get(histname+"_eta_highpt");

    h_sf_lo_up = (TH1F*)file->Get(histUp+"_eta_lowpt");
    h_sf_me_up = (TH1F*)file->Get(histUp+"_eta_midpt");
    h_sf_hi_up = (TH1F*)file->Get(histUp+"_eta_highpt");

    h_sf_lo_down = (TH1F*)file->Get(histDown+"_eta_lowpt");
    h_sf_me_down = (TH1F*)file->Get(histDown+"_eta_midpt");
    h_sf_hi_down = (TH1F*)file->Get(histDown+"_eta_highpt");
  }
}

bool ElecTriggerSF::process(uhh2::Event & event){

  if(event.isRealData || event.electrons->size() < 1){
    event.set(h_ele_weight, 1.);
    event.set(h_ele_weight_up, 1.);
    event.set(h_ele_weight_down, 1.);
    return true;
  }

  if(!isMC) return true;
  if(pteta != "pt" && pteta != "eta" && pteta != "eta_ptbins") return false;

  if(event.electrons->size() < 1) return false;
  double pt = event.electrons->at(0).pt();
  double eta = event.electrons->at(0).eta();
  double UsedVariable;
  if(pteta == "pt") UsedVariable = pt;
  else              UsedVariable = eta;

  // do not set SF if electron out of range
  // this can happen if recsel is not passed
  if(fabs(eta) > 2.4 || pt < 55){
    event.set(h_ele_weight, 1.);
    event.set(h_ele_weight_up, 1.);
    event.set(h_ele_weight_down, 1.);
    return true;
  }

  int bin = 0;
  double sf = 1.0; double sf_up = 0.0; double sf_down = 0.0;
  if(pt < 120){
    bin = h_sf_lo->GetXaxis()->FindBin(UsedVariable);
    sf = h_sf_lo->GetBinContent(bin);
    sf_up = h_sf_lo_up->GetBinContent(bin);
    sf_down = h_sf_lo_down->GetBinContent(bin);
  }
  else if(pt > 120 && pt < 200){
    bin = h_sf_me->GetXaxis()->FindBin(UsedVariable);
    sf = h_sf_me->GetBinContent(bin);
    sf_up = h_sf_me_up->GetBinContent(bin);
    sf_down = h_sf_me_down->GetBinContent(bin);
  }
  else if(pt > 200){
    bin = h_sf_hi->GetXaxis()->FindBin(UsedVariable);
    sf = h_sf_hi->GetBinContent(bin);
    sf_up = h_sf_hi_up->GetBinContent(bin);
    sf_down = h_sf_hi_down->GetBinContent(bin);
  }

  event.set(h_ele_weight, sf);
  event.set(h_ele_weight_up, sf_up);
  event.set(h_ele_weight_down, sf_down);

  if (var == "up") {
    event.weight *= sf_up;
  } else if (var == "down") {
    event.weight *= sf_down;
  } else {
    event.weight *= sf;
  }

  return true;
}

// =============================================================================
// new method from Alex F (Arne & Anna SF)
// =============================================================================

ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString path_, TString syst_direction_): path(path_), syst_direction(syst_direction_) {

  h_ele_weight      = ctx.declare_event_output<float>("weight_sfelec_trigger");
  h_ele_weight_up   = ctx.declare_event_output<float>("weight_sfelec_trigger_up");
  h_ele_weight_down = ctx.declare_event_output<float>("weight_sfelec_trigger_down");

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  year = extract_year(ctx);
  TString yeartag = "2016";
  if(year == Year::is2017v1 || year == Year::is2017v2) yeartag = "2017";
  else if(year == Year::is2018) yeartag = "2018";
  unique_ptr<TFile> file_pt1, file_pt2;
  if(yeartag == "2016"){
    file_pt1.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt30to175_withsyst.root","READ"));
    file_pt2.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt175toInf.root","READ"));
  }
  else if(yeartag == "2017" || yeartag == "2018"){
    file_pt1.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt30to200_withsyst.root","READ"));
    file_pt2.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt200toInf.root","READ"));
  }
  else throw runtime_error("invalid year");

  if(yeartag == "2016") pt_bins = {30, 175};
  else if(yeartag == "2017" || yeartag == "2018") pt_bins = {30, 200};

  g_sf_pt1.reset((TGraphAsymmErrors*)file_pt1->Get("ScaleFactors"));
  g_sf_pt2.reset((TGraphAsymmErrors*)file_pt2->Get("ScaleFactors"));
}

bool ElectronTriggerWeights::process(Event & event){

  if(event.isRealData || event.electrons->size() < 1){
    event.set(h_ele_weight, 1.);
    event.set(h_ele_weight_up, 1.);
    event.set(h_ele_weight_down, 1.);
    return true;
  }

  const Electron ele = event.electrons->at(0);
  double eta = ele.eta();
  double pt = ele.pt();

  // do not set SF if electron out of range
  // this can happen if recsel is not passed
  if(fabs(eta) > 2.4 || pt < pt_bins[0]){
    event.set(h_ele_weight, 1.);
    event.set(h_ele_weight_up, 1.);
    event.set(h_ele_weight_down, 1.);
    return true;
  }

  // find number of correct eta bin
  int idx = 0;
  if(pt < pt_bins[1]){
    bool keep_going = true;
    while(keep_going){
      double x,y;
      g_sf_pt1->GetPoint(idx,x,y);
      keep_going = eta > x + g_sf_pt1->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else {
    bool keep_going = true;
    while(keep_going){
      double x,y;
      g_sf_pt2->GetPoint(idx,x,y);
      keep_going = eta > x + g_sf_pt2->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  //access scale factors, add 2% t&p systematic uncertainty
  double sf, sf_up, sf_down, dummy_x;
  double stat_up = -1., stat_down = -1., tp = 0.0, total_up = -1., total_down = -1.;
  if(pt < pt_bins[1]){
    g_sf_pt1->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_pt1->GetErrorYhigh(idx);
    stat_down = g_sf_pt1->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }
  else {
    g_sf_pt2->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_pt2->GetErrorYhigh(idx);
    stat_down = g_sf_pt2->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }

  event.set(h_ele_weight, sf);
  event.set(h_ele_weight_up, sf_up);
  event.set(h_ele_weight_down, sf_down);

  if (syst_direction == "up") {
    event.weight *= sf_up;
  } else if (syst_direction == "down") {
    event.weight *= sf_down;
  } else {
    event.weight *= sf;
  }

  return true;
}
