#include "UHH2/TstarTstar/include/TstarTstarJetCorrectionHists.h"
#include "UHH2/core/include/Event.h"
#include <UHH2/core/include/Utils.h>
#include <UHH2/common/include/Utils.h>
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

TstarTstarJetCorrectionHists::TstarTstarJetCorrectionHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

    book<TH1D>("ptrec_o_ptgen__AK4_leading", "leading AK4 p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__AK4_second", "second AK4 p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__AK4_all", "all AK4 p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);

    book<TH1D>("ptrec_o_ptgen__HOTVR_leading", "leading HOTVR p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__HOTVR_second", "second HOTVR p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);
    book<TH1D>("ptrec_o_ptgen__HOTVR_all", "all HOTVR p^{REC}_{T} / p^{GEN}_{T}", 100, 0, 2);

    const int n_HOTVRpt_bins = 14;
    double HOTVRpt_bins[n_HOTVRpt_bins] = {0, 200, 225, 250, 275, 300, 350, 400, 500, 600, 800, 1000, 1500, 2000};
    book<TH2D>("2D_jetresponse_genjetpt_HOTVR", ";HOTVR jet response; HOTVR GEN pt", 100, 0, 2, n_HOTVRpt_bins-1, HOTVRpt_bins);

    const int n_AK4pt_bins = 14;
    double AK4pt_bins[n_AK4pt_bins] = {0, 30, 60, 90, 120, 150, 200, 250, 300, 400, 500, 750, 1000, 1500};
    book<TH2D>("2D_jetresponse_genjetpt_AK4", ";HOTVR jet response; HOTVR GEN pt", 100, 0, 2, n_AK4pt_bins-1, AK4pt_bins);

    matching_radius = 0.3;

}


void TstarTstarJetCorrectionHists::fill(const Event & event){

    double weight = event.weight;

    // AK4 part
    int AK4index = 0; // this will be used to determine which is the first and second

    for(const auto & AK4recojet : *event.jets){

        // we have the jet, so we need to find the pt of the matching GEN jet

        double mindR = 99999;
        double pt_of_matched_GENjet = 0;
        for(const auto & AK4genjet : *event.genjets){

            double this_dR = deltaR(AK4recojet, AK4genjet);
            if(this_dR < mindR) {
                mindR = this_dR;
                pt_of_matched_GENjet = AK4genjet.pt();
            }

        }

        if (mindR < matching_radius) {
            if(AK4index == 0) hist("ptrec_o_ptgen__AK4_leading")->Fill( AK4recojet.pt() / pt_of_matched_GENjet, weight);
            if(AK4index == 1) hist("ptrec_o_ptgen__AK4_second")->Fill( AK4recojet.pt() / pt_of_matched_GENjet, weight);
            hist("ptrec_o_ptgen__AK4_all")->Fill( AK4recojet.pt() / pt_of_matched_GENjet, weight);
            ((TH2D*)hist("2D_jetresponse_genjetpt_AK4"))->Fill( AK4recojet.pt() / pt_of_matched_GENjet, pt_of_matched_GENjet, weight);
        }

        AK4index++;
    }

    // HOTVR part
    int HOTVRindex = 0; // this will be used to determine which is the first and second

    for(const auto & HOTVRrecojet : *event.topjets){

        // we have the jet, so we need to find the pt of the matching GEN jet

        double mindR = 99999;
        double pt_of_matched_GENjet = 0;
        for(const auto & HOTVRgenjet : *event.gentopjets){

            double this_dR = deltaR(HOTVRrecojet, HOTVRgenjet);
            if(this_dR < mindR) {
                mindR = this_dR;
                pt_of_matched_GENjet = HOTVRgenjet.pt();
            }

        }

        if (mindR < matching_radius) {
            if(HOTVRindex == 0) hist("ptrec_o_ptgen__HOTVR_leading")->Fill( HOTVRrecojet.pt() / pt_of_matched_GENjet, weight);
            if(HOTVRindex == 1) hist("ptrec_o_ptgen__HOTVR_second")->Fill( HOTVRrecojet.pt() / pt_of_matched_GENjet, weight);
            hist("ptrec_o_ptgen__HOTVR_all")->Fill( HOTVRrecojet.pt() / pt_of_matched_GENjet, weight);
            ((TH2D*)hist("2D_jetresponse_genjetpt_HOTVR"))->Fill( HOTVRrecojet.pt() / pt_of_matched_GENjet, pt_of_matched_GENjet, weight);

        }
        HOTVRindex++;
    }

}


TstarTstarJetCorrectionHists::~TstarTstarJetCorrectionHists(){}
