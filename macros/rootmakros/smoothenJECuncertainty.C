#include <TString.h>
#include <iostream>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TROOT.h>
#include <TKey.h>
#include <TLatex.h>
#include <TClass.h>
#include <fstream>
#include <string>

using namespace std;

void smoothenJECuncertainty(TString year = "UL18", TString channel = "mu", TString region = "ValidationRegion"){

    // method to "smoothen" the JEC uncertainty by fitting a pol0 or pol1 to it

    TString filename_base = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/"+year+"/hadded/uhh2.AnalysisModuleRunner.MC.";
    TString filename_JECup = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/"+year+"/JEC_up/hadded/uhh2.AnalysisModuleRunner.MC.";
    TString filename_JECdown = "/nfs/dust/cms/user/flabe/TstarTstar/data/DNN/"+year+"/JEC_down/hadded/uhh2.AnalysisModuleRunner.MC.";

    TString histfolder = region + "_" + channel;
    TString histname = "pt_ST_nominal";

    // creating sample list
    vector<TString> samples = {"TTbar", "ST"};
    vector<TString> masspoints = {"700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000", "2250", "2500", "2750", "3000"};
    for (auto mass : masspoints) {
        samples.push_back("TstarTstar_M-" + mass);
        samples.push_back("TstarTstar_Spin32_M-" + mass);
    }

    TString fitfuncname = "pol0";

    for (unsigned int i=0; i<samples.size(); i++) {

        std::cout << "Working on sample " << samples.at(i) << std::endl;

        // open the three files (nominal, JEC up and down)
        TFile* f_central = new TFile(filename_base + samples.at(i) + ".root", "READ");
        TFile* f_up = new TFile(filename_JECup + samples.at(i) + ".root", "READ");
        TFile* f_down = new TFile(filename_JECdown + samples.at(i) + ".root", "READ");

        // get the three histograms (nominal, JEC up and down)
        TH1D *h_central = (TH1D*) f_central->Get(histfolder + "/" + histname);
        TH1D *h_up = (TH1D*) f_up->Get(histfolder + "/" + histname);
        TH1D *h_down = (TH1D*) f_down->Get(histfolder + "/" + histname);
        
        // calculate the residual
        h_up->Divide(h_central);
        h_down->Divide(h_central);

        // perform a fit to the residual
        TF1* fit_up = new TF1("fit_up", fitfuncname, 0, 6000);
        h_up->Fit("fit_up", "N", "", 500, 6000);
        TF1* fit_down = new TF1("fit_down", fitfuncname, 0, 6000);
        h_down->Fit("fit_down", "N", "", 500, 6000);

        // take the nominal and vary it up / down by that fit result -> this gives the JEC / JER up / down histogram
        TH1D *h_new_up = (TH1D*)h_central->Clone();
        h_new_up->Multiply(fit_up);
        TH1D *h_new_down = (TH1D*)h_central->Clone();
        h_new_down->Multiply(fit_down);

        // save the resulting histogram
        TFile* f_out = new TFile("/nfs/dust/cms/user/flabe/TstarTstar/ULegacy/CMSSW_10_6_28/src/UHH2/TstarTstar/macros/rootmakros/files/smoothJEC/" + region + "_smoothJEC_" + year + "_" + channel + "_" + samples.at(i) + ".root", "RECREATE");
        h_new_up->SetName(samples.at(i)+"_smoothJEC_up");
        h_new_down->SetName(samples.at(i)+"_smoothJEC_down");
        h_new_up->Write();
        h_new_down->Write();
        delete f_out;

        // and make some plots
        TCanvas *can = new TCanvas("canvas", "c", 600, 600);
        TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
        leg->SetTextSize(0.03);
        leg->SetBorderSize(0);
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

        // plot styling
        h_up->GetYaxis()->SetRangeUser(0.7, 1.3);
        h_up->GetXaxis()->SetRangeUser(600, 6000);
        h_up->GetXaxis()->SetTitle( h_up->GetTitle() );
        h_up->GetYaxis()->SetTitle( "variation / nominal" );

        // plotting the "before" up and down
        h_up->SetLineColor(1);
        h_up->Draw("hist");
        leg->AddEntry(h_up, "JEC variation before", "l");
        h_down->SetLineColor(1);
        h_down->Draw("hist same");
        //leg->AddEntry(h_down, "JEC down", "l");

        // plotting fit functions
        fit_up->Draw("same");
        leg->AddEntry(fit_up, "fit", "l"); // just add one to legend
        fit_down->Draw("same");
        //leg->AddEntry(fit_down, "fit down", "l");

        // and the new functions
        h_new_up->Divide(h_central);
        h_new_down->Divide(h_central);
        h_new_up->SetLineColor(3);
        h_new_up->Draw("hist same");
        leg->AddEntry(h_new_up, "new JEC variation", "l");
        h_new_down->SetLineColor(3);
        h_new_down->Draw("hist same");
        //leg->AddEntry(h_new_down, "new down", "l");

        auto line = TLine(600.1,1,6000,1);
        line.Draw();

        leg->Draw();

        can->SaveAs("plots/JECsmooth_" + samples.at(i) + "_" + region + "_" + year + "_" + channel + ".pdf");

    }

}