#include <TString.h>
#include <iostream>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2F.h>
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
#include <sstream>
#include <iostream>

using namespace std;

void ExtractSignalNorm(){

  TString year = "UL16";
  TString file_dir = "/nfs/dust/cms/user/flabe/TstarTstar/data/Preselection/"+year+"/hadded/";
  vector<TString> masspoints = {"700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000", "2250", "2500", "2750"};
  TString signal_base = "TstarTstar_M-";


  for(unsigned int i=0; i<masspoints.size(); i++){
    TString filename;
    filename = file_dir + "/uhh2.AnalysisModuleRunner.MC." + signal_base + masspoints.at(i) + ".root";
  
    TFile* file_in = new TFile(filename, "READ");

    float nominal = ((TH1F*)(file_in->Get("PDFNorm/nominal")))->Integral(); 
    std::ofstream file_out;
    file_out.open("files/SignalNorm_" + year + "_" + signal_base + masspoints.at(i) + ".txt");
    //cout << "nominal = " << nominal << endl;

    for(int i = 1; i<101; i++){
        stringstream ss_name;
        ss_name << "PDFNorm/PDF_" << i;
        string s_name = ss_name.str();
        const char* char_name = s_name.c_str();
        float sum_pdf =  ((TH1F*)(file_in->Get(char_name)))->Integral();
        //cout << "sum event weights pdf no " << i << "= " << sum_pdf << endl;
        float norm_pdf = nominal / sum_pdf;
        file_out << "pdf" << i << " " << norm_pdf << std::endl;;
    }

    delete file_in;
  
  }

}