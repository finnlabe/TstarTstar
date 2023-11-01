#include "UHH2/TstarTstar/include/TstarTstarPDFNormHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

TstarTstarPDFNormHists::TstarTstarPDFNormHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  needsOtherMCweightHandling = ctx.get("dataset_version").find("TstarTstar") != std::string::npos;
  //if(needsOtherMCweightHandling) std::cout << "We are signal so we need other idices" << std::endl;
  
  // weight handles
  book<TH1F>("nominal", "nominal", 1, 0, 1);

  // 100 histograms for the PDF stuff
  for(int i=0; i<100; i++){
    std::stringstream ss_name;
    ss_name << "PDF_" << i+1;

    stringstream ss_title;
    ss_title << "S_T [GeV] for PDF No. "  << i+1 << " out of 100" ;

    std::string s_name = ss_name.str();
    std::string s_title = ss_title.str();
    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();

    hist_names[i] = s_name;

    book<TH1F>(char_name, char_title, 1, 0, 1);

  }


}


void TstarTstarPDFNormHists::fill(const Event & event){

  bool debug = false;

  double weight = event.weight;

  hist("nominal")->Fill(0.5,weight);

  float orig_weight = event.genInfo->originalXWGTUP();
  int MY_FIRST_INDEX = 9;
  if(needsOtherMCweightHandling) MY_FIRST_INDEX = 47;
  for(int i=0; i<100; i++){
    double pdf_weight = 0;
    if(event.genInfo->systweights().size() > 0) pdf_weight =event.genInfo->systweights().at(i+MY_FIRST_INDEX);
    const char* name = hist_names[i].c_str();
    hist(name)->Fill(0.5, weight * pdf_weight / orig_weight);
  }

}

TstarTstarPDFNormHists::~TstarTstarPDFNormHists(){}
