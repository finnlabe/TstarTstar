// Just a basic macro to create simple plots from given hists.
// author: A.Karavdina (changes by F.Labe)
// date: 24.09.2019
// Run it with following command:
// root -l -b -q FitModelLine.C

void FitModelLine(){
  
  TCanvas *c1 = new TCanvas("chist", "c", 800, 600);
  TLegend *legend = new TLegend(0.2,0.2);

  //Double_t x[] = {700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600};
  //Double_t y[] = {4.92, 1.682, 6.36e-01, 2.62e-1, 1.16e-01, 5.37e-02, 2.61e-02, 1.31e-02, 6.77e-03, 3.59e-03};  
  Double_t x[] = {700, 800, 900, 1100, 1200, 1300, 1400, 1500, 1600};
  Double_t y[] = {3.71e-02, 3.62e-02, 3.47e-02, 3.01e-02, 2.79e-02, 2.58e-02, 2.41e-02, 2.26e-2, 2.14e-02};  

  TGraph *graph = new TGraph(9, x, y);
  TH1 *hist = new TH1D("hist", "hist", 10, 650, 1650);
  for(uint i = 0; i < 10; i++){
    hist->Fill(x[i], y[i]);
  }

  TF1 *f1 = new TF1("f1", "pol 1", 700, 1600);
  TF1 *f2 = new TF1("f2", "pol 2", 700, 1600);
  TF1 *f3 = new TF1("f3", "pol 3", 700, 1600);
  graph->Fit("f1", "N");
  graph->Fit("f2", "N");
  graph->Fit("f3", "N");
  //c1->SetLogy();
  graph->SetTitle("Fit to theory line");
  //graph->SetMaximum(100);
  //graph->SetMinimum(1e-4);
  graph->SetMarkerStyle(2);
  graph->Draw("AP");  

  //hist->Draw("same");  

  f1->SetLineColor(2);
  f1->Draw("same");
  f2->SetLineColor(3);
  f2->Draw("same");
  f3->SetLineColor(4);
  f3->Draw("same");

  legend->AddEntry(f1, "pol 1");
  legend->AddEntry(f2, "pol 2");
  legend->AddEntry(f3, "pol 3");
  legend->Draw();

  c1->SaveAs("Fit_model.pdf");

  cout << "Getting Values..." << endl;
  cout << f3->Eval(1700) << endl;
  cout << f3->Eval(1800) << endl;
  cout << f3->Eval(1900) << endl;
  cout << f3->Eval(2000) << endl;

}
