#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
int plotHist() {

  auto canvas1=new TCanvas("canvas1","canvas1");
  TH1F * betaRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend = new TLegend(0.1,0.99,0.4,0.5);
  gStyle->SetOptStat(0);
  legend->AddEntry((TObject*)0,"gluino mass","");
  for(int i=5; i<1281; i=i*2){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("betaRH",betaRH);
    betaRH->Draw("SAME PLC");
    legend->AddEntry(betaRH,i+std::string(" GeV"),"l");
  }
  TFile inFile("hist2048.root");
  inFile.GetObject("betaRH",betaRH);
  betaRH->Draw("SAME PLC");
  legend->AddEntry(betaRH,"2048 GeV","l");
  gPad->SetLogy();
  betaRH->GetXaxis()->SetTitle("beta");
  betaRH->GetYaxis()->SetTitle("counts");
  legend->Draw();
  canvas1->Print("betaRH.pdf");

  auto canvas2=new TCanvas("canvas2","canvas2");
  TH1F * phiRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend2 = new TLegend(0.8,0.99,0.99,0.5);
  gStyle->SetOptStat(0);
  legend2->AddEntry((TObject*)0,"gluino mass","");
  for(int i=5; i<1281; i=i*2){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("phiRH",phiRH);
    phiRH->Draw("SAME PLC");
    legend2->AddEntry(phiRH,i+std::string(" GeV"),"l");
  }
  TFile inFile2("hist2048.root");
  inFile2.GetObject("phiRH",phiRH);
  phiRH->Draw("SAME PLC");
  legend2->AddEntry(phiRH,"2048 GeV","l");
  phiRH->GetXaxis()->SetTitle("phi [rad]");
  phiRH->GetYaxis()->SetTitle("counts");  
  legend2->Draw();
  canvas2->Print("phiRH.pdf");

  auto canvas3=new TCanvas("canvas3","canvas3");
  TH1F * pTRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend3 = new TLegend(0.8,0.99,0.99,0.5);
  gStyle->SetOptStat(0);
  legend3->AddEntry((TObject*)0,"gluino mass","");
  for(int i=5; i<1281; i=i*2){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("pTRH",pTRH);
    pTRH->Draw("SAME PLC");
    legend3->AddEntry(pTRH,i+std::string(" GeV"),"l");
  }
  TFile inFile3("hist2048.root");
  inFile3.GetObject("pTRH",pTRH);
  pTRH->Draw("SAME PLC");
  legend3->AddEntry(pTRH,"2048 GeV","l");
  gPad->SetLogy();
  gPad->SetLogx();
  pTRH->GetXaxis()->SetTitle("Transverse momentum [GeV]");
  pTRH->GetYaxis()->SetTitle("counts");
  legend3->Draw();
  canvas3->Print("pTRH.pdf");

  auto canvas4=new TCanvas("canvas4","canvas4");
  TH1F * etaRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend4 = new TLegend(0.8,0.99,0.99,0.5);
  gStyle->SetOptStat(0);
  legend4->AddEntry((TObject*)0,"gluino mass","");
   for(int i=5; i<1281; i=i*2){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("etaRH",etaRH);
    etaRH->Draw("SAME PLC");
    etaRH->SetMaximum(22);
    legend4->AddEntry(etaRH,i+std::string(" GeV"),"l");
  }
  TFile inFile4("hist2048.root");
  inFile4.GetObject("etaRH",etaRH);
  etaRH->Draw("SAME PLC");
  legend4->AddEntry(etaRH,"2048 GeV","l");
  etaRH->GetXaxis()->SetTitle("eta");
  etaRH->GetYaxis()->SetTitle("counts");
  legend4->Draw();
  canvas4->Print("etaRH.pdf");
  

  auto canvas5=new TCanvas("canvas5","canvas5");
  TH1F * acceptedPTRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend5 = new TLegend(0.8,0.9,0.99,0.5);
  gStyle->SetOptStat(0);
  legend5->AddEntry((TObject*)0,"gluino mass","");
  for(int i=5; i<1281; i=i*2){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("acceptedPTRH",acceptedPTRH);
    acceptedPTRH->Draw("SAME PLC");
    acceptedPTRH->GetXaxis()->SetTitle("Transverse momentum [GeV]");
    acceptedPTRH->GetYaxis()->SetTitle("counts");
    acceptedPTRH->SetMaximum(45);
    legend5->AddEntry(acceptedPTRH,i+std::string(" GeV"),"l");
  }
  TFile inFile5("hist2048.root");
  inFile5.GetObject("acceptedPTRH",acceptedPTRH);
  acceptedPTRH->Draw("SAME PLC");
  legend5->AddEntry(acceptedPTRH,"2048 GeV","l"); 
  gPad->SetLogx();
  legend5->Draw();
  canvas5->Print("acceptedPTRH.pdf");

  auto canvas6=new TCanvas("canvas6","canvas6");
  TH1F * acceptedBetaRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend6 = new TLegend(0.8,0.9,0.99,0.5);
  gStyle->SetOptStat(0);
  legend6->AddEntry((TObject*)0,"gluino mass","");
  for(int i=5; i<1281; i=i*2){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("acceptedBetaRH",acceptedBetaRH);
    acceptedBetaRH->Draw("SAME PLC");
    acceptedBetaRH->GetXaxis()->SetTitle("beta");
    acceptedBetaRH->GetYaxis()->SetTitle("counts");
    acceptedBetaRH->SetMaximum(20);
    legend6->AddEntry(acceptedBetaRH,i+std::string(" GeV"),"l");
  }
  TFile inFile6("hist2048.root");
  inFile6.GetObject("acceptedBetaRH",acceptedBetaRH);
  acceptedBetaRH->Draw("SAME PLC");
  legend6->AddEntry(acceptedBetaRH,"2048 GeV","l");
  legend6->Draw();
  canvas6->Print("acceptedBetaRH.pdf");

  auto canvas7=new TCanvas("canvas7","canvas7");
  TH1F * efficiencyPTRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend7 = new TLegend(0.8,0.9,0.99,0.5);
  gStyle->SetOptStat(0);
  legend7->AddEntry((TObject*)0,"gluino mass","");
  for(int i=5; i<641; i=i*4){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("efficiencyPTRH",efficiencyPTRH);
    efficiencyPTRH->Draw("SAME PLC");
    efficiencyPTRH->Rebin();    
    efficiencyPTRH->GetXaxis()->SetTitle("Transverse momentum [GeV]");
    efficiencyPTRH->GetYaxis()->SetTitle("Fraction of R-hadrons that hit detector");
    efficiencyPTRH->SetMaximum(1);
    legend7->AddEntry(efficiencyPTRH,i+std::string(" GeV"),"l");
  }
  TFile inFile7("hist2048.root");
  inFile7.GetObject("efficiencyPTRH",efficiencyPTRH);
  efficiencyPTRH->Draw("SAME PLC");
  legend7->AddEntry(efficiencyPTRH,"2048 GeV","l"); 
  efficiencyPTRH->Rebin();    
  gPad->SetLogx();
  gPad->SetLogy();
  legend7->Draw();
  canvas7->Print("efficiencyPTRH.pdf");

  auto canvas8=new TCanvas("canvas8","canvas8");
  TH1F * efficiencyBetaRH;
  TH1::AddDirectory(kFALSE);
  gStyle->SetPalette(kRainBow);
  auto* legend8 = new TLegend(0.8,0.9,0.99,0.5);
  gStyle->SetOptStat(0);
  legend8->AddEntry((TObject*)0,"gluino mass","");
  for(int i=5; i<641; i=i*4){
    TFile inFile(std::string("hist")+i+".root");
    inFile.GetObject("efficiencyBetaRH",efficiencyBetaRH);
    efficiencyBetaRH->Draw("SAME PLC");
    efficiencyBetaRH->Rebin();    
    efficiencyBetaRH->GetXaxis()->SetTitle("beta");
    efficiencyBetaRH->GetYaxis()->SetTitle("Fraction of R-hadrons that hit detector");
    efficiencyBetaRH->SetMaximum(1);
    legend8->AddEntry(efficiencyBetaRH,i+std::string(" GeV"),"l");
  }
  TFile inFile8("hist2048.root");
  inFile8.GetObject("efficiencyBetaRH",efficiencyBetaRH);
  efficiencyBetaRH->Draw("SAME PLC");
  efficiencyBetaRH->Rebin();    
  legend8->AddEntry(efficiencyBetaRH,"2048 GeV","l");
  gPad->SetLogy();
  legend8->Draw();
  canvas8->Print("efficiencyBetaRH.pdf");
  
  return 0;
}
