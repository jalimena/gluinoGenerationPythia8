#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "tdrstyle.C"
int plotHist(char* size, char* position) {

  setTDRStyle(); 
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.05);

  //plot beta, accepted beta, p, acceptred p and angle distributions
  //acceptedBeta and acceptedP are total number of accepted particles, not fraction that are accepted 
  const char* name[7]={"betaRH","pRH","acceptedBetaRH","acceptedPRH","phiRH","thetaRH","etaRH"};//used for filenames
  const char* hist[7]={"betaRH","pRH","acceptedBetaRH","acceptedPRH","phiRH","thetaRH","etaRH"};//used to declare TH1F
  const char* axis[7]={"#beta","Momentum [GeV]","#beta","Momentum [GeV]","#phi","#theta","#eta"};
  int logScale[7] = {0,2,0,2,0,0,0};//1 for logy scale, 2 for logx scale, 3 for logx and logy, 0 for neither
  int mass[10]={2560,1280,640,320,160,80,40,20,10,5};
  //set colour and style for each mass point
  Color_t fillColour[10] = {kRed,kOrange+1,kOrange,kSpring-9,kGreen,kCyan,kCyan-2,kBlue,kViolet-5,kMagenta};
  Color_t lineColour[10] = {kRed+1,kOrange+2,kOrange-3,kSpring-8,kGreen+1,kCyan+1,kCyan-2,kBlue,kViolet-5,kMagenta};
  Style_t fillStyle[10] = {3150,3151,3152,3153,3154,3156,3157,3158,3159,3159};
  double max[7] = {pow(10,5),55000,2700,300,10400,35000,13000};

  for(int i=0; i<7; i++){//for each property being plotted
    auto canvas=new TCanvas("canvas","canvas");
    TH1F* hist[i];
    TH1::AddDirectory(kFALSE);
    auto* legend = new TLegend(0.65,0.45,0.9,0.9); //beta, theta
    auto* legendM = new TLegend(0.25,0.45,0.5,0.9); //momentum
    legend->SetBorderSize(0);
    legendM->SetBorderSize(0);
    legend->AddEntry((TObject*)0,"gluino mass","");
    legendM->AddEntry((TObject*)0,"gluino mass","");

    for(int j=9; j>-1; j-=1){//for each mass point
      TFile inFile(std::string("hist")+mass[j]+"_"+size+"_"+position+".root");
      inFile.GetObject(name[i],hist[i]);//get saved histogram
      if(j==9){
	hist[i]->Draw("");
	hist[i]->SetMaximum(max[i]);
	hist[i]->SetTitle((std::string(";")+axis[i]+";Entries").c_str());
      }
      else{hist[i]->Draw("SAME");}//plot all masses on the same canvas
      hist[i]->SetLineColor(lineColour[j]);
      hist[i]->SetLineWidth(3);
      hist[i]->GetXaxis()->SetTitleOffset(1.1);
      hist[i]->GetXaxis()->SetTitleSize(0.045);
      legend->AddEntry(hist[i],mass[j]+std::string(" GeV"),"l");
      legendM->AddEntry(hist[i],mass[j]+std::string(" GeV"),"l");
    }
    //add log scale if needed
    if(logScale[i]==1){gPad->SetLogy();}
    else if(logScale[i]==2){gPad->SetLogx();}
    else if(logScale[i]==3){gPad->SetLogy();gPad->SetLogx();}

    if(i==1){ //momentum plot
      legendM->Draw();
    }
    else legend->Draw();

    if(i==5){//draw angular acceptance on theta plot
      double pi=3.14159265358979;
      double thetaAcceptance0low=pi/2+atan(1/7.38);
      double thetaAcceptance0up=pi/2-atan(1/7.38);
      double thetaAcceptance1low=atan(2/10.91);
      double thetaAcceptance1up=atan(4/10.91);
      TH1F * thetaAcceptanceHist0 = new TH1F("thetaAcceptanceHist0","thetaAcceptanceHist0",1,thetaAcceptance0low,thetaAcceptance0up);
      thetaAcceptanceHist0->SetBinContent(1,36000);
      thetaAcceptanceHist0->SetFillColor(kGray);
      thetaAcceptanceHist0->SetFillStyle(3154);
      thetaAcceptanceHist0->Draw("B SAME");
      TH1F * thetaAcceptanceHist1 = new TH1F("thetaAcceptanceHist1","thetaAcceptanceHist1",1,thetaAcceptance1low,thetaAcceptance1up);
      thetaAcceptanceHist1->SetBinContent(1,36000);
      thetaAcceptanceHist1->Draw("B SAME");
      thetaAcceptanceHist1->SetFillColor(kGray);
      thetaAcceptanceHist1->SetFillStyle(3154);

      TPaveText *pt1 = new TPaveText(.18,.75,.35,.9,"brNDC");
      pt1->AddText("Absorber");
      pt1->AddText("position 1");
      pt1->AddText("acceptance");
      pt1->Draw();
      pt1->SetBorderSize(0);
      pt1->SetFillStyle(0);

      TPaveText *pt0 = new TPaveText(.46,.75,.63,.9,"brNDC");
      pt0->AddText("Absorber");
      pt0->AddText("position 0");
      pt0->AddText("acceptance");
      pt0->Draw();
      pt0->SetBorderSize(0);
      pt0->SetFillStyle(0);
    }
    canvas->Print((name[i]+std::string(size)+"_"+position+".pdf").c_str());
  }

  //angular acceptance plots (fraction of particles that are accepted)
  auto canvas1=new TCanvas("canvas1","canvas1");
  TEfficiency * betaAcceptance;
  TH1::AddDirectory(kFALSE);
  auto *legend1 = new TLegend(0.7,0.94,0.85,0.75);
  legend1->AddEntry((TObject*)0,"gluino mass","");
  for(int j=1; j<10; j+=4){
    TFile infile1(std::string("hist")+mass[j]+"_"+size+"_"+position+".root");
    infile1.GetObject("betaRH_clone",betaAcceptance);
    if(j==1){
      betaAcceptance->Draw("AC E3");
      betaAcceptance->SetTitle(";#beta;Angular acceptance");
      //change graph limits
      gPad->Update();
      auto graph=betaAcceptance->GetPaintedGraph();
      graph->SetMinimum(0);
      graph->SetMaximum(0.04);
      TAxis *xaxis = graph->GetXaxis();
      xaxis->SetLimits(0.1,1.);
      gPad->Update();
    }
    else{betaAcceptance->Draw("C SAME E3");}
    legend1->AddEntry(betaAcceptance,mass[j]+std::string(" GeV"),"f");
    betaAcceptance->SetFillStyle(fillStyle[j]);
    betaAcceptance->SetFillColor(fillColour[j]);
    betaAcceptance->SetLineColor(lineColour[j]);
  }
  legend1->Draw();
  legend1->SetBorderSize(0);
  //write absorber position and size
  TPaveText *pt0 = new TPaveText(.4,.93,.65,.85,"brNDC");
  pt0->AddText((std::string("Absorber position ")+position).c_str());
  pt0->AddText((size+std::string("m x ")+size+"m x 2m").c_str());
  pt0->Draw();
  pt0->SetBorderSize(0);
  pt0->SetFillStyle(0);
  canvas1->Print((std::string("betaAcceptance_")+size+"_"+position+".pdf").c_str());


  auto canvas2=new TCanvas("canvas2","canvas2");
  TEfficiency * pAcceptance; 
  TH1::AddDirectory(kFALSE);
  auto* legend2 = new TLegend(0.8,0.9,0.99,0.5);
  legend2->AddEntry((TObject*)0,"gluino mass","");
  for(int j=1; j<10; j+=4){
    TFile infile2(std::string("hist")+mass[j]+"_"+size+"_"+position+".root");
    infile2.GetObject("pRH_clone",pAcceptance);
    if(j==1){
      pAcceptance->Draw("AC E3");
      pAcceptance->SetTitle(";Momentum [GeV];Fraction of R-hadrons that hit detector");
      //change graph limits
      gPad->Update();
      auto graph=betaAcceptance->GetPaintedGraph();
      graph->SetMinimum(0);
      graph->SetMaximum(0.04);
      TAxis *xaxis = graph->GetXaxis();
      xaxis->SetLimits(0.1,1000.);
      gPad->Update();
    }
    pAcceptance->Draw("C SAME E3");
    legend2->AddEntry(pAcceptance,mass[j]+std::string(" GeV"),"f");
    pAcceptance->SetFillStyle(fillStyle[j]);
    pAcceptance->SetFillColor(fillColour[j]);
    pAcceptance->SetLineColor(lineColour[j]);
  }
  legend2->Draw();
  gPad->SetLogy();
  gPad->SetLogx();
  canvas2->Print((std::string("pAcceptance_")+size+"_"+position+".pdf").c_str());

  //plot absorption effienciency convolved with kinematic acceptance
  auto canvas3=new TCanvas("canvas3","canvas3");
  TGraphAsymmErrors * pEfficiency; 
  TH1::AddDirectory(kFALSE);
  auto* legend3 = new TLegend(0.2,0.94,0.35,0.6);
  legend3->AddEntry((TObject*)0,"gluino mass","");
  //for(int j=1; j<10; j+=2){//only plot every other mass point so eaier to see 
  for(int j=9; j>-1; j-=2){//for each mass point
    TFile infile3(std::string("hist")+mass[j]+"_"+size+"_"+position+".root");
    infile3.GetObject("Graph;4",pEfficiency);
    if(j==1){  pEfficiency->Draw("AL3");
      pEfficiency->SetMaximum(0.02);//needs changing for different absorber positions/sizes
      pEfficiency->SetMinimum(pow(10,-7));//needs changing for different absorber positions/sizes
      pEfficiency->SetTitle(";Momentum [GeV];Absorption Efficiency");
      TAxis *xaxis = pEfficiency->GetXaxis();
      xaxis->SetLimits(3.,900.); 
    }
    else{pEfficiency->Draw("L3 SAME");}
    legend3->AddEntry(pEfficiency,mass[j]+std::string(" GeV"),"f");
    pEfficiency->SetFillStyle(fillStyle[j]);
    pEfficiency->SetFillColor(fillColour[j]);
    pEfficiency->SetLineColor(lineColour[j]);
  }
  //write detector size and position
  TPaveText *pt1 = new TPaveText(.4,.93,.65,.85,"brNDC");
  pt1->AddText((std::string("Detector position ")+position).c_str());
  pt1->AddText((size+std::string("m x ")+size+"m x 2m").c_str());
  pt1->Draw();
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  legend3->Draw();
  legend3->SetBorderSize(0);
  gPad->SetLogx();
  canvas3->Print((std::string("pEfficiency_")+size+"_"+position+".pdf").c_str());

  auto canvas4=new TCanvas("canvas4","canvas4");
  TGraphAsymmErrors * betaEfficiency; 
  TH1::AddDirectory(kFALSE);
  auto* legend4 = new TLegend(0.69,0.53,0.94,0.93);//use for position 0
  //auto* legend4 = new TLegend(0.2,0.53,0.45,0.93);//use for position 1

  legend4->AddEntry((TObject*)0,"gluino mass","");
  for(int j=9; j>-1; j-=1){//for each mass point
    TFile infile4(std::string("hist")+mass[j]+"_"+size+"_"+position+".root");
    infile4.GetObject("Graph;3",betaEfficiency);
    if(j==9){
      betaEfficiency->Draw("AL3");
      betaEfficiency->SetMaximum(0.002);
      betaEfficiency->SetMinimum(0);
      betaEfficiency->SetTitle(";#beta;Absorption efficiency x acceptance");
      TAxis *xaxis0 = betaEfficiency->GetXaxis();
      xaxis0->SetLimits(0.1,1.);
      xaxis0->SetTitleSize(0.045);
      betaEfficiency->GetYaxis()->SetTitleSize(0.045);
    }
    
    else{betaEfficiency->Draw("L3 SAME");}
    legend4->AddEntry(betaEfficiency,mass[j]+std::string(" GeV"),"f");
    betaEfficiency->SetFillStyle(fillStyle[j]);
    betaEfficiency->SetFillColor(fillColour[j]);
    betaEfficiency->SetLineColor(lineColour[j]);

    //print the total eff x acc per mass point
    double sum = 0.;
    for(int k=2; k<betaEfficiency->GetN(); k++){ //start from k=2, which is 0.1 in beta
      sum += betaEfficiency->GetPointY(k);
      //std::cout<<"  the x value of point "<<k<<" is: "<<betaEfficiency->GetPointX(k)<<std::endl;
      //std::cout<<"  the y value of point "<<k<<" is: "<<betaEfficiency->GetPointY(k)<<std::endl;
    }
    std::cout<<"For mass "<<mass[j]<<", the sum of the y points of the Absorption efficiency x acceptance vs beta histogram is: "<<sum<<std::endl;
  }
  //write detector position and size
  TPaveText *pt = new TPaveText(.30,.85,.70,.93,"brNDC"); //use for position 0
  //TPaveText *pt = new TPaveText(.40,.85,.80,.93,"brNDC"); //use for position 1
  pt->AddText((std::string("Absorber position ")+position).c_str());
  pt->AddText((size+std::string("m x ")+size+"m x 2m").c_str());
  pt->Draw();
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  legend4->Draw();
  legend4->SetBorderSize(0);
  canvas4->Print((std::string("betaEfficiency_")+size+"_"+position+".pdf").c_str());
  
  return 0;
}
