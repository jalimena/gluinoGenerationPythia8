#include "TROOT.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <iostream>
#include "tdrstyle.C"

int energyDeposited() {
  setTDRStyle();
  int deltaMass[5]={10,20,50,100,200};//GeV
  const char* graph[5]={"energyDeposited10","energyDeposited20","energyDeposited50","energyDeposited100","energyDeposited200"};//GeV
  double lambdaArgon=14;//cm radiation length
  double lambdaCopper=1.436;
  double lambdaZinc=1.742;
  double lambdaBrass=2*lambdaCopper/3+lambdaZinc/3;
  double rodWidth=1;//cm
  int vectLength=320;//max rod spacing is vecLength/20
  double rodSpacing[vectLength];
  for(int i=0; i<vectLength; i++){rodSpacing[i]=i/20.;}//fill rod spacing vector
  double areaBrass=rodWidth*rodWidth;
  double effAreaBrass=areaBrass/lambdaBrass;
  double effAreaArgon;
  double areaArgon;
  double energyDeposited[vectLength];
  double ionPairs[vectLength];
  //create energy deposited graph
  auto canvas=new TCanvas("canvas","canvas");
  gStyle->SetPalette(kRainBow);
  auto* legend = new TLegend(0.7,0.4,0.95,0.15);
  legend->AddEntry((TObject*)0,"mass splitting","");
  legend->SetBorderSize(0);
  for(int i=4; i>-1; i-=1){
    for(int j=0; j<vectLength; j++){
      areaArgon=(rodWidth+rodSpacing[j])*(rodWidth+rodSpacing[j])-areaBrass;
      effAreaArgon=areaArgon/lambdaArgon;
      energyDeposited[j]=deltaMass[i]*effAreaArgon/(effAreaBrass+effAreaArgon);
    }
    TGraph *graph[i];
    graph[i] = new TGraph(vectLength,rodSpacing,energyDeposited);
    if(i==4){
      graph[i]->Draw("AL PLC");
      graph[i]->SetTitle(";Rod spacing [cm];Energy deposited [GeV]");
      graph[i]->SetMaximum(300);
      graph[i]->SetMinimum(0);
      gPad->SetLogy();
    }
    else{graph[i]->Draw("L SAME PLC");}
    legend->AddEntry(graph[i],deltaMass[i]+std::string(" GeV"),"l");

  }
  legend->Draw();
  canvas->Print("energyDeposited.pdf");
  //repeat for ions deposited
  auto canvas1=new TCanvas("canvas1","canvas1");
  gStyle->SetPalette(kRainBow);
  auto* legend1 = new TLegend(0.7,0.4,0.95,0.15);
  legend1->AddEntry((TObject*)0,"mass splitting","");
  legend1->SetBorderSize(0);
  for(int i=4; i>-1; i-=1){
    for(int j=0; j<vectLength; j++){
      areaArgon=(rodWidth+rodSpacing[j])*(rodWidth+rodSpacing[j])-areaBrass;
      effAreaArgon=areaArgon/lambdaArgon;
      ionPairs[j]=deltaMass[i]*effAreaArgon/(effAreaBrass+effAreaArgon)*1000000000/23.6;
    }
    TGraph *graph[i];
    graph[i] = new TGraph(vectLength,rodSpacing,ionPairs);
    if(i==4){
      graph[i]->Draw("AL PLC");
      graph[i]->SetTitle(";Rod spacing [cm];Number of ion pairs deposited");
      graph[i]->SetMaximum(1.3*pow(10,10));
      graph[i]->SetMinimum(0);
      gPad->SetLogy();
    }
    else{graph[i]->Draw("L SAME PLC");}
    legend1->AddEntry(graph[i],deltaMass[i]+std::string(" GeV"),"l");
  }
  legend1->Draw();
  canvas1->Print("ionPairs.pdf");

  return 0;
}
