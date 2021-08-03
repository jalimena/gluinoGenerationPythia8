// main28.cc is a part of the PYTHIA event generator.
// Copyright (C) 2021 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: BSM; R-hadron; displaced vertex;

// Example of of R-hadron production.
// Several of the possibilities shown here, like displaced vertices,
// are extras that need not be used for the basic setup.

#include "Pythia8/Pythia.h"
#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TApplication.h" // ROOT, for interactive graphics.
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
using namespace Pythia8;

int main(int argc, char* argv[]) {
  if(argc==1){
    cout<<"input gluino mass e.g. 'nohup ./gluinoGeneration 5 >& out &'";
    return 0;
  }
  if(atoi(argv[3])!=0 && atoi(argv[3])!=1){
    cout<<"input 3rd argument 0 for detector above collision point or 1 for detector diagonal to collsion point e.g 'nohup ./gluinoGeneration 5 2 0>& out &'";
    return 0;
  }
  // Key settings to be used in the main program.
  // nGluino = 0, 1, 2 give stop pair, single gluino or gluino pair.
  // JA: yes, we want 2 gluinos per event
  int nGluino  = 2;
  int nEvent   = 100000;
  int nAbort   = 10000;
  int nList    = 0;
  double eCM   = 13000.; //JA: changed to 13 TeV CME

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Set up beams: p p is default so only need set energy.
  pythia.settings.parm("Beams:eCM", eCM);

  // Squark pair: use stop-antistop as example.
  if (nGluino == 0) {
    pythia.readString("SUSY:gg2squarkantisquark = on");
    pythia.readString("SUSY:idA = 1000006");
    pythia.readString("SUSY:idB = 1000006");
  // Squark-gluino pair: also supersymmetric u has been made long-lived.
  // Stop does not work since then one would need inoming top PDF.
  // Nevertheless R-hadrons are numbered/named as if containing a stop.
  } else if (nGluino == 1) {
    pythia.readString("SUSY:qg2squarkgluino  = on");
    pythia.readString("SUSY:idA = 1000002");
    pythia.readString("RHadrons:idStop = 1000002");
    pythia.readString("SUSY:idB = 1000021");
  // Gluino pair.
  } else {
    pythia.readString("SUSY:gg2gluinogluino  = on");
    pythia.readString("SUSY:qqbar2gluinogluino  = on"); //JA: added this production mode as well
  }

  // Use hacked sps1a file, with stop (+su) and gluino made long-lived.
  // This is based on the width being less than 0.2 GeV by default.
  //pythia.readString("SLHA:file = sps1aNarrowStopGluino.spc");
  // JA: put in HSCP gluino SLHA file with mass of 300 GeV
  pythia.readString(std::string("SLHA:file = HSCP_gluino_")+argv[1]+"_SLHA.spc");

  // Further hacked file, to test R-parity violating gluino decay.
  //pythia.readString("SLHA:file = sps1aNarrowStopGluinoRPV.spc");

  // Allow R-hadron formation.
  pythia.readString("Rhadrons:allow = on");

  // If you want to do the decay separately later,
  // you need to switch off automatic decays.
  // JA: yes, we want r-hadrons to be created but for them not to decay
  pythia.readString("RHadrons:allowDecay = off");

  pythia.readString("RHadrons:setMasses = on");

  // Fraction of gluinoballs.
  // JA: same gluinoball fraction as HSCP search
  pythia.readString("RHadrons:probGluinoball = 0.1");

  // Switch off key components.
  //pythia.readString("PartonLevel:MPI = off");
  //pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("HadronLevel:Hadronize = off");

  // Allow the R-hadrons to have secondary vertices: set c*tau in mm.
  // Note that width and lifetime can be set independently.
  // (Nonzero small widths are needed e.g. to select branching ratios.)
  // JA: commented out tau settings. we want to make the BSM particles completely stable at this stage
  //pythia.readString("1000002:tau0 = 200.");
  //pythia.readString("1000006:tau0 = 250.");
  //pythia.readString("1000021:tau0 = 300.");
  pythia.readString("1000021:tau0 = 1000000"); //JA: set gluino ctau to be 1km

  // Checks. Optionally relax E-p-conservation.
  pythia.readString("Check:nErrList = 2");
  //pythia.readString("Check:epTolErr = 2e-3");

  // Possibility to switch off particle data and event listings.
  // Also to shop location of displaced vertices.
  pythia.readString("Init:showChangedSettings = on");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberShowInfo = 1");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:showScaleAndVertex = on");

  // Initialize.
  pythia.init();

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Create file on which histogram(s) can be saved.
  std::ostringstream fileNameStream("hist");
  fileNameStream << "hist"<< argv[1]<<"_"<< argv[2]<<"_"<<argv[3] << ".root";
  std::string fileName=fileNameStream.str();
  TFile* outFile = new TFile(fileName.c_str(), "RECREATE");

  //Create bin widths for log scale
  const int nbins=25;
  Double_t xbins[nbins+1];
  double dx = 3.5/nbins;
  for (int i=0;i<=nbins;i++) {
    xbins[i] = exp(log(10)*i*dx);
  }
  // Binning for beta histograms
  double betaMin=0.;
  double betaMax=1.;
  double betaBinWidth=(betaMax-betaMin)/nbins;

  // Histograms.
  // JA: Jasmine, you can add more histograms here
  TH1F *nChargedH = new TH1F("nChargedH", "charged multiplicity", 100, -0.5, 799.5);
  TH1F *chargeRH = new TH1F("chargeRH", "charge of R-hadrons", 61, -10, 10);
  TH1F *dndyChargedH = new TH1F("dndyChargedH", "dn/dy charged", 100, -10., 10.);
  TH1F *dndyRH = new TH1F("dndyRH", "dn/dy R-hadrons", 100, -5., 5.);
  TH1F *pTRH = new TH1F("pTRH", "pT of R-hadrons",nbins,xbins);//log scale
  TH1F *etaRH = new TH1F("etaRH", "eta of R-hadrons", 50, -5, 5);
  TH1F *thetaRH = new TH1F("thetaRH", "theta of R-hadrons", 50, -0, 3.2);
  TH1F *xRH = new TH1F("xRH", "p_RHadron / p_sparticle", 100, 0.9, 1.1);
  TH1F *mDiff = new TH1F("mDiff", "m(Rhadron) - m(sparticle)", 100, 0., 5.);
  TH1F *decVtx = new TH1F("decVtx", "R-hadron decay vertex (mm from origin)", 100, 0., 1000.);
  TH1F *phiRH = new TH1F("phiRH", "Phi of R-hadrons", 20, -3.2, 3.2);
  //TH1F *phiPosRH = new TH1F("phiPosRH", "Phi (position) of R-hadrons", 20, -3.2, 3.2);
  //TH1F *phiVelRH = new TH1F("phiVelRH", "Phi (velocity) of R-hadrons", 20, -3.2, 3.2);
  TH1F *betaRH = new TH1F("betaRH", "Beta of R-hadrons", nbins, betaMin,betaMax );
  TH1F *acceptedPTRH = new TH1F("acceptedPTRH", "pT of R-hadrons that hit detector",nbins,xbins);//log scale
  TH1F *acceptedBetaRH = new TH1F("acceptedBetaRH", "Beta of R-hadrons that hit detector", nbins, betaMin, betaMax);
  TH1F *distanceRH = new TH1F("distanceRH", "Distance R-hadron travels through detector", 30, 0.,3 );
  TH1F *absorbedBetaRH = new TH1F("absorbedBetaRH", "Beta of R-hadrons that are absorbed by detector", nbins, betaMin, betaMax);
  TH1F *absorbedPTRH = new TH1F("absorbedPTRH", "pT of R-hadrons that are absorbed by detector",nbins,xbins);//log scale
  
  // Set detector position and size 
  double magneticField1=3.8;
  double magneticField2=1;
  double fieldRadius=2.95; //radius of magetic field
  double fieldLength=6.45;
  double detectorWidth;//side lengths of detector face facing detector 
  double detectorDepth=2;
  double cmsRadius=7.38; 
  double cmsLength=10.91;
  double detectorHeight;//height of detector above beamline for the case of the detector in front of CMS
  double detectorDistance;//distance of detector infront of collision point for the case of the detector in front of CMS
  int detectorPosition;//0 if detctor above collision point, 1 for detctor diagonal to collision point
  if(argc==2){ //if detector size and position is not specified from command line
    detectorWidth=2; //default values
    detectorPosition=0;
  }
  else{ 
    detectorWidth=atof(argv[2]);
    detectorPosition=atoi(argv[3]);
    if(detectorPosition==1){
         if(argc==5||argc==4){//if height or distance isn't set use default values
	detectorHeight=3;
	detectorDistance=cmsLength;
      }
      else{
	detectorHeight=atof(argv[4]);
	detectorDistance=atof(argv[5]);
      }
    }
  }
  //cout<<detectorPosition;
  
  //double acceptanceAngle=atan(detectorWidth/(2*cmsRadius)); //max phi for a particle to hit detector
  double pi = 3.14159265358979323;
  //  double maxAcceptanceEta=log(tan((pi/2+acceptanceAngle)/2)); //max eta for a particle to hit detector
  //  double minAcceptanceEta=log(tan((pi/2-acceptanceAngle)/2)); //min eta for a particle to hit detector
  //  double maxAcceptanceTheta=pi/2+acceptanceAngle; //max theta for a particle to hit detector
  //  double minAcceptanceTheta=pi/2-acceptanceAngle; //min theta for a particle to hit detector
  // R-hadron flavour composition.
  map<int, int> flavours;

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      cout<< "Event aborted\n";
      if (++iAbort < nAbort) continue;
      //      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }
    // Loop over final charged particles in the event.
    // The R-hadrons may not yet have decayed here.

    int nCharged = 0;
    Vec4 pSum;
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal()) {
        pSum += event[i].p();
        if (event[i].isCharged()) {
	  ++nCharged;
          dndyChargedH->Fill( event[i].y() );
        }
      }
    }
    nChargedH->Fill( nCharged );

    // Loop over final R-hadrons in the event: kinematic distribution
    for (int i = 0; i < event.size(); ++i) {
      int idAbs = event[i].idAbs();
      if (idAbs > 1000100 && idAbs < 2000000 && idAbs != 1009002) {
        ++flavours[ event[i].id() ];
        dndyRH->Fill( event[i].y() );
        pTRH->Fill( event[i].pT() );
	double eta = event[i].eta();
	etaRH->Fill(eta);
	double theta = 2*atan(exp(-eta));
	thetaRH->Fill(theta);
        double phi=event[i].phi();
	phiRH->Fill(phi );
	double charge = event[i].charge();
	chargeRH->Fill(charge);
	double p2= event[i].pAbs2();//momentum squared
	double m2= event[i].m()*event[i].m();//mass squared
	double beta=sqrt(p2/(p2+m2));
	betaRH->Fill(beta);
	
	//double phiPos;
	//double phiVel;
	double betaT=event[i].pT()/sqrt(p2+m2);//transverse velocity
	double betaZ=sqrt(beta*beta-betaT*betaT);//Z velocity

	//detector above collision point
	if(detectorPosition==0){
	  double lZ1;//z distance travelled in inner magnetic field
	  double lZ2;//z distance travelled in residual magnetic field
	
	  //For charged events calculate new velocity and position angles after magnetic field
	  if(event[i].isCharged()){
	    double curvatureRadius1=event[i].m()*betaT*1000000000/(299792458*abs(charge)*magneticField1);// radius of curvature due to inner magenetic field
	    double curvatureRadius2=event[i].m()*betaT*1000000000/(299792458*abs(charge)*magneticField2);// radius of curvature due to outer magenetic field
	    double deltaPhiVel1=acos(1-fieldRadius*fieldRadius/(2*curvatureRadius1*curvatureRadius1));// change in angle due to inner field
	    double deltaPhiVel2=acos(1-(cmsRadius-fieldRadius)*(cmsRadius-fieldRadius)/(2*curvatureRadius2*curvatureRadius2));// change in angle due to residual field
	 
	    /*phiVel=phi+(charge/abs(charge))*deltaPhiVel;//phi of velocity vector after magnetic field
	      if(phiVel>pi){phiVel-=2*pi;}//loop phi to keep in range [-pi,pi]
	      if(phiVel<-pi){phiVel+=2*pi;}
	      double deltaPhiPos=asin(detectorWidth/cmsRadius*sin(pi-deltaPhiVel/2));
	      phiPos=phi+(charge/abs(charge))*(deltaPhiVel-deltaPhiPos);//phi of position vector when particle reached detector
	      if(phiPos>pi){phiVel-=2*pi;}
	      if(phiPos<-pi){phiVel+=2*pi;}*/
	    //	  double lT=curvatureRadius*deltaPhiVel; //distance travelled in transverse direction in magnetic field
	    lZ1=curvatureRadius1*deltaPhiVel1*betaZ/betaT; //distance travelled in z direction in inner magnetic field
	    lZ2=curvatureRadius2*deltaPhiVel2*betaZ/betaT;//distance travelled in Z direction in residual magnetic field
	  }
	  else{
	    //phiPos=phi;//if particle not charged angles are not affected by magnetic field
	    //phiVel=phi;
	    lZ1=fieldRadius*betaZ/betaT;
	    lZ2=(cmsRadius-fieldRadius)*betaZ/betaT;
	  }
	  //phiPosRH->Fill(phiPos);
	  //phiVelRH->Fill(phiVel);
	  //if(thetaVel>pi/2){thetaVel=-thetaVel+pi;}

	  //Check if event would hit detector
	  if(lZ1<fieldLength){ //particle leaves magnetic field through curved face
	    if(lZ1+lZ2<detectorWidth/2 && lZ1+lZ2>-detectorWidth/2){  
	      double radius=sqrt((lZ1+lZ2)*(lZ1*lZ2)+(cmsRadius*cmsRadius));
	      double phiMax=asin(detectorWidth/(2*radius));
	      if(phi>=-phiMax && phi<phiMax){
		acceptedPTRH->Fill( event[i].pT() );
		acceptedBetaRH->Fill(sqrt(p2/(p2+m2)) );
	      //Find distance travelled through detector
	      /*double distX=(detectorWidth/2-cmsRadius*sin(abs(phiPos)));
		double distZ=(detectorWidth/2-abs(lZ)-abs(lZMag));
		double distY=min(distZ*tan(thetaVel),distX/tan(abs(phiVel)));
		distY=min(distY,detectorDepth);
		double dist=sqrt(distZ*distZ+distY*distY+distX*distX);
		distanceRH->Fill(dist);*/
	      }
	    }
	  }
	}
	//detector diagonal to collsion point
	if(detectorPosition==1){
	  double lT1;//transverse distance travelled in inner magnetic field
	  double lT2;//transverse distance travelled in outer magnetic field
	  double angle1;//angle swept by particle path in inner magnetic field
	  double angle2;//angle swept by particle in residual magentic field
	  //For charged events calculate new velocity and position angles after magnetic field
	  if(event[i].isCharged()){//find transverse distance travelled if particle is charged
	    double curvatureRadius1=event[i].m()*betaT*1000000000/(299792458*abs(charge)*magneticField1);// radius of curvature due to inner magenetic field
	    double curvatureRadius2=event[i].m()*betaT*1000000000/(299792458*abs(charge)*magneticField2);// radius of curvature due to outer magenetic field
	    angle1=betaT*fieldLength/(betaZ*curvatureRadius1);
	    lT1=curvatureRadius1*sqrt(2*(1-cos(angle1)));
	    angle2=betaT*fieldLength/(betaZ*curvatureRadius2);
	    lT2=curvatureRadius2*sqrt(2*(1-cos(angle2)));
	  }
	  else{//find transverse distance travelled if particle not charged
	    lT1=fieldLength*betaT/betaZ;
	    lT2=(cmsLength-fieldLength)*betaT/betaZ;
	    angle1=0;
	    angle2=0;
	  }
	  //Check if event would hit detector
	  if(lT1<fieldRadius){ //particle leaves magnetic field through flat face
	    double lTTotal=sqrt(lT1*lT1+lT2*lT2-2*lT1*lT2*cos(pi-(angle1+angle2)/2));//total distance travelled in transverse direction before reaching detector
	    if(lTTotal>=detectorHeight-detectorDepth*tan(theta) && lTTotal<detectorHeight+detectorWidth){//particle has transverse displacement to hit detector
	      double phiMax=asin(detectorWidth/(2*lTTotal));
	      if(phi>=-phiMax && phi<phiMax){//particle has phi to hit detector
		acceptedPTRH->Fill( event[i].pT() );
		acceptedBetaRH->Fill(sqrt(p2/(p2+m2)) );
	      }
	    }
	  }
	}

        // Trace back to mother; compare momenta and masses.
        int iMother = i;
        while( event[iMother].statusAbs() > 100)
          iMother = event[iMother].mother1();
        double xFrac = event[i].pAbs() / event[iMother].pAbs();
        xRH->Fill( xFrac);
        double mShift = event[i].m() - event[iMother].m();
        mDiff->Fill( mShift );
        // Separation of R-hadron decay vertex from origin.
        // Don't be fooled by pAbs(); it gives the three-vector length
        // of any Vec4, also one representing spatial coordinates.
        double dist = event[i].vDec().pAbs();
        decVtx->Fill( dist);

        // This is a place where you could allow a R-hadron shift of
        // identity, momentum and decay vertex to allow for detector effects.
        // Identity not illustrated here; requires a change of mass as well.
        // Toy model: assume an exponential energy loss, < > = 1 GeV,
        // but at most half of kinetic energy. Unchanged direction.
        // Note that event will no longer conserve energy and momentum.
	// JA: is this what we want?
        double eLossAvg = 1.;
        double eLoss = 0.;
        do { eLoss = eLossAvg * pythia.rndm.exp(); }
        while (eLoss > 0.5 * (event[i].e() - event[i].m()));
        double eNew = event[i].e() - eLoss;
        Vec4   pNew = event[i].p() * sqrt( pow2(eNew) - pow2(event[i].m()) )
                    / event[i].pAbs();
        pNew.e( eNew);
        event[i].p( pNew);
        // The decay vertex will be calculated based on the production vertex,
        // the proper lifetime tau and the NEW four-momentum, rather than
        // e.g. some average momentum, if you do not set it by hand.
        // This commented-out piece illustrates brute-force setting,
        // but you should provide real numbers from some tracking program.
        // With tau = 0 the decay is right at the chosen point.
        //event[i].tau( 0.);
        //event[i].vProd( 132., 155., 233., 177.);

      // End of loop over final R-hadrons.
      }
    }

    // If you have set R-hadrons stable above,
    // you can still force them to decay at this stage.
    // JA: comment out, we don't want the r-hadrons to decay in pythia
    //pythia.forceRHadronDecays();
    if (iEvent < nList) pythia.event.list(true);

  // End of event loop.
  }

  //Read absorbtion values from files
  string betaLine;
  double betaCentre[320];
  string betaAbsLine;
  double betaAbs[320];
  string pTLine;
  double pTCentre[320];
  string pTAbsLine;
  double pTAbs[320];
  int nBeta;//number of values in beta file
  int nPT;//number of values in pt files
  ifstream betaFile (std::string("totalabs_beta_")+argv[1]+".7GeV.txt");
  if (betaFile.is_open())
    {
      int i=0;
      while ( getline (betaFile,betaLine) )
	{
	  betaCentre[i]=stod(betaLine);
	  i++;
	}
      betaFile.close();
    }
  else {cout << "Unable to open beta file"; return 0;}
  ifstream betaAbsFile (std::string("totalabs_beta_")+argv[1]+".7GeV_abs.txt");
  if (betaAbsFile.is_open())
    {
      int i=0;
      while ( getline (betaAbsFile,betaAbsLine) )
	{
	  betaAbs[i]=stod(betaAbsLine);
	  i++;
	}
      nBeta=i;
      betaAbsFile.close();
    }
  else {cout << "Unable to open beta abs file"; return 0;}
  ifstream pTFile (std::string("totalabs_energy_")+argv[1]+".7GeV.txt");
  if (pTFile.is_open())
    {
      int i=0;
      while ( getline (pTFile,pTLine) )
	{
	  pTCentre[i]=stod(pTLine);
	  i++;
	}
      pTFile.close();
    }
  else {cout << "Unable to open energy file"; return 0;}
  ifstream pTAbsFile (std::string("totalabs_energy_")+argv[1]+".7GeV_abs.txt");
  if (pTAbsFile.is_open())
    {
      int i=0;
      while ( getline (pTAbsFile,pTAbsLine) )
	{
	  pTAbs[i]=stod(pTAbsLine);
	  i++;
	}
      nPT=i;
      pTAbsFile.close();
    }
  else {cout << "Unable to open energy abs file"; return 0;}

  //plot graph of absorbtion efficiency
  TGraph* betaAbsorbed = new TGraph(nBeta,betaCentre,betaAbs);
  betaAbsorbed->SetTitle("Fraction of R-hardrons absorbed by detector;beta;Absorbtion efficiency");
  betaAbsorbed->Draw("AC");
  betaAbsorbed->Write();
  TGraph* pTAbsorbed = new TGraph(nPT,pTCentre,pTAbs);
  pTAbsorbed->SetTitle("Fraction of R-hardrons absorbed by detector;Transverse Momentum [GeV];Absorbtion efficiency");
  pTAbsorbed->Draw("AC");
  pTAbsorbed->Write();

  //Find average absorbtion 
  double nbinsWidth=(betaMax-betaMin)/nbins;
  double averageBetaAbs[nbins];
  double averagePTAbs[nbins];
  double betaEfficiency[nbins];
  double betaErrorUp[nbins];
  double betaErrorLow[nbins];
  double pTEfficiency[nbins];
  double pTErrorUp[nbins];
  double pTErrorLow[nbins];
  for(int i=0; i<nbins; i++){//for all bins in betaRH hist
    int j=0;
    averageBetaAbs[i]=0;
    for(int k=0; k<nBeta; k++){//for all values in absorbtion files
      if(betaCentre[k]>=betaMin+nbinsWidth*i && betaCentre[i]<betaMin+nbinsWidth*(i+1)){
	j++;
	averageBetaAbs[i]+=betaAbs[k];//sum all absorbtion values within each betaRH bin
      }
    }
    if(j==0){averageBetaAbs[i]=0;}
    else{averageBetaAbs[i]/=j;}//find average absorption efficiency for each betaRH bin
  }
  for(int i=0; i<nbins; i++){//for all bins in pTRH hist
    int j=0;
    averagePTAbs[i]=0;
    for(int k=0; k<nPT; k++){//for all values of absorbtion files
      if(pTCentre[k]>=xbins[i] && pTCentre[i]<xbins[i+1]){
	j++;
	averagePTAbs[i]+=pTAbs[k];//summ all absorbtion values within each pTRH bin
      }
    }
    if(j==0){averagePTAbs[i]=0;}
    else{averagePTAbs[i]/=j;}//find average absorption efficiency for each pTRH bin
  }

  //Make efficiency plot
  TEfficiency* betaAcceptedEff = 0;
  if(TEfficiency::CheckConsistency(*acceptedBetaRH,*betaRH))
    {
      betaAcceptedEff = new TEfficiency(*acceptedBetaRH,*betaRH);
      betaAcceptedEff->Write();
      for(int i=0;i<nbins;i++){
	betaEfficiency[i]=betaAcceptedEff->GetEfficiency(i+1);//get efficiency values
	betaErrorUp[i]=betaAcceptedEff->GetEfficiencyErrorUp(i+1);
	betaErrorLow[i]=betaAcceptedEff->GetEfficiencyErrorLow(i+1);
	betaEfficiency[i]*=averageBetaAbs[i];//get combine absorbtion efficiency and angular acceptance
	betaErrorLow[i]*=averageBetaAbs[i];
	betaErrorUp[i]*=averageBetaAbs[i];
      }
    }
  TEfficiency* pTAcceptedEff = 0;
  if(TEfficiency::CheckConsistency(*acceptedPTRH,*pTRH))
    {
      pTAcceptedEff = new TEfficiency(*acceptedPTRH,*pTRH);
      pTAcceptedEff->Write();
      for(int i=0;i<nbins;i++){
	pTEfficiency[i]=pTAcceptedEff->GetEfficiency(i+1);//get efficiency values
	pTErrorUp[i]=pTAcceptedEff->GetEfficiencyErrorUp(i+1);
	pTErrorLow[i]=pTAcceptedEff->GetEfficiencyErrorLow(i+1);
	pTEfficiency[i]*=averagePTAbs[i];//get combine absorbtion efficiency and angular acceptance
	pTErrorLow[i]*=averagePTAbs[i];
	pTErrorUp[i]*=averagePTAbs[i];
      }
    }

  //Plot absorption graphs
  double beta[nbins];//beta bin centres for x-axis
  double pT[nbins];//pT bin centres for x-axis
  for(int i=0;i<nbins;i++){
    pT[i]=(xbins[i]+xbins[i+1])/2;
    if(i==0){beta[i]=betaBinWidth/2;}
    else{beta[i]=beta[i-1]+betaBinWidth;}
  }
  auto betaEff = new TGraphAsymmErrors(nbins,beta,betaEfficiency,0,0,betaErrorLow,betaErrorUp);
  betaEff->SetTitle("Fraction of R-hardrons absorbed by detector;beta;Absorbtion efficiency");
  betaEff->Draw("AC");
  betaEff->Write();
  auto pTEff = new TGraphAsymmErrors(nbins,pT,pTEfficiency,0,0,pTErrorLow,pTErrorUp);
  pTEff->SetTitle("Fraction of R-hardrons absorbed by detector;Transverse Momentum [GeV];Absorbtion efficiency");
  pTEff->Draw("AC");
  pTEff->Write();

  // Final statistics, flavour composition and histogram output.
  pythia.stat();
  cout << "\n Composition of produced R-hadrons \n    code            "
       << "name   times " << endl;
  for (map<int, int>::iterator flavNow = flavours.begin();
    flavNow != flavours.end(); ++flavNow)  cout << setw(8)
    << flavNow->first << setw(16) << pythia.particleData.name(flavNow->first)
    << setw(8) << flavNow->second << endl;
  
  //write histograms to output file and close it
  nChargedH->Write();
  dndyChargedH->Write();
  chargeRH->Write();
  dndyRH->Write();
  pTRH->Write();
  etaRH->Write();
  thetaRH->Write();
  xRH->Write();
  mDiff->Write();
  decVtx->Write();
  phiRH->Write();
  //phiPosRH->Write();
  //phiVelRH->Write();
  betaRH->Write();
  acceptedBetaRH->Write();
  acceptedPTRH->Write();
  absorbedBetaRH->Write();
  absorbedPTRH->Write();
  distanceRH->Write();
  
  delete outFile;

  // Done.
  return 0;
}
