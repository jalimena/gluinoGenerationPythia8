// main28.cc is a part of the PYTHIA event generator.
// Copyright (C) 2021 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: BSM; R-hadron; displaced vertex;

// Example of of R-hadron production.
// Several of the possibilities shown here, like displaced vertices,
// are extras that need not be used for the basic setup.

#include "Pythia8/Pythia.h"

//#include "TH1.h"
//#include "TFile.h"
//#include "TApplication.h" // ROOT, for interactive graphics.

using namespace Pythia8;

int main(int argc, char* argv[]) {

  // Key settings to be used in the main program.
  // nGluino = 0, 1, 2 give stop pair, single gluino or gluino pair.
  // JA: yes, we want 2 gluinos per event
  int nGluino  = 2;
  int nEvent   = 200;
  int nAbort   = 3;
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
  pythia.readString("SLHA:file = HSCP_gluino_300_SLHA.spc");

  // Further hacked file, to test R-parity violating gluino decay.
  //pythia.readString("SLHA:file = sps1aNarrowStopGluinoRPV.spc");

  // Allow R-hadron formation.
  pythia.readString("Rhadrons:allow = on");

  // If you want to do the decay separately later,
  // you need to switch off automatic decays.
  // JA: yes, we want r-hadrons to be created but for them not to decay
  pythia.readString("RHadrons:allowDecay = off");

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
  //TApplication theApp("hist", &argc, argv);

  // Create file on which histogram(s) can be saved.
  //TFile* outFile = new TFile("hist.root", "RECREATE");

  // Histograms.
  // JA: Jasmine, you can add more histograms here
  /*
  TH1F *nChargedH = new TH1F("nChargedH", "charged multiplicity", 100, -0.5, 799.5);
  TH1F *dndyChargedH = new TH1F("dndyChargedH", "dn/dy charged", 100, -10., 10.);
  TH1F *dndyRH = new TH1F("dndyRH", "dn/dy R-hadrons", 100, -5., 5.);
  TH1F *pTRH = new TH1F("pTRH", "pT R-hadrons", 100, 0., 1000.);
  TH1F *xRH = new TH1F("xRH", "p_RHadron / p_sparticle", 100, 0.9, 1.1);
  TH1F *mDiff = new TH1F("mDiff", "m(Rhadron) - m(sparticle)", 100, 0., 5.);
  TH1F *decVtx = new TH1F("decVtx", "R-hadron decay vertex (mm from origin)", 100, 0., 1000.);
  TH1F *eta = new TH1F("eta", "eta", 100, -5, 5);
  */

  Hist nChargedH("charged multiplicity", 100, -0.5, 799.5);
  Hist dndyChargedH("dn/dy charged", 100, -10., 10.);
  Hist dndyRH("dn/dy R-hadrons", 100, -5., 5.);
  Hist pTRH("pT R-hadrons", 100, 0., 1000.);
  Hist xRH("p_RHadron / p_sparticle", 100, 0.9, 1.1);
  Hist mDiff("m(Rhadron) - m(sparticle)", 100, 0., 5.);
  Hist decVtx("R-hadron decay vertex (mm from origin)", 100, 0., 1000.);
  Hist eta("eta", 100, -5, 5);

  // R-hadron flavour composition.
  map<int, int> flavours;

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
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
          dndyChargedH.fill( event[i].y() );
        }
      }
    }
    nChargedH.fill( nCharged );

    // Loop over final R-hadrons in the event: kinematic distribution
    for (int i = 0; i < event.size(); ++i) {
      int idAbs = event[i].idAbs();
      if (idAbs > 1000100 && idAbs < 2000000 && idAbs != 1009002) {
        ++flavours[ event[i].id() ];
        dndyRH.fill( event[i].y() );
        pTRH.fill( event[i].pT() );
        eta.fill( event[i].eta() );
        // Trace back to mother; compare momenta and masses.
        int iMother = i;
        while( event[iMother].statusAbs() > 100)
          iMother = event[iMother].mother1();
        double xFrac = event[i].pAbs() / event[iMother].pAbs();
        xRH.fill( xFrac);
        double mShift = event[i].m() - event[iMother].m();
        mDiff.fill( mShift );
        // Separation of R-hadron decay vertex from origin.
        // Don't be fooled by pAbs(); it gives the three-vector length
        // of any Vec4, also one representing spatial coordinates.
        double dist = event[i].vDec().pAbs();
        decVtx.fill( dist);

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

  // Final statistics, flavour composition and histogram output.
  pythia.stat();
  cout << "\n Composition of produced R-hadrons \n    code            "
       << "name   times " << endl;
  for (map<int, int>::iterator flavNow = flavours.begin();
    flavNow != flavours.end(); ++flavNow)  cout << setw(8)
    << flavNow->first << setw(16) << pythia.particleData.name(flavNow->first)
    << setw(8) << flavNow->second << endl;
  cout << nChargedH << dndyChargedH << dndyRH << pTRH << xRH << mDiff
       << decVtx << eta;

  //write histograms to output file and close it
  /*
  nChargedH->Write();
  dndyChargedH->Write();
  dndyRH->Write();
  pTRH->Write();
  xRH->Write();
  mDiff->Write();
  decVtx->Write();
  delete outFile;
  */

  // Done.
  return 0;
}