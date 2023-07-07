///Modifying the the program to generate diffractive dijet events in ultra-peripheral
// Pb+p collisions. from  Author: Ilkka Helenius


// My Main94.cc to produce photo-nuclear dijet production set up 
//Modified from  test79.cc 



#include "Pythia8/Pythia.h"

using namespace Pythia8;

// Photon flux from lead-ions. Integrated over impact parameters > 2*r_Pb.
// Suitable for photo-nuclear processes but not for photon-photon.
// This should be considered as an experimental setup and used with caution.

class Nucleus2gamma2 : public PDF {

public:

  // Constructor.
  Nucleus2gamma2(int idBeamIn) : PDF(idBeamIn) {}

  // Update the photon flux.
  void xfUpdate(int , double x, double ) {

    // Minimum impact parameter (~2*radius) [fm].
    // double bmin = 2 * 6.636;
    double bmin = 6.636 + 0.7;

    // Charge of the nucleus.
    double z = 82.;

    // Per-nucleon mass for lead.
    double m2 = pow2(0.9314);
    double alphaEM = 0.007297353080;
    double hbarc = 0.197;
    double xi = x * sqrt(m2) * bmin / hbarc;
    double bK0 = besselK0(xi);
    double bK1 = besselK1(xi);
    double intB = xi * bK1 * bK0 - 0.5 * pow2(xi) * ( pow2(bK1) - pow2(bK0) );
    xgamma = 2. * alphaEM * pow2(z) / M_PI * intB;
  }


};

// The main program.

int main() {

  // Generator.
  Pythia pythia;
  Settings settings = pythia.settings;

  // Decrease the output.
  pythia.readString("Init:showChangedSettings = off");
  pythia.readString("Init:showChangedParticleData = off");
  pythia.readString("Next:numberCount = 100");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 1");


  // Beam settings.
  pythia.readString("Beams:frameType = 2");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eA = 4080."); /// Beam "A" will be the one to associate the photon flux
  pythia.readString("Beams:eB = 4080.");///Check energy****
  pythia.readString("PDF:beamA2gamma = on"); /// Enables photon sub-beam from beam particle "A"
  pythia.readString("PDF:Proton2gammaSet = 0"); /// this is to use an externally provided photon flux. we'll associate the lead flux defined above to the proton beam A, because there is not flux from heavy ions internally defined.  
  pythia.readString("PDF:beam2gammaApprox = 2"); /// Estimate optimized for ultraperipheral heavy-ion collisions for photon flux sampling. it is optimized for pPb collisions where the lead emits the photon. 
  pythia.readString("Photon:SampleQ2 = off"); /// because we're not using pronton flux,"Q2 integrated flux".
  
  
  // Set up the photon flux.
  PDFPtr photonFlux = make_shared<Nucleus2gamma2>(2212);
  pythia.setPhotonFluxPtr(photonFlux, 0);   ///this set up the photon-flux to beam "A", nothing to beam B (second entry)


  // Use optimized pT0ref for photon-hadron, a more detailed fit for this?
  pythia.readString("MultipartonInteractions:pT0Ref = 3.0");

  // Photoproduction and relevant hard processes.
  pythia.readString("Photon:ProcessType = 0");
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhotonParton:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 15.0");

// ************** for difractive events include the lines below ********** ///

// Example settings for hard diffraction. 
// Comment these out if after for inclusive jets.
//  pythia.readString("PDF:PomSet = 8");
//  pythia.readString("Diffraction:hardDiffSide = 1");
//  pythia.readString("SigmaDiffractive:PomFlux = 7");
// pythia.readString("Diffraction:doHard = on");

  //Comment: sampleType == 3 => PDF selection, sampleType == 4 => MPI selection.
  //this means: option 3 Generate an exclusive diffractive sample with no MPI. //
  //option4 Generate an exclusive diffractive sample with MPI. 

//  pythia.readString("Diffraction:sampleType = 3");

// **********************************************************************///

  // Initialize the generator.
  pythia.init();

  // Counters.
  int nDiffA = 0;
 // int nDiffB = 0;

  // Number of events.
  int nEvent = 1000;

  // Begin event loop. Skip if fails.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate next event.
    if (!pythia.next()) continue;

    // Calculate number of diffractive events in side B.
    if (pythia.info.isHardDiffractiveA()) ++nDiffA;
    //if (pythia.info.isHardDiffractiveB()) ++nDiffB;

    /// Analysis here...

  } // End of event loop.

  // Show statistics.
  pythia.stat();

  cout << "Number of diffractive events in side A = " << nDiffA << endl << endl;
  //cout << "Number of diffractive events in side B = " << nDiffB << endl << endl;

  // Done.
  return 0;
}
