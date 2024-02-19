///Modifying the the program to generate diffractive dijet events in ultra-peripheral
// Pb+p collisions. from  Author: Ilkka Helenius


// My Main94.cc to produce photo-nuclear dijet production set up 
//Modified from  test78.cc Ilkka Helenius


#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

#include <math.h>
#include <cmath>

#include "Pythia8/Pythia.h"

using namespace Pythia8;

using namespace std;

///Adding photon-flux def

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


int main() {
    
    Pythia pythia;
    
    // Decrease the output.
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberCount = 1000");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 1");  // default value
    pythia.readString("Next:numberShowEvent = 0");   //  default value

    // Beam parameters with heavy-ion optimized photon flux
    pythia.readString("Beams:frameType = 2");    // to identify two beams
    pythia.readString("Beams:idA = 2212"); //Proton(2212)  
    pythia.readString("Beams:idB = 2212"); //Proton (2212) 
    
   

    // p Pb energy.  energy per nucleon  
    pythia.readString("Beams:eA = 6500");  // Beam going +z. In this case proton 
    pythia.readString("Beams:eB = 2560");  // Beam going -z . Pb going left  
    
//  pythia.readString("PDF:beamA2gamma = on");
    pythia.readString("PDF:beamB2gamma = on"); // Enables photon sub-beam from beam particle B so the photon flux is associatetd to the lead. 
    pythia.readString("PDF:proton2gammaSet = 0");
    pythia.readString("PDF:beam2gammaApprox = 2");// Estimate optimized for ultraperipheral heavy-ion collisions for photon flux sampling. it is optimized for pPb collisions where the lead emits the photon.
    pythia.readString("Photon:sampleQ2 = off");/// This is because we're not using pronton flux,"Q2 integrated flux".
    
    
    
   // Set up the photon flux,
  
   PDFPtr photonFlux = make_shared<Nucleus2gamma2>(2212);
// pythia.setPhotonFluxPtr(photonFlux, 0);   ///this set up the photon-flux to beam "A", nothing to beam B (second entry)
   
   pythia.setPhotonFluxPtr(0, photonFlux); // this set up the photon-flux to beam "B" , so the lead in our case.
    
    

 //  pythia.readString("MultipartonInteractions:pT0Ref = 3.0"); do we need pT0ref?
    
     // Photoproduction and relevant hard processes.
    pythia.settings.mode("Photon:ProcessType", 0); // 0 is Mix of resolved-resolved, resolved- direct, direct-resolved,  direct- direct .   3 is the case direct-resolved. 

    pythia.readString("PhotonParton:all = on"); // all q, qbar
    
    pythia.readString("PartonLevel:MPI = off"); // 

  // Initialize the generator.
//    pythia.readString("Random:setSeed = on");
//    pythia.readString("Random:seed = 1");
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 10.");
    
   
/********* Example settings for hard diffraction. ***********************//////
// Comment these out if after for inclusive jets.
/*  
   pythia.readString("PDF:PomSet = 8");
   pythia.readString("Diffraction:hardDiffSide = 0"); // 1 Check for diffraction on side A only. 2 para B, and 0 Check for diffraction on both side
   pythia.readString("SigmaDiffractive:PomFlux = 7");
   pythia.readString("Diffraction:doHard = on"); //can be any hard process (e.g. jets)
    
 
 //Comment: sampleType == 3 => PDF selection, sampleType == 4 => MPI selection.
  //this means: option 3 Generate an exclusive diffractive sample with no MPI. //
  //option4 Generate an exclusive diffractive sample with MPI. 

  pythia.readString("Diffraction:sampleType = 3");
*/    
  
    
    
//    int nDiffA=0;
//    int nDiffB=0; 
    
  // Number of events.    
    int numEvent = 10000;
    
    pythia.init();
    
    // Create and open file for LHEF 3.0 output
    LHEF3FromPythia8 myLHEF3(&pythia.event, &pythia.info);
    // Open a file on which LHEF events should be stored, and write header
    myLHEF3.openLHEF("lhe_files_before_rm_headers/gamma_proton/gamma_proton_pPb_8p16TeV_test_file10.lhe");
    myLHEF3.setInit();
    myLHEF3.initLHEF();
    /* ---- ---- ---- ---- ---- ---- ---- ---- 
    ---- ---- ---- ---- ---- ---- ---- ---- */
            // Begin event loop. Generate event.
    for (int iEvent = 0; iEvent < numEvent; ++iEvent) {

        if (!pythia.next()) continue;
        // Store event info in the LHAEF3 object
        myLHEF3.setEvent();
        // Write out this event info on the file
        //myLHEF3.eventLHEF();
        
//       if (pythia.info.isHardDiffractiveA()) ++nDiffA;
//       if (pythia.info.isHardDiffractiveB()) ++nDiffB;
    
    }
    // Statistics: full printout.
    pythia.stat();
    // Write endtag. Overwrite initialization info with new cross sections.
    myLHEF3.closeLHEF(true);
    cout << "Finishing ... " << endl;
    
//    cout << "Number of diffractive events in side A = " << nDiffA << endl << endl;
//    cout << "Number of diffractive events in side B = " << nDiffB << endl << endl;
    
    // Done.
    return 0;
}
