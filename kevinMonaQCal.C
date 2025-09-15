// MoNA light calibration using cosmics
// Kevin Eisenberg, September 2025
// For use in FRIB e23033

#include <fstream>
#include "TStopwatch.h"

void kevinMonaQCal(int x, int y, int z, bool graphing = true) {
		gROOT->SetBatch(!graphing);
		
    // random run
    TFile *f = new TFile("run-0034-00--raw.root");
    TTree *t = (TTree*)f->Get("t");

    // output file
    ofstream csvFile("calibratedPeaks.csv");
    csvFile << "x,y,z,peak\n";

    UShort_t light[18][16][2];
    t->SetBranchAddress("light[18][16][2]", light);
    Long64_t nEntries = t->GetEntries(); 
	
    TStopwatch timer; timer.Start();

    // declare hists
    TH1F *h1 = new TH1F(Form("raw[%d][%d][%d]",x,y,z),"Raw charge;Channel;Counts",275,0,4096);
    TH1F *h2 = new TH1F(Form("cal[%d][%d][%d]",x,y,z),"Calibrated charge;Energy (MeV);Counts",130,0,80);

    // raw data
    for (Long64_t i=0; i<nEntries; i++) { 
    	t->GetEntry(i); 
	h1->Fill(light[x][y][z]); //x<=18 (but only has data up to 8), y<=16, z 1 or 2 (L/R)
    }
    
    // find peaks
    TSpectrum *spec = new TSpectrum(5);
    Int_t nRawPeaks = spec->Search(h1,2,"",0.0005);

    Double_t *xRawPeaks = spec->GetPositionX();
    cout << "Calibrating bar " << x << ", " << y << ", " << z << endl;

    if (nRawPeaks >= 2) {
	double x1 = xRawPeaks[0];
	double x2 = xRawPeaks[1];
	std::cout << "Raw peaks: " << x1 << " , " << x2 << std::endl;

	// linear transformation params
	double a = 20.5/(x2-x1);
	double b = -20.5*x1/(x2-x1);

	// calibrated hist
	for (Long64_t i=0; i<nEntries; i++) {
	    t->GetEntry(i);
	    double calEntry = a * light[x][y][z] + b;
	    h2->Fill(calEntry);
	}
    }
	
    // output: location of muon peak in calibrated spectrum
    TSpectrum *calSpec = new TSpectrum(3);
    Int_t nCalPeaks = calSpec->Search(h2,2,"",0.1);
    Double_t *xPeaks = calSpec->GetPositionX();

    if (nCalPeaks > 0) {
	timer.Stop();
	double realTime = timer.RealTime();
	std::cout << "Time: " << realTime << "s. Calibrated peak: " << xPeaks[0] << std::endl;
	csvFile << x << "," << y << "," << z << "," << xPeaks[0] << "\n";
    }

    if (graphing) {
    TCanvas *c1 = new TCanvas("c1","raw",800,600);
    h1->Draw();
    gPad->SetLogy();

    TCanvas *c2 = new TCanvas("c2","cal",800,600);
    h2->Draw("HIST");
    gPad->SetLogy();

    }
    csvFile.close();
}
