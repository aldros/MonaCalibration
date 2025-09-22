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
    TSpectrum *spec = new TSpectrum(4);
    Int_t nRawPeaks = spec->Search(h1,2,"",0.0001);

    Double_t *xRawPeaks = spec->GetPositionX();
    cout << "Calibrating bar " << x << ", " << y << ", " << z << endl;

    if (nRawPeaks >= 2) {
	double x1 = xRawPeaks[0];
	double x2 = xRawPeaks[1];
	
	double fitMin = x2-100;
	double fitMax = x2+200;

	TF1 *fRaw = new TF1("fRaw","landau",fitMin,fitMax);
	h1->Fit(fRaw,"R"); 
	
	double pedestal = x1;
	double muonPeak = fRaw->GetMaximumX();

	cout << "Raw peaks: " << pedestal << " , " << muonPeak << endl;

	// linear transformation params
	double a = 20.5/(muonPeak-pedestal);
	double b = -20.5*pedestal/(muonPeak-pedestal);
	 
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

    double calMin = xPeaks[0] - 6;
    double calMax = xPeaks[0] + 10;

    TF1 *fCal = new TF1("fCal","landau",calMin,calMax);
    h2->Fit(fCal,"R");

    double calibratedPeak = fCal->GetMaximumX();

    if (nCalPeaks > 0) {
	timer.Stop();
	double realTime = timer.RealTime();
	std::cout << "Time: " << realTime << "s. Calibrated peak: " << calibratedPeak << std::endl;
	csvFile << x << "," << y << "," << z << "," << calibratedPeak << "\n";
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
