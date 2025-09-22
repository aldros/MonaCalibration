// MoNA light calibration using cosmics
// Kevin Eisenberg, September 2025
// For use in FRIB e23033

#include <fstream>
#include "TStopwatch.h"


void kevinMonaQCalLOOP() {
    gROOT->SetBatch();

    TStopwatch timer; timer.Start();
		
    // random run
    TFile *f = new TFile("run-0034-00--raw.root");
    TTree *t = (TTree*)f->Get("t");

    // output file
    ofstream csvFile("calibratedPeaksFULL.csv");
    csvFile << "x,y,z,a,b,peak\n";

    UShort_t light[18][16][2];
    t->SetBranchAddress("light[18][16][2]", light);
    Long64_t nEntries = t->GetEntries(); 

    const int nWalls = 9; //9
    const int nVert = 16; //16
    const int nSides = 2; // L/R

    const int nPMTs = nWalls * nVert * nSides;

    vector<TH1F*> hRaw(nPMTs); //arrays of histograms
    vector<TH1F*> hCal(nPMTs);

    vector<double> a(nPMTs); //calibration constants
    vector<double> b(nPMTs);

    cout << "Beginning mass calibration..." << endl;

    for (int x=0; x<nWalls; x++){
	for (int y=0; y<nVert; y++){
	    for (int z=0; z<nSides; z++){
    		// declare hists
		int Index = x*nVert*nSides + y*nSides + z; // master index
		hRaw[Index] = new TH1F(Form("raw[%d][%d][%d]",x,y,z),"Raw charge;Channel;Counts",350,0,4096);
		hCal[Index] = new TH1F(Form("cal[%d][%d][%d]",x,y,z),"Calibrated charge;Energy (MeV);Counts",140,0,80);
	    }
	}
    }
    cout << "Histograms declared" << endl;
    // raw data FIRST LOOP
    for (Long64_t i=0; i<nEntries; i++){ 
	t->GetEntry(i);
	for (int x=0; x<nWalls; x++){
	     for (int y=0; y<nVert; y++){
		for (int z=0; z<nSides; z++){
			int Index = x*nVert*nSides + y*nSides + z; 
			hRaw[Index]->Fill(light[x][y][z]);
		}
	    }
	}
    }
    cout << "Raw histograms filled" << endl;

    TSpectrum *rawSpec = new TSpectrum(4);

    double x1;
    double x2;
    double fitMin;
    double fitMax;
    double pedestal;
    double muonPeak;

    for (int i=0; i<nPMTs; i++){
	
	TH1F* h1 = hRaw[i];

	int x = i / (nVert*nSides); //anti-conversion from master index
	int y = (i / nSides) % nVert;
	int z = i % nSides;

	Int_t nRawPeaks;
		
	    if (x==3 && y==4 && z==1){nRawPeaks = rawSpec->Search(h1,0.9,"",0.0001);} //bad pmt
	    else if (x==2 && y==5 && z==0){nRawPeaks = rawSpec->Search(h1,1,"",0.0001);} //bad pmt
	    else if (x==5 && y==15 && z==0){nRawPeaks = rawSpec->Search(h1,1,"",0.0001);} //bad pmt
	

	    else {
	    nRawPeaks = rawSpec->Search(h1,2,"",0.0001);
	    }

	    Double_t *xRawPeaks = rawSpec->GetPositionX();
	    cout << "Calibrating bar " << x << ", " << y << ", " << z << endl;

	    if (nRawPeaks >= 2) {
		x1 = xRawPeaks[0];
		x2 = xRawPeaks[1];
		
		if (x==5 && y==15 && z==0){fitMin = x2-30; fitMax = x2+70;}//bad pmt
		
		else {
		fitMin = x2-120;
		fitMax = x2+200;
		}

		TF1 *fRaw = new TF1("fRaw","landau",fitMin,fitMax);
		h1->Fit(fRaw,"R"); 
		
		pedestal = x1;
		muonPeak = fRaw->GetMaximumX();

		cout << "Raw peaks: " << pedestal << " , " << muonPeak << endl;

		// linear transformation params
		a[i] = 20.5/(muonPeak-pedestal);
		b[i] = -20.5*pedestal/(muonPeak-pedestal);
	    }
	    else if (nRawPeaks < 2) {cout << "1 peak: x,y,z = " << x << "," << y << "," << z << endl;}
    }

    // calibrated hist SECOND LOOP
    cout << "Beginning to fill calibrated histograms..." << endl;
    for (Long64_t i=0; i<nEntries; i++){
	t->GetEntry(i);
	for (int x=0; x<nWalls; x++){
	     for (int y=0; y<nVert; y++){
		for (int z=0; z<nSides; z++){
			int Index = x*nVert*nSides + y*nSides + z; 
			double calEntry = a[Index] * light[x][y][z] + b[Index];
			hCal[Index]->Fill(calEntry);
		}
    	     }
	}
    }
    cout << "Calibrated histograms filled. Writing to csv..." << endl;

    // output: location of muon peak in calibrated spectrum
    TSpectrum *calSpec = new TSpectrum(3);

    double calMin;
    double calMax;
    for (int i=0; i<nPMTs; i++){
	    TH1F* h2 = hCal[i];
	    Int_t nCalPeaks = calSpec->Search(h2,2,"",0.2);
	    Double_t *xPeaks = calSpec->GetPositionX();

	    calMin = xPeaks[0] - 6;
	    calMax = xPeaks[0] + 10;

	    TF1 *fCal = new TF1("fCal","landau",calMin,calMax);
	    h2->Fit(fCal,"R");

	    double calibratedPeak = fCal->GetMaximumX();

		
	    int x = i / (nVert*nSides);
	    int y = (i / nSides) % nVert;
	    int z = i % nSides;

	    if (nCalPeaks > 0) {
		std::cout << "Calibrated peak: " << calibratedPeak << "MeV" << std::endl;
		csvFile << x << "," << y << "," << z << "," << a[i] << "," << b[i] << "," << calibratedPeak << "\n";
    	    }
    }
    timer.Stop();
    double realTime = timer.RealTime();
    cout << "Completed in " << realTime << "s" << endl;
    csvFile.close();
}
