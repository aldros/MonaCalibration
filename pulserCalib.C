//Kevin Eisenberg, 2025
//Pulser calibration of MoNA, run5212 pulser run to extract ns/ch for each PMT

#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>

void pulserCalib(){
  gROOT->SetBatch();

	TFile *f = new TFile("/mnt/analysis/e23033/analysis/kevin/rootfiles/run5212.root");
	TTree *t = (TTree*)f->Get("t");

	ofstream csvFile("/mnt/analysis/e23033/analysis/kevin/outputs/pulserParams.csv");
    	csvFile << "x,y,z,slope\n";

	UShort_t time[2][18][16]; 
	t->SetBranchAddress("t[2][18][16]",time);

	Long64_t nEvents = t->GetEntries();

	const int nWalls = 9;
	const int nVert  = 16;
	const int nSides = 2;
	const int nPMTs = 288;

	vector<TH1F*> rawHist(nPMTs);

	for (int x = 0; x < nWalls; ++x) {
        	for (int y = 0; y < nVert; ++y) {
            		for (int z = 0; z < nSides; ++z) {
				int Index = x*nVert*nSides + y*nSides + z;
				cout << "Calibrating PMT " << x << "," << y << "," << z << endl;
				rawHist[Index] = new TH1F(Form("raw[%d][%d][%d]",x,y,z),"Pulser;TDC;Counts",512,0,4096);
				
			}
		}
	}

	for (int i=0; i<nEvents; i++) {
		t->GetEntry(i);
		for (int x=0; x<nWalls; x++) {
      			for (int y=0; y<nVert; y++) {
        			for (int z=0; z<nSides; z++) {
          				int Index = x*nVert*nSides + y*nSides + z;
          				rawHist[Index]->Fill(time[z][x][y]);
        			}
      			}
    		}
	}


	TSpectrum *spec = new TSpectrum(9);
	TF1 *fit = new TF1("fit", "pol1",0,4096);

	for (int i=0; i<nPMTs; i++){

		int x = i / (nVert*nSides);
	    	int y = (i / nSides) % nVert;
	    	int z = i % nSides;

		Int_t peaks = spec->Search(rawHist[i],1,"",0.2);
		Double_t *xPeaks = spec->GetPositionX();
		Int_t nPeaks = spec->GetNPeaks();
		//rawHist->Draw();
		vector<double> xArr(xPeaks, xPeaks+nPeaks);
		sort(xArr.begin(),xArr.end());
	
		vector<double> yArr;
		for (int i=0; i<nPeaks; i++){
			yArr.push_back(20.0*(i+1));
		}

	
		TGraph *g = new TGraph(xArr.size(), xArr.data(), yArr.data());
		g->Fit(fit,"Q");
	
		double slope = fit->GetParameter(1);

		cout << "slope (ns/ch): " << slope << endl;
		csvFile << x << "," << y << "," << z << "," << slope << "\n";
	}
	csvFile.close();

}
