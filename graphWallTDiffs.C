//Kevin Eisenberg, 2025
//This file takes a wall of MoNA from a (paddle) run and calculates offsets in ns for relative alignment of each bar to the top bar. Each bar should be ~0.34 ns later than the one above it, so we offset the calculated time difference to be multiples of 0.34ns from the top bar.

#include <fstream>

void graphWallTDiffs(int wallNumber, string runNumber){

	//GRAPH A LARGE GRID OF TIME DIFFERENCE GRAPHS WRT THE TOP BAR
	TTree *params = new TTree("params","params");
	params->ReadFile("/mnt/analysis/e23033/analysis/kevin/outputs/pulserParams.csv","x/I:y/I:z/I:slope/D"); // PARAMETERS TO CONVERT FROM CHANNEL TO ns

	Long64_t nParams = params->GetEntries();

	int x; int y; int z; double slope;
	params->SetBranchAddress("x",&x);
	params->SetBranchAddress("y",&y);
	params->SetBranchAddress("z",&z);
	params->SetBranchAddress("slope",&slope);

	TFile *f = new TFile(Form("/mnt/analysis/e23033/analysis/kevin/rootfiles/run%s.root",runNumber.c_str()));
	TTree *t = (TTree*)f->Get("t");

	UShort_t time[2][18][16];
	t->SetBranchAddress("time[2][18][16]",time);

	Long64_t nEvents = t->GetEntries();

	const int nWalls = 9;
	const int nVert  = 16;
	const int nSides = 2;
	const int nPMTs  = nWalls * nVert * nSides;

	double topLeftSlope;
	double topRightSlope;

	//MAKE DICTIONARY OF SLOPES
	vector<double> slopes(nPMTs);
	for (int i=0; i<nParams; i++) {
		params->GetEntry(i);
		int Index = x*nVert*nSides + y*nSides + z;
		slopes[Index] = slope;
		if (x==wallNumber && y==15 && z==0){
			topLeftSlope = slope;
		}
		else if (x==wallNumber && y==15 && z==1){
			topRightSlope = slope;
		}

	}
	//DECLARE HISTOGRAMS (since each bar is compared to the top bar, there's 15 histograms)
	vector<TH1F*> hWall(nVert-1);
	for (int j=0; j<(nVert-1); j++){
		hWall[j] = new TH1F(Form("h%d",j),Form("Time Diff %d to 15;ns;counts",j),300,-50,50);
	}
	double topAVG;
	double barAVG;
	double timeDiff;
	for (int i=0; i<nEvents; i++){
		t->GetEntry(i);
		if (time[0][wallNumber][15]==0 || time[1][wallNumber][15]==0)
			continue;
		for (int y=0; y<nVert-1; y++){
			if (time[0][wallNumber][0]>0 && time[1][wallNumber][0]>0 && time[0][wallNumber][y]>0 && time[1][wallNumber][y]>0){ //ensure all times nonzero, and also ensure entry is collinear, vertical through ALL bars
				int leftIndex = wallNumber*nVert*nSides + y*nSides;
				int rightIndex = wallNumber*nVert*nSides + y*nSides + 1;
				double leftSlope = slopes[leftIndex];
				double rightSlope = slopes[rightIndex];
				barAVG = ((time[0][wallNumber][y]*leftSlope) + (time[1][wallNumber][y]*rightSlope))/2.0;

				topAVG = ((time[0][wallNumber][15]*topLeftSlope) + (time[1][wallNumber][15]*topRightSlope))/2.0;
				timeDiff = barAVG - topAVG;

				hWall[y]->Fill(timeDiff); //quantity in question
			}	
			
		}
	}
	// GRAPH THEM ALL?
	TCanvas *c = new TCanvas("c","Wall time diffs");
	c->Divide(5,3);

	TF1 *gaus = new TF1("gaus","gaus");
	double rawMax;
	int binMax;
	
	for (int j = 0; j < nVert - 1; j++) {
        	c->cd(j + 1);
		binMax = hWall[j]->GetMaximumBin();
		rawMax = hWall[j]->GetBinCenter(binMax);
		hWall[j]->Fit("gaus","Q","",rawMax-2,rawMax+2);
		gaus->SetLineWidth(1);
	
		double centroid = gaus->GetParameter(1);
		//cout << j << " centroid: " << centroid << endl;
		cout << wallNumber << "," << j << "," << (16-(j+1))*0.34 - centroid << endl;
		hWall[j]->Draw();
		TLine *target = new TLine((16-(j+1))*0.34,0,(16-(j+1))*0.34,180);
		target->SetLineColor(kGreen);
		target->Draw();

	}

	c->Update();

}
