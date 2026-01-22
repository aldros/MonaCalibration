// Position calibration  for a single PMT
// Kevin Eisenberg, September 2025
// For use in FRIB e23033

void singleTimeCalib(int x, int y){

	TFile *f = new TFile("run5200.root");
    	TTree *t = (TTree*)f->Get("t");

	UShort_t time[2][18][16]; 
	t->SetBranchAddress("time[2][18][16]", time); 
	UShort_t light[2][18][16];
	t->SetBranchAddress("light[2][18][16]", light);
	
	//declare hists
	TH1F *hRaw = new TH1F("hRaw","raw time diff",200,-400,400);
	TH1F *hCal = new TH1F("hCal","calibrated position;x (cm);counts",100,-200,200); 

	//raw data (don't need all entries)
	Long64_t nEntries = t->GetEntries();
	double dt;
	for (int i=0; i<1000000; i++){
	t->GetEntry(i); 
		if (time[0][x][y]>0 && time[1][x][y]>0){
			dt = static_cast<double>(time[0][x][y]) - static_cast<double>(time[1][x][y]);
			hRaw->Fill( dt );
		}
	}
	
	TF1 *fRaw = new TF1("fRaw","[p0]*(tanh((x-[p1])/[p2])-tanh((x-[p3])/[p4]))");
	fRaw->SetParameters(hRaw->GetMaximum(),-150,10,150,10);
	
	double fitmin = -400.0;
	double fitmax = 400.0;

	hRaw->Fit("fRaw","R","",fitmin,fitmax);

	double amplitude,xL,tauL,xR,tauR;
	amplitude = fRaw->GetParameter(0);
	xL = fRaw->GetParameter(1);
	tauL = fRaw->GetParameter(2);
	xR = fRaw->GetParameter(3);
	tauR = fRaw->GetParameter(4);

	//optional logarithmic correction 
	//double alpha = 0.01;
	//xL += tauL*log(alpha/(1-alpha));
	//xR -= tauR*log(alpha/(1-alpha));

	//fit parameters
	double a,b;
	
	a = 200/(xR-xL);
	b = -100 - 200*xL/(xR-xL);
	
	//calibrated histogram
	double calEntry;
	for (int i=0; i<1000000; i++){
	t->GetEntry(i); 
		if (time[0][x][y]>0 && time[1][x][y]>0){
			calEntry = b + a*(time[0][x][y]-time[1][x][y]);
			hCal->Fill( calEntry );
		}
	}	
	cout << "a: " << a << endl;
	cout << "b: " << b << endl;
	TCanvas *c = new TCanvas("c","MoNA time calibration",800,600);
	c->Divide(2,1);

	c->cd(1);
	hRaw->Draw(); 
	//gPad->SetLogy();
	c->cd(2);
	hCal->Draw();
	//gPad->SetLogy();

}
