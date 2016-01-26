#include "Discriminator.h"
#include "TGraph.h"

float Discriminator::sfVersion = 1.07;
// 0.0093634,0.0115473,0.0092957 before 160125
Discriminator::Discriminator() : gain{0.0098938,0.011679,0.010655} {

	strcpy(cDiscrimNames[d_CCM_t], "CCM");
	strcpy(cDiscrimNames[d_DFT_t], "DFT");
	strcpy(cDiscrimNames[d_NGM_t], "NGM");
	strcpy(cDiscrimNames[d_LAP_t], "LAP");
	strcpy(cDiscrimNames[d_PGA_t], "PGA");
	strcpy(cDiscrimNames[d_WBS_t], "WBS");

	strcpy(cBandNames[mean], "mean");
	strcpy(cBandNames[sig1], "sig1");
	strcpy(cBandNames[sig3], "sig3");

	memset(sbCutPass, 0, sizeof(sbCutPass));
	memset(dDiscrim, 0, sizeof(dDiscrim));
	memset(dDiscrimBand, 0, sizeof(dDiscrimBand));

	TGraph* gDiscrimCut = nullptr;

	unique_ptr<TFile> discrim_file;

	double* dDiscrimValue (nullptr);

	char cBand[16];

	try{discrim_file.reset(new TFile((sConfigDir+"/config/discrimination_bands.root").c_str(), "READ"));}
	catch (bad_alloc& ba) {
		cout << error_message[alloc_error] << "Bands\n";
		iFailed = 1;
		return;
	}
	if (!discrim_file->IsOpen()) {
		cout << error_message[file_error] << "Bands\n";
		iFailed = 1;
		return;
	}
	for (int d = 0; d < NUM_DISCRIMS; d++) {
		for (int ch = 0; ch < NUM_CHANS; ch++) {
			for (int b = 0; b < NUM_BANDS; b++) {
				sprintf(cBand, "%s_%d_%s", cDiscrimNames[d], ch ,cBandNames[b]);
				gDiscrimCut = (TGraph*)discrim_file->Get(cBand);
				if (gDiscrimCut == nullptr) {
					cout << error_message[root_error] << "Band " << cBand << "\n";
					iFailed = 1;
					return;
				}
				dDiscrimValue = gDiscrimCut->GetY();
				for (int i = 0; i < iDiscrimBins; i++) dDiscrimBand[d][ch][b][i] = dDiscrimValue[i];
			}
		}
	}
	gDiscrimCut = nullptr;
	dDiscrimValue = nullptr;
	discrim_file->Close();
	discrim_file = nullptr;
}

Discriminator::~Discriminator(){
	if (g_verbose > 1) cout << "Discriminator c'tor\n";
	pPeak2 = nullptr;
	T0.reset();
	T2.reset();
	if (f) f->Close();
	f.reset();
}

void Discriminator::Setup(string filein) {
	unsigned short mask;
	bNewPeak2 = false;
	pPeak2 = &vPeak2;
	f.reset(new TFile((sWorkingDir + "/prodata/" + filein + ".root").c_str(), "UPDATE"));
	if (f->IsZombie()) {
		cout << error_message[root_error] << "TFile\n";
		iFailed = 1;
		return;
	} else {
		f->cd();
	}
	iNChan = 0;
	memset(iChan, 0, sizeof(iChan));
	T0 = unique_ptr<TTree>((TTree*)f->Get("TI"));
	if (!T0) {
		cout << error_message[root_error] << "TI\n";
		iFailed = 1;
		return;
	} else {
		T0->SetBranchAddress("Mask", &mask);
		T0->GetEntry(0);
		for (int i = 0; i < NUM_CHANS; i++) if (mask & (1 << i)) iChan[iNChan++] = i;
		T0.reset();
	}

	T0 = unique_ptr<TTree>((TTree*)f->Get("T0"));
	if (!T0) {
		cout << error_message[root_error] << "T0\n";
		iFailed = 1;
		return;
	} else {
		T0->SetBranchStatus("*",0);
		T0->SetBranchStatus("Integral",			1);
		T0->SetBranchStatus("FastInt",			1);
		T0->SetBranchStatus("SlowInt",			1);
		T0->SetBranchStatus("Sample",			1);
		T0->SetBranchStatus("Peakheight2",		1);
		T0->SetBranchStatus("Xsq_n_f",			1);
		T0->SetBranchStatus("Xsq_y_f",			1);
		T0->SetBranchStatus("Peakscale_n_f",	1);
		T0->SetBranchStatus("Peakscale_y_f",	1);
		T0->SetBranchStatus("Base_shift_n",		1);
		T0->SetBranchStatus("Base_shift_y",		1);
		T0->SetBranchStatus("Even",				1);
		T0->SetBranchStatus("Odd",				1);
		T0->SetBranchStatus("LapHi",			1);
		T0->SetBranchStatus("LapLow",			1);

		T0->SetBranchAddress("Integral",		integral);
		T0->SetBranchAddress("FastInt",			fastint);
		T0->SetBranchAddress("SlowInt",			slowint);
		T0->SetBranchAddress("Sample",			sample);
		T0->SetBranchAddress("Xsq_n_f",			xsq_n);
		T0->SetBranchAddress("Xsq_y_f",			xsq_y);
		T0->SetBranchAddress("Peakscale_n_f",	peakscale_n);
		T0->SetBranchAddress("Peakscale_y_f",	peakscale_y);
		T0->SetBranchAddress("Base_shift_n",	baseshift_n);
		T0->SetBranchAddress("Base_shift_y",	baseshift_y);
		T0->SetBranchAddress("Even",			dft_even);
		T0->SetBranchAddress("Odd",				dft_odd);
		T0->SetBranchAddress("LapHi",			lap_high);
		T0->SetBranchAddress("LapLow",			lap_low);

		if (T0->SetBranchAddress("Peakheight2",		dPeak2)) {
			bNewPeak2 = true;
			T0->SetBranchAddress("Peakheight2",		&pPeak2);
		}

		lNumEvents = T0->GetEntries();
	}

	try {T2 = unique_ptr<TTree>(new TTree("T2","Distriminator"));}
	catch (bad_alloc& ba) {
		cout << error_message[alloc_error] << "T2\n";
		iFailed = 1;
		return;
	}
	if (T2->IsZombie()) {
		cout << error_message[root_error] << "T2\n";
		iFailed = 1;
		return;
	} else {
		T2->Branch("CCM", sbCutPass[d_CCM_t], "ccmPass[4][3]/O");
		T2->Branch("PGA", sbCutPass[d_PGA_t], "pgaPass[4][3]/O");
		T2->Branch("NGM", sbCutPass[d_NGM_t], "ngmPass[4][3]/O");
		T2->Branch("WBS", sbCutPass[d_WBS_t], "wbsPass[4][3]/O");
		T2->Branch("LAP", sbCutPass[d_LAP_t], "lapPass[4][3]/O");
		T2->Branch("DFT", sbCutPass[d_DFT_t], "dftPass[4][3]/O");
	}
}

void Discriminator::Discriminate() {
	if (g_verbose) cout << "Discriminating...\n";
	int iBinNumber(0);
	for (long e = 0; e < lNumEvents; e++) {
		T0->GetEntry(e);
		for (int ch = 0; ch < iNChan; ch++) {
			dDiscrim[d_CCM_t][iChan[ch]] = fastint[iChan[ch]] == 0 ?									-1 : slowint[iChan[ch]]/fastint[iChan[ch]];
			dDiscrim[d_DFT_t][iChan[ch]] = dft_odd[iChan[ch]] == 0 ?									-1 : dft_even[iChan[ch]]/dft_odd[iChan[ch]];
			dDiscrim[d_NGM_t][iChan[ch]] = peakscale_y[iChan[ch]] == 0 || peakscale_n[iChan[ch]] == 0 ?	-1 : xsq_y[iChan[ch]]/peakscale_y[iChan[ch]]-xsq_n[iChan[ch]]/peakscale_n[iChan[ch]];
			dDiscrim[d_LAP_t][iChan[ch]] = lap_low[iChan[ch]] == 0 ?									-1 : lap_high[iChan[ch]]/lap_low[iChan[ch]];
			if (bNewPeak2) dDiscrim[d_PGA_t][iChan[ch]] = pPeak2->at(iChan[ch])[0] == 0 ?				-1 : sample[iChan[ch]]/pPeak2->at(iChan[ch])[0];
			else dDiscrim[d_PGA_t][iChan[ch]] = dPeak2[iChan[ch]] == 0 ?								-1 : sample[iChan[ch]]/dPeak2[iChan[ch]];
			dDiscrim[d_WBS_t][iChan[ch]] = baseshift_n[iChan[ch]]+baseshift_y[iChan[ch]];
			iBinNumber = integral[iChan[ch]]/(gain[ch]*iDiscrimBins); // gain isn't handled the same way as the stuff from the tree
			for (int d = 0; d < NUM_DISCRIMS; d++) {
				for (int b = 0; b < NUM_BANDS; b++) {
					sbCutPass[d][iChan[ch]][b] = false;
					if ((d == d_WBS_t) || (d == d_DFT_t)) {
						if (dDiscrim[d][iChan[ch]] < dDiscrimBand[d][iChan[ch]][b][iBinNumber]) sbCutPass[d][iChan[ch]][b] = true;
					} else if (dDiscrim[d][iChan[ch]] > dDiscrimBand[d][iChan[ch]][b][iBinNumber]) sbCutPass[d][iChan[ch]][b] = true;
				} // b
			} // d
		} // ch
		T2->Fill();
	} // e
	if (g_verbose) cout << "Discriminated\nMaking friends";
	T2->AddFriend("TS");
	T2->AddFriend("T0");
	T2->AddFriend("T1");
	T0->AddFriend("T2");
	T0->Write("T0", TObject::kOverwrite);
	T2->Write("T2", TObject::kOverwrite);
	T0.reset((TTree*)f->Get("T1"));
	if (T0) {
		T0->AddFriend("T2");
		T0->Write("T1", TObject::kOverwrite);
	}
	T0.reset();
	T2.reset();
	f->Close();
	f.reset();
}
