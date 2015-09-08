#include "Discriminator.h"
#include "TGraph.h"
#include "TFile.h"

float Discriminator::sfVersion = 1.0;

bool Discriminator::sbCCMCutPass[4][NUM_BANDS] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
bool Discriminator::sbPGACutPass[4][NUM_BANDS] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
bool Discriminator::sbNGMCutPass[4][NUM_BANDS] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
bool Discriminator::sbWBSCutPass[4][NUM_BANDS] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
bool Discriminator::sbLAPCutPass[4][NUM_BANDS] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
bool Discriminator::sbDFTCutPass[4][NUM_BANDS] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0}};

bool Discriminator::sbCutPass[NUM_DISCRIMS][4][NUM_BANDS] = {{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
															{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
															{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
															{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
															{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},
															{{0,0,0},{0,0,0},{0,0,0},{0,0,0}}};

double Discriminator::dDiscrim[NUM_DISCRIMS][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

unique_ptr<TTree> Discriminator::tree = nullptr;

Discriminator::Discriminator(int channel): gain{0.0093634,0.0115473,0.0092957}, ch(channel), ciChan(4) {
	if (g_verbose) cout << "Discriminator " << ch << " c'tor\n";
	if ((ch >= ciChan) || (ch < 0)) {
		cout << error_message[method_error] << "Channel\n";
		return;
	}
	int d(0), b(0), i(0), chan(0);

	strcpy(cDiscrimNames[d_CCM_t], "CCM");
	strcpy(cDiscrimNames[d_DFT_t], "DFT");
	strcpy(cDiscrimNames[d_NGM_t], "NGM");
	strcpy(cDiscrimNames[d_LAP_t], "LAP");
	strcpy(cDiscrimNames[d_PGA_t], "PGA");
	strcpy(cDiscrimNames[d_WBS_t], "WBS");

	strcpy(cBandNames[mean], "mean");
	strcpy(cBandNames[sig1], "sig1");
	strcpy(cBandNames[sig3], "sig3");

	iDiscrimBins = 500;

	TGraph* gDiscrimCut = nullptr;

	unique_ptr<TFile> discrim_file;

	double* dDiscrimValue (nullptr);

	char cBand[16];

	try{discrim_file.reset(new TFile((sWorkingDir+"/Data_Processor/config/discrimination_bands.root").c_str(), "READ"));}
	catch (bad_alloc& ba) {cout << error_message[alloc_error] << "Bands\n"; return;}
	if (!discrim_file->IsOpen()){cout << error_message[file_error] << "Bands\n"; return;}
	for (d=0; d<NUM_DISCRIMS; d++){
		for (chan =0; chan<ciChan; chan++){
			for (b=0; b<NUM_BANDS; b++){
				try{
					vDiscrimBand[d][chan][b].reserve(iDiscrimBins);
				}
				catch (bad_alloc& ba) {
					iFailed = (1<<alloc_error);
					return;
				}
				sprintf(cBand, "%s_%d_%s", cDiscrimNames[d], chan ,cBandNames[b]);
				gDiscrimCut = (TGraph*)discrim_file->Get(cBand);
				if (gDiscrimCut==nullptr) {cout << error_message[root_error] << "Bands\n"; return;}
				dDiscrimValue = gDiscrimCut->GetY();
				for (i=0; i<iDiscrimBins; i++)vDiscrimBand[d][chan][b].push_back(dDiscrimValue[i]);
			}
		}
	}
	gDiscrimCut = nullptr;
	dDiscrimValue = nullptr;
	discrim_file->Close();
	discrim_file = nullptr;

}


Discriminator::~Discriminator(){
	if (g_verbose) cout << "Discriminator " << ch << " c'tor\n";
	for (int d=0; d<NUM_DISCRIMS; d++){
		for (int chan=0; chan<ciChan; chan++){
			for (int b=0; b<NUM_BANDS; b++){
//				vDiscrimBand[d][chan][b].reset();
			}
		}
	}
}

void Discriminator::CutsTree_init(TTree* tree_cuts){
	Discriminator::tree = unique_ptr<TTree>(tree_cuts);

	Discriminator::tree->Branch("CCM", Discriminator::sbCCMCutPass , "ccmPass[4][3]/O");
	Discriminator::tree->Branch("PGA", Discriminator::sbPGACutPass , "pgaPass[4][3]/O");
	Discriminator::tree->Branch("NGM", Discriminator::sbNGMCutPass , "ngmPass[4][3]/O");
	Discriminator::tree->Branch("WBS", Discriminator::sbWBSCutPass , "wbsPass[4][3]/O");
	Discriminator::tree->Branch("LAP", Discriminator::sbLAPCutPass , "lapPass[4][3]/O");
	Discriminator::tree->Branch("DFT", Discriminator::sbDFTCutPass , "dftPass[4][3]/O");
}

void Discriminator::SetDiscriminationValue(){
	dDiscrim[d_CCM_t][ch] = *slowint/(*fastint);
	dDiscrim[d_PGA_t][ch] = *sample/(*peakheight2);
	dDiscrim[d_NGM_t][ch] = *xsq_y/(*peakscale_y)-*xsq_n/(*peakscale_n);
	dDiscrim[d_WBS_t][ch] = *baseshift_n+*baseshift_y;
	dDiscrim[d_DFT_t][ch] = *dft_even/(*dft_odd);
	dDiscrim[d_LAP_t][ch] = *lap_high/(*lap_low);
}

void Discriminator::SetAddresses(vector<void*> add) {
	int i(0);
	fastint = (double*)add[i++];
	slowint = (double*)add[i++];
	sample = (double*)add[i++];
	peakheight2 = (double*)add[i++];
	xsq_n = (double*)add[i++];
	xsq_y = (double*)add[i++];
	peakscale_n = (double*)add[i++];
	peakscale_y = (double*)add[i++];
	baseshift_n = (double*)add[i++];
	baseshift_y = (double*)add[i++];
	dft_even = (double*)add[i++];
	dft_odd = (double*)add[i++];
	lap_high = (double*)add[i++];
	lap_low = (double*)add[i++];
	integral = (double*)add[i++];
}

void Discriminator::Discriminate(){

	for (int d=0; d<NUM_DISCRIMS; d++){
		for (int b=0; b<NUM_BANDS; b++) sbCutPass[d][ch][b] = false;
	}

	for (int d=0; d<NUM_DISCRIMS; d++){
		int iBinNumber = *integral/(gain[ch]*iDiscrimBins);
		for (int b=0; b<NUM_BANDS; b++){
			if (d_WBS_t||d_DFT_t){
				if (dDiscrim[d][ch]<vDiscrimBand[d][ch][b].at(iBinNumber)) sbCutPass[d][ch][b] = true;
			}
			else{
				if (dDiscrim[d][ch]>vDiscrimBand[d][ch][b].at(iBinNumber)) sbCutPass[d][ch][b] = true;
			}
		}
	}
}

void Discriminator::Cuts_fill(){
	for (int ch=0; ch<4; ch++){
		for (int b=0; b<NUM_BANDS; b++){
			sbCCMCutPass[ch][b] = sbCutPass[0][ch][b];
			sbDFTCutPass[ch][b] = sbCutPass[1][ch][b];
			sbNGMCutPass[ch][b] = sbCutPass[2][ch][b];
			sbLAPCutPass[ch][b] = sbCutPass[3][ch][b];
			sbPGACutPass[ch][b] = sbCutPass[4][ch][b];
			sbWBSCutPass[ch][b] = sbCutPass[5][ch][b];
		}
	}
	Discriminator::tree->Fill();
}

void Discriminator::FriendshipIsMagic() {
	tree->AddFriend("TS");
	tree->AddFriend("T0");
	tree->AddFriend("T1");
}