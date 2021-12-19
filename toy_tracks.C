#include <cstdio>
#define toy_tracks_cxx

#include "hits.h"

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <vector>

#include <algorithm>
#include <fstream>

using std::vector;

struct pseudo_track
{
	int straw_id;
	int scint_id;
	int gem_id;
};

std::string prd(const double x, const int decDigits, const int width) {
    stringstream ss;
    ss << fixed << right;
    ss.fill(' ');        // fill space around displayed #
    ss.width(width);     // set  width around displayed #
    ss.precision(decDigits); // set # places after decimal
    ss << x;
    return ss.str();
};

std::string center(const string s, const int w) {
    stringstream ss, spaces;
    int padding = w - s.size();                 // count excess room to pad
    for(int i=0; i<padding/2; ++i)
        spaces << " ";
    ss << spaces.str() << s << spaces.str();    // format with padding
    if(padding>0 && padding%2!=0)               // if odd #, add 1 space
        ss << " ";
    return ss.str();
};

// bool compareByTime(const pseudo_track &a, const pseudo_track &b) { return a.time < b.time; }

// bool comp(const pseudo_track &lhs, const pseudo_track &rhs) { return lhs.id == rhs.id; }

void hits::Loop()
{
	TH2D *GemY_vs_StrawTypeA = new TH2D(
		"GemY_vs_StrawTypeA", "GemY_vs_StrawTypeA; StrawTypeA, ch; GEM3 Y-plane, ch",
		32, 0, 64, 128, 0, 256);

	TH2D *GemY_vs_StrawTypeB = new TH2D(
		"GemY_vs_StrawTypeB", "GemY_vs_StrawTypeB; StrawTypeB, ch; GEM3 Y-plane, ch",
		32, 0, 64, 128, 0, 256);

	// double tStart = 18;
	// double tEnd = 38;

	double sigma_gem_straw = 21.69;
	double mean_gem_straw = 28.75;

	double sigma_gem_scint = 16.95;
	double mean_gem_scint = -136.5;

	vector<pseudo_track> typeA;
	vector<pseudo_track> typeB;

	if (fChain == 0)
		return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry = 0; jentry < nentries; jentry++)
	{
		Long64_t ientry = LoadTree(jentry);

		if (ientry < 0)
			break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		for (int i = 0; i < hits_; i++)
		{
			if (!hits_over_threshold[i])
				continue;

			if (hits_det[i] == 3 && hits_plane[i] == 1)
			{
				long double gemTime = (long double)hits_time[i];
				int scinId = 0;

				for (int j = 0; j < hits_; j++)
				{
					if (i == j)
						continue;

					if (abs(gemTime - (long double)hits_time[j]) > 300)
						continue;

					if (hits_fec[j] == 2 && hits_vmm[j] == 9 && hits_ch[j] == 0)
					{
						if ((long double)hits_time[j] - gemTime > mean_gem_scint + 4 * sigma_gem_scint ||
							(long double)hits_time[j] - gemTime < mean_gem_scint - 4 * sigma_gem_scint)
						{
							continue;
						}
						else
						{
							scinId = (int)hits_id[j];
						}
					}

					if (hits_fec[j] == 2 && hits_vmm[j] == 10)
					{
						if (gemTime - (long double)hits_time[j] > mean_gem_straw + 4 * sigma_gem_straw ||
							gemTime - (long double)hits_time[j] < mean_gem_straw - 4 * sigma_gem_straw)
						{
							continue;
						}
						else
						{
							long double trawTime = (long double)hits_time[j];
							if (hits_vmm[i] > 1)
							{
								hits_vmm[i] -= 2;
							}
							else
							{
								hits_vmm[i] += 2;
							}
							if (hits_id[j] % 4 == 0)
							{
								GemY_vs_StrawTypeA->Fill((int)hits_ch[j], 64 * (int)hits_vmm[i] + (int)hits_ch[i]);
								if (scinId == 0)
								{
									typeA.push_back({(int)hits_id[j], -1 ,(int)hits_id[i]});
								}
								else
								{
									typeA.push_back({(int)hits_id[j], scinId ,(int)hits_id[i]});
								}
							}
							else if (hits_id[j] % 4 == 3)
							{
								GemY_vs_StrawTypeB->Fill((int)hits_ch[j], 64 * (int)hits_vmm[i] + (int)hits_ch[i]);
								if (scinId == 0)
								{
									typeB.push_back({(int)hits_id[j], -1 ,(int)hits_id[i]});
								}
								else
								{
									typeB.push_back({(int)hits_id[j], scinId ,(int)hits_id[i]});
								}
							}
							else
							{
								continue;
							}
						}
					}
				}
			}
		}
	}

	ofstream file_TypeA;
	file_TypeA.open ("TypeA.txt");


	file_TypeA << center("STRAW ID",15) << " | "
          	  << center("SCINT ID",15) << " | "
         	  << center("GEM3Y ID",15) << "\n";
	file_TypeA << std::string(15*3 + 2*3, '-') << "\n";
	
	for (int i = 0; i < typeA.size(); i++)
	{
		file_TypeA << prd(typeA[i].straw_id,0,15) << " | "
          	      << prd(typeA[i].scint_id,0,15) << " | "
         	      << prd(typeA[i].gem_id,0,15)   << "\n";
	}

	ofstream file_TypeB;
	file_TypeB.open ("TypeB.txt");


	file_TypeB << center("STRAW ID",15) << " | "
          	  << center("SCINT ID",15) << " | "
         	  << center("GEM3Y ID",15) << "\n";
	file_TypeB << std::string(15*3 + 2*3, '-') << "\n";
	
	for (int i = 0; i < typeB.size(); i++)
	{
		file_TypeB << prd(typeB[i].straw_id,0,15) << " | "
          	      << prd(typeB[i].scint_id,0,15) << " | "
         	      << prd(typeB[i].gem_id,0,15)   << "\n";
	}
	

	gStyle->SetOptFit();

	TFile *out = new TFile("toy_track_1522_full_4sigma.root", "RECREATE");
	GemY_vs_StrawTypeA->Write("GemY_vs_StrawTypeA");
	GemY_vs_StrawTypeB->Write("GemY_vs_StrawTypeB");
	out->Close();
}
