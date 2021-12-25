#include <cstdio>
#define toy_tracks_cxx

#include "hits.h"
#include "tqdm/tqdm.h"

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
	int gem_ch;
	double gem_weight;
	long double gem_scint;
	long double gem_straw;
};

std::string prd(const double x, const int decDigits, const int width)
{
	stringstream ss;
	ss << fixed << right;
	ss.fill(' ');			 // fill space around displayed #
	ss.width(width);		 // set  width around displayed #
	ss.precision(decDigits); // set # places after decimal
	ss << x;
	return ss.str();
};

std::string center(const string s, const int w)
{
	stringstream ss, spaces;
	int padding = w - s.size(); // count excess room to pad
	for (int i = 0; i < padding / 2; ++i)
		spaces << " ";
	ss << spaces.str() << s << spaces.str(); // format with padding
	if (padding > 0 && padding % 2 != 0)	 // if odd #, add 1 space
		ss << " ";
	return ss.str();
};

void hits::Loop()
{
	TH2D *GemY_vs_StrawTypeA = new TH2D(
		"GemY_vs_StrawTypeA", "GemY_vs_StrawTypeA; StrawTypeA, ch; GEM3 Y-plane, ch",
		32, 0, 64, 128, 0, 256);

	TH2D *GemY_vs_StrawTypeB = new TH2D(
		"GemY_vs_StrawTypeB", "GemY_vs_StrawTypeB; StrawTypeB, ch; GEM3 Y-plane, ch",
		32, 0, 64, 128, 0, 256);

	double sigma_gem_straw = 21.69;
	double mean_gem_straw = 28.75;

	double sigma_gem_scint = 16.95;
	double mean_gem_scint = -136.5;

	vector<pseudo_track> typeA;
	vector<pseudo_track> typeB;

	vector<vector< pseudo_track >> typeAbasedTracks;

	int gem_ch = 0;
	double gem_weight = 0;

	long double gem_scint_deltaT = 0;
	long double gem_straw_deltaT = 0;

	if (fChain == 0)
		return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	// for (Long64_t jentry : tqdm::range(nentries))
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
			{
				// if (hits_det[i] == 3 && hits_plane[i] == 1)
				// {
				// 	long double gemTime = (long double)hits_time[i];
				// 	for (int j = 0; j < hits_; j++)
				// 	{
				// 		if (!hits_over_threshold[j])
				// 			continue;
				// 		if (i == j)
				// 			continue;

				// 		if (abs(gemTime - (long double)hits_time[j]) < 300)
				// 			std::cout << "Det det " << (int)hits_det[j] << " det vmm (int)hits_vmm[j] " << (int)hits_vmm[j] << "\n";
				// 		else
				// 		{
				// 			continue;
				// 		}
				// 	}
				// }
				
				continue;
			}
			// else 
			// {
			// 	continue;
			// }
			if (hits_det[i] == 3 && hits_plane[i] == 1)
			{
				long double gemTime = (long double)hits_time[i];

				if (hits_vmm[i] > 1)
				{
					hits_vmm[i] -= 2;
				}
				else
				{
					hits_vmm[i] += 2;
				}
				gem_ch = 64 * (int)hits_vmm[i] + (int)hits_ch[i];
				gem_weight = (int)hits_adc[i] / 1024.0;
				int scinId = -1;

				for (int j = 0; j < hits_; j++)
				{
					if (!hits_over_threshold[j])
						continue;
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
							gem_scint_deltaT = gemTime - (long double)hits_time[j];
							scinId = (int)hits_id[j];
						}
					}

					if (scinId == -1)
						gem_scint_deltaT = -1.0;

					if (hits_fec[j] == 2 && hits_vmm[j] == 10)
					{
						if (gemTime - (long double)hits_time[j] > mean_gem_straw + 4 * sigma_gem_straw ||
							gemTime - (long double)hits_time[j] < mean_gem_straw - 4 * sigma_gem_straw)
						{
							continue;
						}
						else
						{
							gem_straw_deltaT = gemTime - (long double)hits_time[j];
							if (hits_id[j] % 4 == 0)
							{
								GemY_vs_StrawTypeA->Fill((int)hits_ch[j], gem_ch);
								typeA.push_back({(int)hits_id[j], scinId, (int)hits_id[i], gem_ch, gem_weight, gem_scint_deltaT, gem_straw_deltaT});
							}
							else if (hits_id[j] % 4 == 3)
							{
								GemY_vs_StrawTypeB->Fill((int)hits_ch[j], gem_ch);
								typeB.push_back({(int)hits_id[j], scinId, (int)hits_id[i], gem_ch, gem_weight, gem_scint_deltaT, gem_straw_deltaT});
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

	int typeA_all = 0, typeA_wScint = 0, typeA_woScint = 0;
	int typeB_all = 0, typeB_wScint = 0, typeB_woScint = 0;

	ofstream file_TypeA;
	file_TypeA.open("TypeA_w_utt.txt");

	file_TypeA << center("STRAW ID", 15) << " | "
			   << center("SCINT ID", 15) << " | "
			   << center("GEM3Y ID", 15) << " | "
			   << center("GEM3Y CH", 15) << " | "
			   << center("GEM - SCINT dT, ns", 20) << " | "
			   << center("GEM - STRAW dT, ns", 20) << "\n";
	file_TypeA << std::string(15 * 6 + 3 * 6 + 10, '-') << "\n";

	for (int i = 1; i < typeA.size(); i++)
	{
		if (typeA[i].scint_id != -1)
			typeA_wScint++;
		else
			typeA_woScint++;
		typeA_all++;

		if (typeA[i].straw_id == typeA[i - 1].straw_id && typeA[i].scint_id == typeA[i - 1].scint_id)
		{
			int j = i;
			int counter = 0;
			file_TypeA << prd(typeA[i - 1].straw_id, 0, 15) << " | "
					   << prd(typeA[i - 1].scint_id, 0, 15) << " | "
					   << prd(typeA[i - 1].gem_id, 0, 15) << " | "
					   << prd(typeA[i - 1].gem_ch, 0, 15) << " | "
					   << prd(typeA[i - 1].gem_scint, 0, 20) << " | "
					   << prd(typeA[i - 1].gem_straw, 0, 20) << "\n";
			while (typeA[j].straw_id == typeA[j - 1].straw_id && typeA[j].scint_id == typeA[j - 1].scint_id)
			{
				file_TypeA << std::string(15, ' ') << " | "
						   << std::string(15, ' ') << " | "
						   << prd(typeA[j].gem_id, 0, 15) << " | "
						   << prd(typeA[j].gem_ch, 0, 15) << " | "
						   << prd(typeA[j].gem_scint, 0, 20) << " | "
						   << prd(typeA[j].gem_straw, 0, 20) << "\n";
				j++;
				counter++;
			}
			i += counter - 1;
			file_TypeA << std::string(15 * 6 + 3 * 6 + 10, '-') << "\n";
		}
		else if (typeA[i].straw_id == typeA[i + 1].straw_id && typeA[i].scint_id == typeA[i + 1].scint_id)
		{
			continue;
		}
		else
		{
			file_TypeA << prd(typeA[i].straw_id, 0, 15) << " | "
					   << prd(typeA[i].scint_id, 0, 15) << " | "
					   << prd(typeA[i].gem_id, 0, 15) << " | "
					   << prd(typeA[i].gem_ch, 0, 15) << " | "
					   << prd(typeA[i].gem_scint, 0, 20) << " | "
					   << prd(typeA[i].gem_straw, 0, 20) << "\n";
			file_TypeA << std::string(15 * 6 + 3 * 6 + 10, '-') << "\n";
		}
	}

	printf("%d Pseudo tracks for TypeA straws: %d with SCINT and %d without \n", typeA_all, typeA_wScint, typeA_woScint);

	ofstream file_TypeB;
	file_TypeB.open("TypeB_w_utt.txt");

	file_TypeB << center("STRAW ID", 15) << " | "
			   << center("SCINT ID", 15) << " | "
			   << center("GEM3Y ID", 15) << " | "
			   << center("GEM3Y CH", 15) << " | "
			   << center("GEM - SCINT dT, ns", 20) << " | "
			   << center("GEM - STRAW dT, ns", 20) << "\n";
	file_TypeB << std::string(15 * 6 + 3 * 6 + 10, '-') << "\n";

	for (int i = 1; i < typeB.size(); i++)
	{
		if (typeB[i].scint_id != -1)
			typeB_wScint++;
		else
			typeB_woScint++;
		typeB_all++;

		if (typeB[i].straw_id == typeB[i - 1].straw_id && typeB[i].scint_id == typeB[i - 1].scint_id)
		{
			int j = i;
			int counter = 0;
			file_TypeB << prd(typeB[i - 1].straw_id, 0, 15) << " | "
					   << prd(typeB[i - 1].scint_id, 0, 15) << " | "
					   << prd(typeB[i - 1].gem_id, 0, 15) << " | "
					   << prd(typeB[i - 1].gem_ch, 0, 15) << " | "
					   << prd(typeB[i - 1].gem_scint, 0, 20) << " | "
					   << prd(typeB[i - 1].gem_straw, 0, 20) << "\n";
			while (typeB[j].straw_id == typeB[j - 1].straw_id && typeB[j].scint_id == typeB[j - 1].scint_id)
			{
				file_TypeB << std::string(15, ' ') << " | "
						   << std::string(15, ' ') << " | "
						   << prd(typeB[j].gem_id, 0, 15) << " | "
						   << prd(typeB[j].gem_ch, 0, 15) << " | "
						   << prd(typeB[j].gem_scint, 0, 20) << " | "
						   << prd(typeB[j].gem_straw, 0, 20) << "\n";
				j++;
				counter++;
			}
			i += counter - 1;
			file_TypeB << std::string(15 * 6 + 3 * 6 + 10, '-') << "\n";
		}
		else if (typeB[i].straw_id == typeB[i + 1].straw_id && typeB[i].scint_id == typeB[i + 1].scint_id)
		{
			continue;
		}
		else
		{
			file_TypeB << prd(typeB[i].straw_id, 0, 15) << " | "
					   << prd(typeB[i].scint_id, 0, 15) << " | "
					   << prd(typeB[i].gem_id, 0, 15) << " | "
					   << prd(typeB[i].gem_ch, 0, 15) << " | "
					   << prd(typeB[i].gem_scint, 0, 20) << " | "
					   << prd(typeB[i].gem_straw, 0, 20) << "\n";
			file_TypeB << std::string(15 * 6 + 3 * 6 + 10, '-') << "\n";
		}
	}

	printf("%d Pseudo tracks for typeB straws: %d with SCINT and %d without \n", typeB_all, typeB_wScint, typeB_woScint);

	for (int i = 1; i < typeA.size(); i++)
	{
		if (typeA[i].straw_id == typeA[i - 1].straw_id && typeA[i].scint_id == typeA[i - 1].scint_id)
		{
			int j = i;
			int counter = 0;
			vector< pseudo_track > tempTrack;
			tempTrack.push_back(typeA[j-1]);
			while (typeA[j].straw_id == typeA[j - 1].straw_id && typeA[j].scint_id == typeA[j - 1].scint_id)
			{
				tempTrack.push_back(typeA[j]);
				j++;
				counter++;
			}
			typeAbasedTracks.push_back(tempTrack);
			i += counter - 1;
		}
		else if (typeA[i].straw_id == typeA[i + 1].straw_id && typeA[i].scint_id == typeA[i + 1].scint_id)
		{
			continue;
		}
		else
		{
			vector< pseudo_track > tempTrack;
			tempTrack.push_back(typeA[i]);
			typeAbasedTracks.push_back(tempTrack);
		}
	}

	for (int i = 0; i < typeAbasedTracks.size(); i++) 
	{
		double sum = 0;
		double w_sum = 0;
		for (int j = 0; j < typeAbasedTracks[i].size(); j++)
		{
			sum+= typeAbasedTracks[i][j].gem_ch * typeAbasedTracks[i][j].gem_weight;
			w_sum += typeAbasedTracks[i][j].gem_weight;
		}
		double w_mean = (int)(sum / w_sum * 100.0) / 100.0 ;

		printf("Straw id \t %d \t and w_mean \t %f \n", (int)typeAbasedTracks[i][0].straw_id, w_mean);
		
	}



	gStyle->SetOptFit();

	TFile *out = new TFile("toy_track_1522_full_4sigma.root", "RECREATE");
	GemY_vs_StrawTypeA->Write("GemY_vs_StrawTypeA");
	GemY_vs_StrawTypeB->Write("GemY_vs_StrawTypeB");
	out->Close();
}
