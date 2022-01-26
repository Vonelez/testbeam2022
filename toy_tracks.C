#include <cstdio>
#define toy_tracks_cxx

#include "hits.h"

#include <TCanvas.h>
#include <TString.h>
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
	long double sci_straw_deltaT;
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
		64, 0, 64, 256, 0, 256);

	TH2D *GemY_vs_StrawTypeB = new TH2D(
		"GemY_vs_StrawTypeB", "GemY_vs_StrawTypeB; StrawTypeB, ch; GEM3 Y-plane, ch",
		64, 0, 64, 256, 0, 256);

	TH2D *RT_curve_GEM = new TH2D("RT_curve_GEM", "RT with T_0 from GEM and coord from GEM; R, mm; #Delta t, ns", 15, 0, 0.6, 250, 0., 0.);
	TH2D *RT_curve_SCI = new TH2D("RT_curve_SCI", "RT with T_0 from SCI and coord from GEM; R, mm; #Delta t, ns", 15, 0, 0.6, 250, 0., 0.);

	TH2D *GEM1 = new TH2D(
		"GEM1", "GEM1; CH; VMM", 32, 0, 64, 10, 0, 10);
	TH2D *GEM2 = new TH2D(
		"GEM2", "GEM2; CH; VMM", 32, 0, 64, 10, 0, 10);
	TH2D *GEM3 = new TH2D(
		"GEM3", "GEM3; CH; VMM", 32, 0, 64, 10, 0, 10);

	TH1D *bcid = new TH1D("bcid", "bcid; BCID", 4200, 0, 4200);
	TH1D *tdc = new TH1D("tdc", "tdc", 500, 0, 500);
	TH1D *adc = new TH1D("adc", "adc", 1200, 0, 1200);
	TH1D *spills = new TH1D("spills", "Num of spills; t, sec", 300, 0., 300.);

	TH1D *gem_straw_timeRes_slice1 = new TH1D("gem_straw_timeRes_slice1", "GEM3 Y-plane CHs:32-42; #Delta t, ns", 250, -500, 500);
	TH1D *gem_straw_timeRes_slice2 = new TH1D("gem_straw_timeRes_slice2", "GEM3 Y-plane CHs:82-92; #Delta t, ns", 250, -500, 500);
	// 1505 -- 238.4; 19.7
	// 1522 -- 31.1; 21.4
	// 1606 -- 147.3; 55.8
	// 1430 -- 238.9; 20.8
	double sigma_gem_straw = 20.8;
	double mean_gem_straw = 238.9;
	// 1505 -- -193.8; 18.9
	// 1522 -- -136.4; 17.1
	// 1606 -- -136.6; 17.0
	// 1430 -- -193.2; 19.1
	double sigma_gem_scint = 19.1;
	double mean_gem_scint = -193.2;
	int scint_check = 8;

	int vmm_check = 0;
	int ch_check = 0;

	if (scint_check == 8)
	{
		vmm_check = 8;
		ch_check = 63;
	}
	else
	{
		vmm_check = 9;
		ch_check = 0;
	}
	vector<pseudo_track> typeA;
	vector<pseudo_track> typeB;

	vector<vector<pseudo_track> > typeAbasedTracks;

	int gem_ch = 0;
	double gem_weight = 0;

	long double gem_scint_deltaT = 0;
	long double gem_straw_deltaT = 0;
	long double sci_straw_deltaT = 0;
	long double sciT = 0;

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
				// std::cout << "Det det " << (int)hits_det[i] << " det vmm (int)hits_vmm[i] " << (int)hits_vmm[i] << "\n";
				if (hits_det[i] == 1 || hits_det[i] == 2 || hits_det[i] == 3)
				{
					bcid->Fill((int)hits_bcid[i]);
					tdc->Fill((int)hits_tdc[i]);
					adc->Fill((int)hits_adc[i]);
					if (hits_det[i] == 1)
					{
						GEM1->Fill((int)hits_ch[i], (int)hits_vmm[i]);
					}
					else if (hits_det[i] == 2)
					{
						GEM2->Fill((int)hits_ch[i], (int)hits_vmm[i]);
					}
					else
					{
						GEM3->Fill((int)hits_ch[i], (int)hits_vmm[i]);
					}
				}
				else
				{
					continue;
				}
			}
			else
			{
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

						if (hits_fec[j] == 2 && hits_vmm[j] == vmm_check && hits_ch[j] == ch_check)
						{
							if ((long double)hits_time[j] - gemTime > mean_gem_scint + 4 * sigma_gem_scint ||
								(long double)hits_time[j] - gemTime < mean_gem_scint - 4 * sigma_gem_scint)
							{
								continue;
							}
							else
							{
								gem_scint_deltaT = gemTime - (long double)hits_time[j];
								sciT = (long double)hits_time[j];
								scinId = (int)hits_id[j];
								spills->Fill(sciT / 1e9);
							}
						}

						if (scinId == -1)
						{
							gem_scint_deltaT = -1.0;
							sci_straw_deltaT = -1.0;
							sciT = -1.0;
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
								gem_straw_deltaT = gemTime - (long double)hits_time[j];
								sci_straw_deltaT = sciT - (long double)hits_time[j];

								if (gem_ch > 32 && gem_ch < 42)
								{
									if (hits_ch[j] > 25 && hits_ch[j] < 30)
									{
										gem_straw_timeRes_slice1->Fill(gem_straw_deltaT);
									}
								}
								if (gem_ch > 82 && gem_ch < 92)
								{
									if (hits_ch[j] > 42 && hits_ch[j] < 44)
									{
										gem_straw_timeRes_slice2->Fill(gem_straw_deltaT);
									}
								}
								if (hits_ch[j] % 4 == 0)
								{
									int first_point_A = 14;
									GemY_vs_StrawTypeA->Fill((int)hits_ch[j], gem_ch);
									typeA.push_back({(int)hits_id[j], scinId, (int)hits_id[i], gem_ch, gem_weight, gem_scint_deltaT, gem_straw_deltaT, sci_straw_deltaT});
									// if (gem_ch > first_point_A)
									// {
									// 	for (int k = 0; k < 15; k++)
									// 	{
									// 		if ((gem_ch - first_point_A) % 15 == k )
									// 		{
									// 			RT_curve_GEM->Fill(k * 0.4, gem_straw_deltaT);
									// 			if (scinId > 0)
									// 			{
									// 				RT_curve_SCI->Fill(k * 0.4, sci_straw_deltaT);
									// 			}
									// 		}
									// 	}
									// }
								}
								else if (hits_ch[j] % 4 == 3)
								{
									GemY_vs_StrawTypeB->Fill((int)hits_ch[j], gem_ch);
									typeB.push_back({(int)hits_id[j], scinId, (int)hits_id[i], gem_ch, gem_weight, gem_scint_deltaT, gem_straw_deltaT, sci_straw_deltaT});
									int first_point_B = 4;
									if (gem_ch > first_point_B)
									{
										for (int k = 0; k < 15; k++)
										{
											if ((gem_ch - first_point_B) % 15 == k)
											{
												RT_curve_GEM->Fill(k * 0.4, gem_straw_deltaT);
												if (scinId > 0)
												{
													RT_curve_SCI->Fill(k * 0.4, sci_straw_deltaT);
												}
											}
										}
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
	}

	int typeA_all = 0, typeA_wScint = 0, typeA_woScint = 0;
	int typeB_all = 0, typeB_wScint = 0, typeB_woScint = 0;

	ofstream file_TypeA;
	file_TypeA.open("TypeA_" + file + ".txt");

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
	file_TypeB.open("TypeB_" + file + ".txt");

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
			vector<pseudo_track> tempTrack;
			tempTrack.push_back(typeA[j - 1]);
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
			vector<pseudo_track> tempTrack;
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
			sum += typeAbasedTracks[i][j].gem_ch * typeAbasedTracks[i][j].gem_weight;
			w_sum += typeAbasedTracks[i][j].gem_weight;
		}
		double w_mean = (int)(sum / w_sum * 100.0) / 100.0;

		// int first_point_A = 14;

		// if (w_mean > first_point_A)
		// {
		// 	for (int k = 0; k < 15; k++)
		// 	{
		// 		if (((int)w_mean - first_point_A) % 15 == k)
		// 		{
		// 			RT_curve_GEM->Fill(k * 0.4, gem_straw_deltaT);
		// 			if (scinId > 0)
		// 			{
		// 				RT_curve_SCI->Fill(k * 0.4, sci_straw_deltaT);
		// 			}
		// 		}
		// 	}
		// }
	}

	// printf("Straw id \t %d \t and w_mean \t %f \n", (int)typeAbasedTracks[i][0].straw_id, w_mean);
}

gStyle->SetOptFit();
gStyle->SetOptStat();

gem_straw_timeRes_slice1->Fit("gaus", "", "", 0., 0.);
gem_straw_timeRes_slice2->Fit("gaus", "", "", 0., 0.);

TCanvas *result_plots = new TCanvas("result_plots", "result_plots", 1440, 900);
result_plots->Divide(2, 2);
result_plots->cd(1);
GemY_vs_StrawTypeA->Draw("COLZ");
TLine l;
l.SetLineColor(kRed);
l.SetLineWidth(2);
l.DrawLine(0, 32, 64, 32);
l.DrawLine(0, 42, 64, 42);
l.DrawLine(0, 82, 64, 82);
l.DrawLine(0, 92, 64, 92);
result_plots->cd(2);
GemY_vs_StrawTypeB->Draw("COLZ");
TLine k;
k.SetLineColor(kRed);
k.SetLineWidth(2);
k.DrawLine(0, 32, 64, 32);
k.DrawLine(0, 42, 64, 42);
k.DrawLine(0, 82, 64, 82);
k.DrawLine(0, 92, 64, 92);
result_plots->cd(3);
gem_straw_timeRes_slice1->Draw();
result_plots->cd(4);
gem_straw_timeRes_slice2->Draw();
result_plots->SaveAs("resultPlots" + file + ".pdf");

TCanvas *curve_plots = new TCanvas("curve_plots", "curve_plots", 1440, 900);
curve_plots->Divide(2, 2);
curve_plots->cd(1);
GemY_vs_StrawTypeA->Draw("COLZ");
curve_plots->cd(2);
spills->Draw();
curve_plots->cd(3);
RT_curve_GEM->Draw("COLZ");
curve_plots->cd(4);
RT_curve_SCI->Draw("COLZ");
curve_plots->SaveAs("curve" + file + ".pdf");

TFile *out = new TFile("evBuildResults" + file + ending, "RECREATE");
GemY_vs_StrawTypeA->Write("GemY_vs_StrawTypeA");
GemY_vs_StrawTypeB->Write("GemY_vs_StrawTypeB");
bcid->Write("bcid");
tdc->Write("tdc");
adc->Write("adc");
GEM1->Write("GEM1");
GEM2->Write("GEM2");
GEM3->Write("GEM3");
gem_straw_timeRes_slice1->Write("slice1");
gem_straw_timeRes_slice2->Write("slice2");
RT_curve_SCI->Write("RT_curve_SCI");
RT_curve_GEM->Write("RT_curve_GEM");
spills->Write("spills");
out->Close();
}
