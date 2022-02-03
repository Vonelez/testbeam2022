#include <cstdio>
#define evbuilder_cxx

#include "hits.h"

#include <TCanvas.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <vector>
#include <array>

#include <algorithm>
#include <fstream>

using std::vector;
using std::array;

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

struct track
{
    char straw_type;
    int straw_ch;
    int straw_hit_id;
    vector <array<int, 3>> gemCluster;
    double gem_wm_ch; //the weighted mean ch of cluster in GEM
    bool scintillator;
    long double strawT;
    long double sciT;
    int sci_hit_id;
};

int hits::gemChConverter(int ch, int vmm)
{
    int out_vmm = 0;
    if (vmm > 1)
    {
        out_vmm = vmm - 2;
    }
    else
    {
        out_vmm = vmm + 2;
    }
    int out = 64 * out_vmm + ch;
    return out;
}

void hits::Loop()
{
    // 1522 -- 31.1; 21.4
	// 1606 -- 147.3; 55.8
    double meanStrawGem = 147.3;
    double sigmaStrawGem = 55.8;

    // 1522 -- -103.4; 14.2
	// 1606 -- 11.0; 52.5
    double meanStrawScint = 11.0;
    double sigmaStrawScint = 52.5;

    bool jinrScint = false; // vmm8 ch63 = true

    int vmm_check = 0;
    int ch_check = 0;

    int strawCh = 0;
    int strawHitId = 0;
    long double strawT = 0;
    long double gemT = 0;
    long double sciT = 0;
    int gemCh = 0;
    int sciHitId = 0;

    if (jinrScint)
    {
        vmm_check = 8;
        ch_check = 63;
    }
    else
    {
        vmm_check = 9;
        ch_check = 0;
    }

    vector<track> tracks;

    auto *gem_strawA_correlarion_all = new TH2D("gem_strawA_correlarion_all", "gem_strawA_correlarion_all; StrawTypeA, ch; GEM3 Y-plane, ch", 64, 0, 64, 256, 0, 256);
    auto *gem_strawB_correlarion_all = new TH2D("gem_strawB_correlarion_all", "gem_strawB_correlarion_all; StrawTypeB, ch; GEM3 Y-plane, ch", 64, 0, 64, 256, 0, 256);
    auto *gem_strawA_correlarion_sci_only = new TH2D("gem_strawA_correlarion_sci_only", "gem_strawA_correlarion_sci_only; StrawTypeA, ch; GEM3 Y-plane, ch", 64, 0, 64, 256, 0, 256);
    auto *gem_strawB_correlarion_sci_only = new TH2D("gem_strawB_correlarion_sci_only", "gem_strawB_correlarion_sci_only; StrawTypeB, ch; GEM3 Y-plane, ch", 64, 0, 64, 256, 0, 256);

    auto *two_det_trackA_rate = new TH1D("two_det_trackA_rate", "two_det_trackA_rate; spills, sec; N", 900, 0, 900);
    auto *two_det_trackB_rate = new TH1D("two_det_trackB_rate", "two_det_trackB_rate; spills, sec; N", 900, 0, 900);

    auto *three_det_trackA_rate = new TH1D("three_det_trackA_rate", "three_det_trackA_rate; spills, sec; N", 900, 0, 900);
    auto *three_det_trackB_rate = new TH1D("three_det_trackB_rate", "three_det_trackB_rate; spills, sec; N", 900, 0, 900);

    auto *RT_curve_SCI = new TH2D("RT_curve_SCI", "RT with t_{0} from SCI and coord from GEM; R, mm; #Delta t, ns", 15, 0, 6, 500, -200., 300.);

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

            bool TypeA = false;
            bool TypeB = false;

            if (hits_fec[i] == 2 && hits_vmm[i] == 10)
            {
                if (hits_ch[i] % 4 == 0)
                {
                    TypeA = true;
                }
                else if (hits_ch[i] % 4 == 3)
                {
                    TypeB = true;
                }
                else
                    continue;

                strawCh = (int)hits_ch[i];
                strawT = (long double)hits_time[i];
                strawHitId = (int)hits_id[i];

                gemT = 0;
                gemCh = 0;

                vector <array<int, 3>> GemClusterA;
                vector <array<int, 3>> GemClusterB;

                for (int j = 0; j < hits_; j++)
                {
                    if (!hits_over_threshold[j])
                        continue;

                    if (i == j)
                        continue;

                    if (abs(strawT - (long double)hits_time[j]) > 300)
                        continue;

                    if (hits_det[j] == 3 && hits_plane[j] == 1)
                    {
                        gemT = (long double)hits_time[j];
                        if (gemT - strawT > meanStrawGem + 4 * sigmaStrawGem ||
                            gemT - strawT < meanStrawGem - 4 * sigmaStrawGem)
                        {
                            continue;
                        }
                        else
                        {
                            int gemCh = gemChConverter((int)hits_ch[j], (int)hits_vmm[j]);
                            if (TypeA)
                            {
                                if (abs(strawCh * 3.75 - 69.25 - gemCh) > 30)
                                    continue;
                                else
                                {
                                    array<int, 3> pair = { {gemCh, (int)hits_adc[j], (int)hits_id[j]} };
                                    GemClusterA.push_back(pair);
                                } 
                                
                            }
                            else if (TypeB)
                            {   
                                if (abs(strawCh * 3.75 - 76.25 - gemCh) > 30) 
                                    continue;
                                else
                                {
                                    array<int, 3> pair = { {gemCh, (int)hits_adc[j], (int)hits_id[j]} };
                                    GemClusterB.push_back(pair);
                                } 
                            }
                            else
                                continue;
                        }
                    }
                }
                double gem_wmCh;
                if (GemClusterA.size() > 0)
                {
                    double sum = 0;
		            double w_sum = 0;
                    for (int l = 0; l < GemClusterA.size(); l++)
                    {
                        double weight = GemClusterA[l][1] / 1024.0;
                        sum += GemClusterA[l][0] * weight;
			            w_sum += weight;
                    }
                    gem_wmCh = (int)(sum / w_sum * 100.0) / 100.0;
                }
                else if (GemClusterB.size() > 0)
                {
                    double sum = 0;
		            double w_sum = 0;
                    for (int l = 0; l < GemClusterB.size(); l++)
                    {
                        double weight = GemClusterB[l][1] / 1024.0;
                        sum += GemClusterB[l][0] * weight;
			            w_sum += weight;
                    }
                    gem_wmCh = sum / w_sum;
                }
                else 
                    continue;

                sciT = 0;
                int sciCount = 0;

                for (int j = 0; j < hits_; j++)
                {
                    if (!hits_over_threshold[j])
                        continue;

                    if (i == j)
                        continue;

                    if (abs(strawT - (long double)hits_time[j]) > 300)
                        continue;

                    if (hits_vmm[j] == vmm_check && hits_ch[j] == ch_check)
                    {
                        if ((long double)hits_time[j] - strawT > meanStrawScint + 4 * sigmaStrawScint ||
                            (long double)hits_time[j] - strawT < meanStrawScint - 4 * sigmaStrawScint)
                        {
                            continue;
                        }
                        else
                        {
                            sciCount++;
                            sciT += (long double)hits_time[j];
                            sciHitId = (int)hits_id[j];
                        }
                    }
                }
                if (sciCount > 0)
                {
                    sciT /= sciCount;
                    if (TypeA)
                    {
                        tracks.push_back({'A', strawCh, strawHitId, GemClusterA, gem_wmCh, true, strawT, sciT, sciHitId});
                    }
                    else if (TypeB)
                    {
                        tracks.push_back({'B', strawCh, strawHitId, GemClusterB, gem_wmCh, true, strawT, sciT, sciHitId});
                    }
                    else
                        continue;
                }
                else
                {
                    if (TypeA)
                    {
                        tracks.push_back({'A', strawCh, strawHitId, GemClusterA, gem_wmCh, false, strawT, -1, -1});
                    }
                    else if (TypeB)
                    {
                        tracks.push_back({'B', strawCh, strawHitId, GemClusterB, gem_wmCh, false, strawT, -1, -1});
                    }
                    else
                        continue;
                }
            }   
            else
                continue;
        }
    }

    std::cout << "\t -----> " << tracks.size() << "\n";
    int countA = 0;
    int countB = 0;
    for (int i = 0; i < tracks.size(); i++)
    {
        track tmpTrack = tracks[i];

        if (tmpTrack.straw_type == 'A')
        {
            gem_strawA_correlarion_all->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
        }
        else
        {
            gem_strawB_correlarion_all->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
        }
        if (!tmpTrack.scintillator)
        {   
            if (tmpTrack.straw_type == 'A')
                two_det_trackA_rate->Fill(tmpTrack.strawT / 1e9);
            else
                two_det_trackB_rate->Fill(tmpTrack.strawT / 1e9);
            continue;
        }
        double hitCoord = 0;
        double u = 0;
        int v = (int)tmpTrack.gem_wm_ch;
        if (tmpTrack.straw_type == 'A')
        {
            countA++;
            u = ((v - 0) % 15);
            gem_strawA_correlarion_sci_only->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
            three_det_trackA_rate->Fill(tmpTrack.strawT / 1e9);
        }
        else
        {
            countB++;
            u = ((v - 4) % 15);
            gem_strawB_correlarion_sci_only->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
            three_det_trackB_rate->Fill(tmpTrack.strawT / 1e9);
        }
        hitCoord = u * 0.4;
        if (hitCoord < 0 || hitCoord > 6)
            continue;
        RT_curve_SCI->Fill(hitCoord, tmpTrack.strawT - tmpTrack.sciT);
    }

    std::cout << "Type A with SCINT " << countA << "\n";
    std::cout << "Type B with SCINT " << countB << "\n";

    ofstream file_tracks;
	file_tracks.open("txtFiles/tracks_" + file + ".txt");

	file_tracks << center("N", 5)                   << " | "
                << center("STRAW TYPE", 10)         << " | "
                << center("STRAW CH", 10)           << " | "
                << center("STRAW ID", 10)           << " | "
                << center("SCI ID", 10)             << " | "
                << center("GEM3Y claster CHs", 20)  << " | "
                << center("GEM3Y claster IDs", 20)  << " | "
                << center("GEM3Y wmCH", 15)         << "\n";

	file_tracks << std::string(105 + 3 * 6, '-') << "\n";

    for (int i = 0; i < tracks.size(); i++)
    {
        track tmpTrack = tracks[i];
        file_tracks << prd(i+1, 0, 5)                           << " | "
                    << prd(66 - tmpTrack.straw_type, 0, 10)    << " | "
                    << prd(tmpTrack.straw_ch, 0, 10)           << " | "
                    << prd(tmpTrack.straw_hit_id, 0, 10)       << " | "
                    << prd(tmpTrack.sci_hit_id, 0, 10)         << " | "
                    << prd(tmpTrack.gemCluster[0][0], 0, 20)   << " | "
                    << prd(tmpTrack.gemCluster[0][2], 0, 20)   << " | "
                    << prd(tmpTrack.gem_wm_ch, 0, 15)          << "\n";
        
        for (int j = 1; j < tmpTrack.gemCluster.size(); j++)
        {
            file_tracks << std::string(5, ' ')                      << " | "
                        << std::string(10, ' ')                     << " | "
                        << std::string(10, ' ')                     << " | "
                        << std::string(10, ' ')                     << " | "
                        << std::string(10, ' ')                     << " | "
                        << prd(tmpTrack.gemCluster[j][0], 0, 20)   << " | "
                        << prd(tmpTrack.gemCluster[j][2], 0, 20)   << " | "
                        << std::string(15, ' ')                     << "\n";
        }
        file_tracks << std::string(105 + 3 * 6, '-') << "\n";
    }

	// char straw_type;
    // int straw_ch;
    // int straw_hit_id;
    // vector <array<int, 3>> gemCluster;
    // double gem_wm_ch; //the weighted mean ch of cluster in GEM
    // bool scintillator;
    // long double strawT;
    // long double sciT;
    // int sci_hit_id;

    TFile *out = new TFile("out/evBuildPrelim" + file + ending, "RECREATE");
    gem_strawA_correlarion_all->Write("gem_strawA_correlarion_all");
    gem_strawB_correlarion_all->Write("gem_strawB_correlarion_all");
    gem_strawA_correlarion_sci_only->Write("gem_strawA_correlarion_sci_only");
    gem_strawB_correlarion_sci_only->Write("gem_strawB_correlarion_sci_only");
    two_det_trackA_rate->Write("two_det_trackA_rate");
    two_det_trackB_rate->Write("two_det_trackB_rate");
    three_det_trackA_rate->Write("three_det_trackA_rate");
    three_det_trackB_rate->Write("three_det_trackB_rate");
    RT_curve_SCI->Write("RT_curve_SCI");
    out->Close();

    TCanvas *c1 = new TCanvas("c1", "Type-A straws", 1440, 900);
	c1->Divide(2, 2);
	c1->cd(1);
    gem_strawA_correlarion_all->Draw();
    c1->cd(2);
    gem_strawA_correlarion_sci_only->Draw();
    c1->cd(3);
    gPad->SetLogy();
    two_det_trackA_rate->Draw();
    c1->cd(4);
    gPad->SetLogy();
    three_det_trackA_rate->Draw();
	c1->SaveAs("img/TypeA_evB_" + file + ".pdf");

    TCanvas *c2 = new TCanvas("c2", "Type-B straws", 1440, 900);
	c2->Divide(2, 2);
	c2->cd(1);
    gem_strawB_correlarion_all->Draw();
    c2->cd(2);
    gem_strawB_correlarion_sci_only->Draw();
    c2->cd(3);
    gPad->SetLogy();
    two_det_trackB_rate->Draw();
    c2->cd(4);
    gPad->SetLogy();
    three_det_trackB_rate->Draw();
	c2->SaveAs("img/TypeB_evB_" + file + ".pdf");
}