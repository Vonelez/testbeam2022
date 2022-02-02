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

struct pseudo_track
{
    int straw_id;
    int scint_id;
    int gem_id;
    int gem_ch;
    double gem_weight;
    long double gem_scint;
    long double gem_straw;
    long double sci_straw;
};

struct track
{
    char straw_type;
    int straw_ch;
    double gem_wm_ch; //the weighted mean ch of cluster in GEM
    bool scintillator;
    long double strawT;
    long double sciT;
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
    long double strawT = 0;
    long double gemT = 0;
    long double sciT = 0;
    int gemCh = 0;

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

    // auto *gem_strawA_correlarion = new TH2D("h2", "h2; StrawTypeA, ch; GEM3 Y-plane, ch", 64, 0, 64, 512, 0, 256);
    auto *RT_curve_SCI = new TH2D("RT_curve_SCI", "RT with t_{0} from SCI and coord from GEM; R, mm; #Delta t, ns", 60, 0, 6, 500, -200., 300.);

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

                gemT = 0;
                gemCh = 0;

                vector <std::array<int, 2>> GemClusterA;
                vector <std::array<int, 2>> GemClusterB;

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
                                if (abs(strawCh * 2.81 - 15.01 - gemCh) > 30)
                                // p0                        =     -15.0052   +/-   2.39285     
                                // p1                        =      2.80836   +/-   0.0455132
                                    continue;
                                else
                                {
                                    array<int, 2> pair = { {gemCh, (int)hits_adc[j]} };
                                    GemClusterA.push_back(pair);
                                } 
                                
                            }
                            else if (TypeB)
                            {   if (abs(strawCh * 2.81 + 3.72 - gemCh) > 30)
                                // p0                        =      3.72347   +/-   2.33071     
                                // p1                        =      2.39773   +/-   0.0427727  
                                    continue;
                                else
                                {
                                    array<int, 2> pair = { {gemCh, (int)hits_adc[j]} };
                                    GemClusterB.push_back(pair);
                                } 
                            }
                            else
                                continue;
                        }
                    }
                }
                char strawType;
                double gem_wmCh;
                if (GemClusterA.size() > 0)
                {
                    // std::cout << "-----> STARW CH \t" << strawCh << std::endl;
                    double sum = 0;
		            double w_sum = 0;
                    for (int l = 0; l < GemClusterA.size(); l++)
                    {
                        // std::cout << "--> GEM CH \t" << GemClusterA[l][0] << std::endl;
                        double weight = GemClusterA[l][1] / 1024.0;
                        sum += GemClusterA[l][0] * weight;
			            w_sum += weight;
                    }
                    gem_wmCh = (int)(sum / w_sum * 100.0) / 100.0;
                    // std::cout << "-> GEM AV \t" << (int)(sum / w_sum * 100.0) / 100.0 << std::endl;
                    strawType = 'A';
                    // gem_strawA_correlarion->Fill(strawCh, gem_wmCh);
                    // gem_strawA_correlarion->Draw("COLZ");
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
                    gem_wmCh = (int)(sum / w_sum * 100.0) / 100.0;
                    strawType = 'B';
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
                        }
                    }
                }
                if (sciCount > 0)
                {
                    sciT /= sciCount;
                    tracks.push_back({strawType, strawCh, gem_wmCh, true, strawT, sciT});

                }
                else
                    tracks.push_back({strawType, strawCh, gem_wmCh, false, strawT, -1.0});
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
        if (!tmpTrack.scintillator)
            continue;
        double hitCoord = 0;
        double u = 0;
        if (tmpTrack.straw_type == 'A')
        {
            countA++;
            u = (((int)tmpTrack.gem_wm_ch - 14) % 15) * 1.0;
        }
        else
        {
            countB++;
            u = (((int)tmpTrack.gem_wm_ch - 4) % 15) * 1.0;
        }
        if (hitCoord < 0 || hitCoord > 6)
            continue;
        hitCoord = (u + tmpTrack.gem_wm_ch - (int)tmpTrack.gem_wm_ch) * 0.4;
        RT_curve_SCI->Fill(hitCoord, tmpTrack.strawT - tmpTrack.sciT);
    }

    std::cout << "Type A with SCINT " << countA << "\n";
    std::cout << "Type B with SCINT " << countB << "\n";
    RT_curve_SCI->Draw("COLZ");

}