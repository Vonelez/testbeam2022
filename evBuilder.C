#include <cstdio>
#define evbuilder_cxx

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
    long double sci_straw;
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
    double meanStrawGem = 31.1;
    double sigmaStrawGem = 21.4;

    double meanStrawScint = -103.4;
    double sigmaStrawScint = 14.2;

    bool jinrScint = true; // vmm8 ch63 = true

    int vmm_check = 0;
    int ch_check = 0;

    int strawCh = 0;
    long double strawT = 0;
    long double gemT = 0;
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

                vector <int> GemClusterA;
                vector <int> GemClusterB;

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
                                if (abs(strawCh * 2.62 - 4.83 - gemCh) > 20)
                                //p0                        =     -4.83457   +/-   1.40857     
                                //p1                        =      2.62137   +/-   0.0282988 
                                    continue;
                                else
                                {
                                    GemCluster.push_back(gemCh);
                                } 
                                
                            }
                            else if (TypeB)
                            {   
                                continue;
                                //p0                        =    -0.631874   +/-   1.54827     
                                //p1                        =       2.4745   +/-   0.0293748 
                            }
                            else
                                continue;

                        }
                    }
                }
                if (GemCluster.size() > 0)
                {
                    std::cout << "-----> STARW CH \t" << strawCh << std::endl;
                    for (int l = 0; l < GemCluster.size(); l++)
                    {
                        std::cout << "--> GEM CH \t" << GemCluster[l] << std::endl;
                    }
                }
            }
            else
                continue;
        }
    }
}