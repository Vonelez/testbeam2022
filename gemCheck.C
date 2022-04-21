#define gemCheck_cxx

#include "hits.h"

#include <TCanvas.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <cstdio>
#include <fstream>

using std::array;
using std::vector;

std::string prd(const double x, const int decDigits, const int width)
{
    stringstream ss;
    ss << fixed << right;
    ss.fill(' ');            // fill space around displayed #
    ss.width(width);         // set  width around displayed #
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
    if (padding > 0 && padding % 2 != 0)     // if odd #, add 1 space
        ss << " ";
    return ss.str();
};


void hits::Loop()
{
    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    auto *spills = new TH1D("spills", "Num of spills; t, sec", 300, 0., 3200.);
    auto *palneXvsYdeltaT = new TH1D("palneXvsYdeltaT", "Plane X vs Plane Y <#Delta t>; #Delta t, nsec", 300, -100., 200.);
    auto *gemXYgeomCheck = new TH2D("gemXYgeomCheck", "Plane X vs Plane Y; Plane X, ch; Plane Y, ch", 256, 0., 128., 256, 0., 128.);


    vector< vector <int> > clusterX;
    vector< vector <int> > clusterY;



    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        
        Long64_t ientry = LoadTree(jentry);

        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        vector <int> prelimX;
        vector <int> prelimY;

        
        for (int i = 0; i < hits_; i++)
        {
            prelimY.clear();

            if (!hits_over_threshold[i])
                continue;
            
            if (hits_vmm[i] == 1)
            {
                spills->Fill(hits_time[i] / 1e9);
            }
            int clusterLastIndex = 0;

            if (hits_vmm[i] == 2 ||  hits_vmm[i] == 3)
            {
                int gemYCh = 64 * (hits_vmm[i] - 2) + hits_ch[i];
                if (gemYCh == 0) continue;
                int gemXCh = 0;

                prelimY.push_back(gemYCh);
                for (int j = 0; j < hits_; j++)
                {
                    if (!hits_over_threshold[i])
                    continue;

                    if (i == j)
                    continue;

                    if (hits_vmm[j] == 2 ||  hits_vmm[j] == 3)
                    {
                        int gemYCh_cl = 64 * (hits_vmm[j] - 2) + hits_ch[j];
                        if (abs(hits_time[i] - hits_time[j]) > 5)
                        {
                            continue;
                        }
                        if (abs(gemYCh_cl - gemYCh) > 5)
                        {
                            continue;
                        }
                        
                        if (gemYCh_cl == 0) continue;
                        clusterLastIndex = j;
                        prelimY.push_back(gemYCh_cl);
                    }
                }


                
                for (int j = 0; j < hits_; j++)
                {
                    prelimX.clear();
                    
                    if (!hits_over_threshold[i] || i == j)
                    continue;

                    if (hits_vmm[j] == 0 ||  hits_vmm[j] == 1)
                    {
                        if (abs(hits_time[i] - hits_time[j]) > 10)
                        {
                            continue;
                        }
                        int gemXCh = 64 * hits_vmm[j] + hits_ch[j];
                        if (gemXCh == 0) continue;
                        int clusterXLastIndex = 0;
                        prelimX.push_back(gemXCh);

                        for (int k = 0; k < hits_; k++)
                        {
                            if (!hits_over_threshold[k])
                            {
                                continue;
                            }

                            if (j == k || i == k)
                            continue;

                            if (hits_vmm[k] == 0 ||  hits_vmm[k] == 1)
                            {
                                int gemXCh_cl = 64 * hits_vmm[k] + hits_ch[k];
                                if (abs(hits_time[k] - hits_time[j]) > 5)
                                {
                                    continue;
                                }
                                if (abs(gemXCh_cl - gemXCh) > 5)
                                {
                                    continue;
                                }
                                if (gemXCh_cl == 0) continue;
                                clusterXLastIndex = k;
                                prelimX.push_back(gemXCh_cl);
                            }
                        }
                        if (clusterXLastIndex != 0 && clusterXLastIndex > j){
                            j = clusterXLastIndex;
                        }

                    }
                }

                if (prelimX.size() != 0)
                {
                    clusterX.push_back(prelimX);
                    clusterY.push_back(prelimY);
                }


                if (clusterLastIndex != 0 && clusterLastIndex > i){
                    i = clusterLastIndex;
                }
            }
        }
    }


    for (int i = 0; i < clusterX.size(); i++)
    {
        vector <int> fX = clusterX.at(i);
        vector <int> fY = clusterY.at(i);

        double x_mean = 0, y_mean = 0;

        for (int j = 0; j < fX.size(); j++)
        {
            x_mean += fX[j];
        }
        x_mean /= fX.size();

        for (int j = 0; j < fY.size(); j++)
        {
            y_mean += fY[j];
        }
        y_mean /= fY.size();

        gemXYgeomCheck->Fill(x_mean, y_mean);
    }
    

    TFile *out = new TFile(file+ "_CHECK" + ending, "RECREATE");
    spills->Write("spills");
    gemXYgeomCheck->Write("gemXYgeomCheck");
    palneXvsYdeltaT->Write("palneXvsYdeltaT");
    out->Close();


    ofstream file_tracks;

    file_tracks.open(file + "_clusters" +  + ".txt");

    file_tracks << center("N", 5) << " | "
                << center("GEM Y CH", 10) << " | "
                << center("GEM X CH", 10) << "\n";
    
    file_tracks << std::string(25 + 2 * 3, '-') << "\n";

    for (int i = 0; i < clusterX.size(); i++)
    {
        vector <int> fX = clusterX.at(i);
        vector <int> fY = clusterY.at(i);

        file_tracks << prd(i + 1, 0, 5) << " | " 
                    << prd(fY[0], 0, 10) << " | "
                    << prd(fX[0], 0, 10) << "\n";

        for (int j = 1; j < fY.size(); j++)
        {
            file_tracks << std::string(5, ' ')  << " | "
                        << prd(fY[j], 0, 10) << " | "
                        << std::string(10, ' ') << "\n";
        }
        file_tracks << std::string(25 + 2 * 3, '-') << "\n";
    }
}
