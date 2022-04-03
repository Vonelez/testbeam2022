#include <cstdio>
#define evbuilder_cxx

#include "hits.h"

#include <TCanvas.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <THStack.h>
#include <vector>
#include <array>

#include <algorithm>
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

struct track
{
    char straw_type;
    int straw_ch;
    int strawVmm;
    int straw_hit_id;
    vector<array<int, 3> > gemCluster;
    double gem_wm_ch; // the weighted mean ch of cluster in GEM
    long double gem_av_T;  // the average T of cluster in GEM
    bool scintillator;
    long double strawT;
    long double sciT;
    int sci_hit_id;
};

struct strawAB_plusSci_withGem
{
    int straw_ch_a;
    int strawVmm_a;
    int straw_hit_id_a;
    long double strawT_a;
    int straw_ch_b;
    int strawVmm_b;
    int straw_hit_id_b;
    long double strawT_b;
    bool up;
    vector<array<int, 3> > gemCluster;
    double gem_wm_ch;
    long double gem_av_T; 
    bool scintillator;
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

void hits::threePlotDrawF(TH1D *h1, TH1D *h2, TH1D *h3, TString name)
{
    h1->SetLineColor(kGreen - 2);
    h2->SetLineColor(kMagenta);
    h3->SetLineColor(kBlack);
    h1->SetStats(0);
    h2->SetStats(0);
    h3->SetStats(0);



    double m1 = meanGemScint;
    double s1 = sigmaGemScint;

    double m2 = meanStrawScint;
    double s2 = sigmaStrawScint;

    double m3 = meanStrawGem;
    double s3 = sigmaStrawGem;

    TCanvas *three_plots = new TCanvas(name, name);
    three_plots->cd();
    if (s1 < s2)
    {
        h1->Draw();
        h2->Draw("SAME");
        h3->Draw("SAME");
    }
    else
    {
        h2->Draw();
        h3->Draw("SAME");
        h1->Draw("SAME");
    }


    h1->Fit("gaus", "", "", m1 - 4 * s1, m1 + 4 * s1); // TScint - TGEM3
    h2->Fit("gaus", "", "", m2 - 4 * s2, m2 + 4 * s2); // TScint - TStraw
    h3->Fit("gaus", "", "", m3 - 4 * s3, m3 + 4 * s3); // TGEM3 - TStraw
    gStyle->SetOptFit(1111);
    // draw fit parameters as legends:
    Char_t ndf[80];
    Char_t sigma[80];
    Char_t mean[80];
    Char_t constant[80];
    auto legend = new TLegend(0.65, 0.9, 1.0, 0.75, "TScint - TGEM3");
    sprintf(ndf, "#chi^{2}/NDF = %.2f / %.2i", h1->GetFunction("gaus")->GetChisquare(), h1->GetFunction("gaus")->GetNDF());
    legend->AddEntry(h1, ndf);
    sprintf(sigma, "#sigma = %.1f", h1->GetFunction("gaus")->GetParameter(2));
    legend->AddEntry(h1, sigma);
    sprintf(mean, "Mean = %.1f", h1->GetFunction("gaus")->GetParameter(1));
    legend->AddEntry(h1, mean);
    sprintf(constant, "Events under peak: %.f", h1->GetFunction("gaus")->Integral(m1 - 4 * s1, m1 + 4 * s1) / h1->GetBinWidth(1));
    legend->AddEntry(h1, constant);
    legend->Draw("same");

    Char_t ndf1[80];
    Char_t sigma1[80];
    Char_t mean1[80];
    Char_t constant1[80];
    auto legend1 = new TLegend(0.65, 0.75, 1.0, 0.60, "TScint - TStraw");
    sprintf(ndf1, "#chi^{2}/NDF = %.2f / %.2i", h2->GetFunction("gaus")->GetChisquare(), h2->GetFunction("gaus")->GetNDF());
    legend1->AddEntry(h2, ndf1);
    sprintf(sigma1, "#sigma = %.1f", h2->GetFunction("gaus")->GetParameter(2));
    legend1->AddEntry(h2, sigma1);
    sprintf(mean1, "Mean = %.1f", h2->GetFunction("gaus")->GetParameter(1));
    legend1->AddEntry(h2, mean1);
    sprintf(constant1, "Events under peak: %.f", h2->GetFunction("gaus")->Integral(m2 - 4 * s2, m2 + 4 * s2) / h2->GetBinWidth(1));
    legend1->AddEntry(h2, constant1);
    legend1->Draw("same");

    Char_t ndf2[80];
    Char_t sigma2[80];
    Char_t mean2[80];
    Char_t constant2[80];
    auto legend2 = new TLegend(0.65, 0.60, 1.0, 0.45, "TGEM3 - TStraw");
    sprintf(ndf2, "#chi^{2}/NDF = %.2f / %.2i", h3->GetFunction("gaus")->GetChisquare(), h3->GetFunction("gaus")->GetNDF());
    legend2->AddEntry(h3, ndf2);
    sprintf(sigma2, "#sigma = %.1f", h3->GetFunction("gaus")->GetParameter(2));
    legend2->AddEntry(h3, sigma2);
    sprintf(mean2, "Mean = %.1f", h3->GetFunction("gaus")->GetParameter(1));
    legend2->AddEntry(h3, mean2);
    sprintf(constant2, "Events under peak: %.f", h3->GetFunction("gaus")->Integral(m3 - 4 * s3, m3 + 4 * s3) / h3->GetBinWidth(1));
    legend2->AddEntry(h3, constant2);
    legend2->Draw("same");

    three_plots->SaveAs("img/3plots_" + name + "_" + file + ".pdf");
}

void hits::Loop()
{

    int vmm_check = 0;
    int ch_check = 0;

    int strawCh = 0;
    int strawVmm = 0;
    int strawHitId = 0;
    long double strawT = 0;
    long double gemT = 0;
    long double sciT = 0;
    int gemCh = 0;
    int sciHitId = 0;

    bool ChCheck = false;

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
    vector<strawAB_plusSci_withGem> abCorrVect;

    auto *h_adc_straw = new TH1D("h_adc_straw", "h_adc_straw", 1200, 0, 1200);

    auto *deltaTstrawAvsB = new TH1D("deltaTstrawAvsB", "deltaTstrawAvsB",1000, -500, 500);
    auto *deltaTstrawAvsB_15ns = new TH1D("deltaTstrawAvsB_15ns", "deltaTstrawAvsB_15ns", 1000, -500, 500);
    auto *strawAstrawB_correlarion = new TH2D("strawAstrawB_correlarion", "strawAstrawB_correlarion; StrawTypeA, ch; StrawTypeB, ch", 128, 0, 128, 128, 0, 128);
    auto *bananaPlotUp = new TH2D("bananaPlotUp", "bananaPlotUp; StrawTypeA #Delta t, ns; StrawTypeB #Delta t, ns", 400, -200, 200, 400, -200, 200);
    auto *bananaPlotDown = new TH2D("bananaPlotDown", "bananaPlotDown; StrawTypeA #Delta t, ns; StrawTypeB #Delta t, ns", 400, -200, 200, 400, -200, 200);


    auto *h1_3hit = new TH1D("h1_3hit", "3 detectors track; #Deltat, ns; N", 1000, -1000, 1000);

    auto *h2_3hit = new TH1D("h2_3hit", "3 detectors track; #Deltat, ns; N", 1000, -1000, 1000);

    auto *h3_3hit = new TH1D("h3_3hit", "3 detectors track; #Deltat, ns; N", 1000, -1000, 1000);

    auto *gem_strawA_correlarion_all = new TH2D("gem_strawA_correlarion_all", "gem_strawA_correlarion_all; StrawTypeA, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);
    auto *gem_strawB_correlarion_all = new TH2D("gem_strawB_correlarion_all", "gem_strawB_correlarion_all; StrawTypeB, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);
    auto *gem_strawA_correlarion_sci_only = new TH2D("gem_strawA_correlarion_sci_only", "gem_strawA_correlarion_sci_only; StrawTypeA, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);
    auto *gem_strawB_correlarion_sci_only = new TH2D("gem_strawB_correlarion_sci_only", "gem_strawB_correlarion_sci_only; StrawTypeB, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);
    auto *gem_strawAB_correlarion_sci_only = new TH2D("gem_strawAB_correlarion_sci_only", "gem_strawAB_correlarion_sci_only; Straw, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);


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

            // if ((long double)hits_time[i] / 1e9 > 200)
                // continue;

            if (hits_fec[i] == 2 && (hits_vmm[i] == 10 || hits_vmm[i] == 11))
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

                if (hits_vmm[i] == 11)
                {
                    strawCh = (int)hits_ch[i] + 64;
                }
                else
                {
                    strawCh = (int)hits_ch[i];
                }
                strawT = (long double)hits_time[i];
                strawHitId = (int)hits_id[i];
                strawVmm = (int)hits_vmm[i];
                h_adc_straw->Fill((int)hits_adc[i]);

                gemT = 0;
                gemCh = 0;

                vector<array<int, 3> > GemClusterA;
                vector<array<int, 3> > GemClusterB;

                // vector <long double> GemTvector;
                double gemAvT = 0;


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
                            gemAvT += gemT;

                            if (TypeA)
                            {
                                if (abs(strawCh * 3.75 - 69.25 - gemCh) > 30)
                                    continue;
                                else
                                {
                                    array<int, 3> pair = {{gemCh, (int)hits_adc[j], (int)hits_id[j]}};
                                    GemClusterA.push_back(pair);
                                }
                            }
                            else if (TypeB)
                            {
                                if (abs(strawCh * 3.75 - 76.25 - gemCh) > 30)
                                    continue;
                                else
                                {
                                    array<int, 3> pair = {{gemCh, (int)hits_adc[j], (int)hits_id[j]}};
                                    GemClusterB.push_back(pair);
                                }
                            }
                            else
                                continue;
                        }
                    }
                }

                double gem_wmCh;
                
                int straw_ch_b_temp;
                int strawVmm_b_temp;
                int straw_hit_id_b_temp;
                long double strawT_b_temp;
                bool up_temp;
                double minTime;

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
                    gem_wmCh = sum / w_sum;
                    gemAvT /= GemClusterA.size();

                    // loking for typeB straw for existing typeA straw
                    minTime = 300;
                    straw_ch_b_temp = 0;
                    strawVmm_b_temp = 0;
                    straw_hit_id_b_temp = 0;
                    strawT_b_temp = 0;
                    up_temp = false;
                    for (int j = 0; j < hits_; j++)
                    {
                        if (!hits_over_threshold[j])
                            continue;

                        if (i == j)
                            continue;

                        if (abs(strawT - (long double)hits_time[j]) > 300)
                            continue;

                        if (hits_fec[j] == 2 && (hits_vmm[j] == 10 || hits_vmm[j] == 11))
                        {
                            int ch = 0;
                            if (hits_vmm[j] == 11)
                            {
                                ch = (int)hits_ch[j] + 64;
                            }
                            else
                            {
                                ch = (int)hits_ch[j];
                            }
                            if (hits_ch[j] % 4 == 3 && abs(strawCh - ch) < 12)
                            {
                                deltaTstrawAvsB->Fill(strawT - (long double)hits_time[j]);
                                if (abs(strawCh - ch) < 12)
                                {
                                    strawAstrawB_correlarion->Fill(strawCh, ch);
                                    deltaTstrawAvsB_15ns->Fill(strawT - (long double)hits_time[j]);

                                    if (abs(strawT - (long double)hits_time[j]) < minTime)
                                    {
                                        minTime = abs(strawT - (long double)hits_time[j]);
                                        straw_ch_b_temp = ch;
                                        strawVmm_b_temp = (int)hits_vmm[j];
                                        straw_hit_id_b_temp = (int)hits_id[j];
                                        strawT_b_temp = (long double)hits_time[j];
                                        
                                        if (strawCh - ch < 0)
                                        {
                                            up_temp = true;
                                        }
                                        else
                                        {
                                            up_temp = false;
                                        }
                                    }
                                }

                            }
                        }
                    }
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
                    gemAvT /= GemClusterB.size();
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
                        tracks.push_back({'A', strawCh, strawVmm, strawHitId, GemClusterA, gem_wmCh, gemAvT, true, strawT, sciT, sciHitId});
                        if (straw_ch_b_temp != 0) 
                        {
                            abCorrVect.push_back({strawCh, strawVmm, strawHitId, strawT, straw_ch_b_temp, strawVmm_b_temp, straw_hit_id_b_temp, strawT_b_temp, up_temp, GemClusterA, gem_wmCh, gemAvT, true, sciT, sciHitId});
                        }
                    }
                    else if (TypeB)
                    {
                        tracks.push_back({'B', strawCh, strawVmm, strawHitId, GemClusterB, gem_wmCh, gemAvT, true, strawT, sciT, sciHitId});
                    }
                    else
                        continue;
                }
                else
                {
                    if (TypeA)
                    {
                        tracks.push_back({'A', strawCh, strawVmm, strawHitId, GemClusterA, gem_wmCh, gemAvT, false, strawT, -1, -1});
                        if (straw_ch_b_temp != 0)
                        {
                            abCorrVect.push_back({strawCh, strawVmm, strawHitId, strawT, straw_ch_b_temp, strawVmm_b_temp, straw_hit_id_b_temp, strawT_b_temp, up_temp, GemClusterA, gem_wmCh, gemAvT, false, -1, -1});
                        }
                    }
                    else if (TypeB)
                    {
                        tracks.push_back({'B', strawCh, strawVmm, strawHitId, GemClusterB, gem_wmCh, gemAvT, false, strawT, -1, -1});
                    }
                    else
                        continue;
                }
            }
            else
                continue;
        }
    }
    // int straw_ch_a;
    // int strawVmm_a;
    // int straw_hit_id_a;
    // long double strawT_a;
    // int straw_ch_b;
    // int strawVmm_b;
    // int straw_hit_id_b;
    // long double strawT_b;
    // bool up;
    // vector<array<int, 3> > gemCluster;
    // long double gem_av_T;
    // bool scintillator;
    // long double sciT;

    ofstream abCorFile;
    abCorFile.open("txtFiles/abCorFile" + file + ".txt");

    abCorFile << center("N", 5) << " | "
              << center("STRAW A CH", 10) << " | "
              << center("STRAW A ID", 10) << " | "
              << center("STRAW B CH", 10) << " | "
              << center("STRAW B ID", 10) << " | "
              << center("STRAW B UP", 10) << " | "
              << center("SCI ID", 10) << " | "
              << center("GEM3Y claster CHs", 20) << " | "
              << center("GEM3Y claster IDs", 20) << " | "
              << center("GEM3Y wmCH", 15) << "\n";

    abCorFile << std::string(120 + 10 * 3, '-') << "\n";

    for (int i = 0; i < abCorrVect.size(); i++)
    {
        strawAB_plusSci_withGem tmpTrack = abCorrVect[i];
        abCorFile << prd(i + 1, 0, 5) << " | "
                  << prd(tmpTrack.straw_ch_a, 0, 10) << " | "
                  << prd(tmpTrack.straw_hit_id_a, 0, 10) << " | "
                  << prd(tmpTrack.straw_ch_b, 0, 10) << " | "
                  << prd(tmpTrack.straw_hit_id_b, 0, 10) << " | "
                  << prd(tmpTrack.up, 0, 10) << " | "
                  << prd(tmpTrack.sci_hit_id, 0, 10) << " | "
                  << prd(tmpTrack.gemCluster[0][0], 0, 20) << " | "
                  << prd(tmpTrack.gemCluster[0][2], 0, 20) << " | "
                  << prd(tmpTrack.gem_wm_ch, 0, 15) << "\n";

        for (int j = 1; j < tmpTrack.gemCluster.size(); j++)
        {
            abCorFile << std::string(5, ' ') << " | "
                      << std::string(10, ' ') << " | "
                      << std::string(10, ' ') << " | "
                      << std::string(10, ' ') << " | "
                      << std::string(10, ' ') << " | "
                      << std::string(10, ' ') << " | "
                      << std::string(10, ' ') << " | "
                      << prd(tmpTrack.gemCluster[j][0], 0, 20) << " | "
                      << prd(tmpTrack.gemCluster[j][2], 0, 20) << " | "
                      << std::string(15, ' ') << "\n";
        }
        abCorFile << std::string(120 + 10 * 3, '-') << "\n";
    }

    for (int j = 0; j < abCorrVect.size(); j++)
    {
        strawAB_plusSci_withGem tmpTrack = abCorrVect[j];
        if (tmpTrack.scintillator)
        {
            if (tmpTrack.up)
            {
                bananaPlotUp->Fill(tmpTrack.strawT_a - tmpTrack.sciT, tmpTrack.strawT_b - tmpTrack.sciT);
            }
            else
            {
                bananaPlotDown->Fill(tmpTrack.strawT_a - tmpTrack.sciT, tmpTrack.strawT_b - tmpTrack.sciT);
            }
        }
    }

    std::cout << "\t -----> " << tracks.size() << "\n";
    int countA = 0;
    int countB = 0;
    for (int i = 0; i < tracks.size(); i++)
    {
        track tmpTrack = tracks[i];

        if (tmpTrack.strawVmm == 11 && tmpTrack.straw_ch == 72)
        {
            ChCheck = true;
        }
        else
        {
            ChCheck = true;
        }

        if (tmpTrack.straw_type == 'A')
        {
            if (ChCheck)
                gem_strawA_correlarion_all->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
        }
        else
        {
            if (ChCheck)
                gem_strawB_correlarion_all->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
        }
        if (!tmpTrack.scintillator)
        {
            if (tmpTrack.straw_type == 'A')
            {
                if (ChCheck)
                    two_det_trackA_rate->Fill(tmpTrack.strawT / 1e9);
            }
            else
            {
                if (ChCheck)
                    two_det_trackB_rate->Fill(tmpTrack.strawT / 1e9);
            }

            continue;
        }

        double hitCoord = 0;
        double u = 0;
        int v = (int)tmpTrack.gem_wm_ch;

        h1_3hit->Fill(tmpTrack.sciT - tmpTrack.gem_av_T);
        h2_3hit->Fill(tmpTrack.sciT - tmpTrack.strawT);
        h3_3hit->Fill(tmpTrack.gem_av_T - tmpTrack.strawT);

        if (tmpTrack.straw_type == 'A')
        {
            countA++;
            u = ((v - 0) % 15);
            if (ChCheck)
            {
                gem_strawA_correlarion_sci_only->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
                three_det_trackA_rate->Fill(tmpTrack.strawT / 1e9);
            }
        }
        else
        {
            countB++;
            u = ((v - 4) % 15);
            if (ChCheck)
            {
                gem_strawB_correlarion_sci_only->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
                three_det_trackB_rate->Fill(tmpTrack.strawT / 1e9);
            }
        }
        gem_strawAB_correlarion_sci_only->Fill(tmpTrack.straw_ch, tmpTrack.gem_wm_ch);
        hitCoord = u * 0.4;
        if (hitCoord < 0 || hitCoord > 6)
            continue;
        if (ChCheck)
            RT_curve_SCI->Fill(hitCoord, tmpTrack.strawT - tmpTrack.sciT);
    }
    
    threePlotDrawF(h1_3hit, h2_3hit, h3_3hit, "three_hits_correlations");


    std::cout << "Type A with SCINT " << countA << "\n";
    std::cout << "Type B with SCINT " << countB << "\n";

    ofstream file_tracks;
    file_tracks.open("txtFiles/tracks_" + file + ".txt");

    file_tracks << center("N", 5) << " | "
                << center("STRAW TYPE", 10) << " | "
                << center("STRAW CH", 10) << " | "
                << center("STRAW ID", 10) << " | "
                << center("SCI ID", 10) << " | "
                << center("GEM3Y claster CHs", 20) << " | "
                << center("GEM3Y claster IDs", 20) << " | "
                << center("GEM3Y wmCH", 15) << "\n";

    file_tracks << std::string(105 + 3 * 6, '-') << "\n";

    for (int i = 0; i < tracks.size(); i++)
    {
        track tmpTrack = tracks[i];
        file_tracks << prd(i + 1, 0, 5) << " | "
                    << prd(66 - tmpTrack.straw_type, 0, 10) << " | "
                    << prd(tmpTrack.straw_ch, 0, 10) << " | "
                    << prd(tmpTrack.straw_hit_id, 0, 10) << " | "
                    << prd(tmpTrack.sci_hit_id, 0, 10) << " | "
                    << prd(tmpTrack.gemCluster[0][0], 0, 20) << " | "
                    << prd(tmpTrack.gemCluster[0][2], 0, 20) << " | "
                    << prd(tmpTrack.gem_wm_ch, 0, 15) << "\n";

        for (int j = 1; j < tmpTrack.gemCluster.size(); j++)
        {
            file_tracks << std::string(5, ' ') << " | "
                        << std::string(10, ' ') << " | "
                        << std::string(10, ' ') << " | "
                        << std::string(10, ' ') << " | "
                        << std::string(10, ' ') << " | "
                        << prd(tmpTrack.gemCluster[j][0], 0, 20) << " | "
                        << prd(tmpTrack.gemCluster[j][2], 0, 20) << " | "
                        << std::string(15, ' ') << "\n";
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
    h_adc_straw->Write("straw_adc");
    deltaTstrawAvsB->Write("deltaTstrawAvsB");
    strawAstrawB_correlarion->Write("strawAstrawB_correlarion");
    gem_strawAB_correlarion_sci_only->Write("gem_strawAB_correlarion_sci_only");
    bananaPlotUp->Write("bananaPlotUp"); 
    bananaPlotDown->Write("bananaPlotDown"); 
    h2_3hit->Write("straw_scint");
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
    c1->SaveAs("img/TypeA_evB_vmm9" + file + ".pdf");

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
    c2->SaveAs("img/TypeB_evB_vmm9" + file + ".pdf");


    deltaTstrawAvsB->Fit("gaus", "", "", -200, 200); 
    deltaTstrawAvsB_15ns->Fit("gaus", "", "", -200, 200); 
    deltaTstrawAvsB->SetLineColor(kGreen - 2);
    deltaTstrawAvsB_15ns->SetLineColor(kMagenta);

    
    TCanvas *c3 = new TCanvas("c3", "Type-B straws", 900, 1400);
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    c3->Divide(1, 2);
    c3->cd(1);
    deltaTstrawAvsB->Draw();
    deltaTstrawAvsB_15ns->Draw("SAME");
    Char_t ndf[80];
    Char_t sigma[80];
    Char_t mean[80];
    auto legend = new TLegend(0.65, 0.9, 1.0, 0.75, "Any time and 4 neighbors");
    sprintf(ndf, "#chi^{2}/NDF = %.2f / %.2i", deltaTstrawAvsB->GetFunction("gaus")->GetChisquare(), deltaTstrawAvsB->GetFunction("gaus")->GetNDF());
    legend->AddEntry(deltaTstrawAvsB, ndf);
    sprintf(sigma, "#sigma = %.1f", deltaTstrawAvsB->GetFunction("gaus")->GetParameter(2));
    legend->AddEntry(deltaTstrawAvsB, sigma);
    sprintf(mean, "Mean = %.1f", deltaTstrawAvsB->GetFunction("gaus")->GetParameter(1));
    legend->AddEntry(deltaTstrawAvsB, mean);
    legend->Draw("same");
    Char_t ndf1[80];
    Char_t sigma1[80];
    Char_t mean1[80];
    auto legend1 = new TLegend(0.65, 0.75, 1.0, 0.60, "Any time and 2 neighbor");
    sprintf(ndf1, "#chi^{2}/NDF = %.2f / %.2i", deltaTstrawAvsB_15ns->GetFunction("gaus")->GetChisquare(), deltaTstrawAvsB_15ns->GetFunction("gaus")->GetNDF());
    legend1->AddEntry(deltaTstrawAvsB_15ns, ndf1);
    sprintf(sigma1, "#sigma = %.1f", deltaTstrawAvsB_15ns->GetFunction("gaus")->GetParameter(2));
    legend1->AddEntry(deltaTstrawAvsB_15ns, sigma1);
    sprintf(mean1, "Mean = %.1f", deltaTstrawAvsB_15ns->GetFunction("gaus")->GetParameter(1));
    legend1->AddEntry(deltaTstrawAvsB_15ns, mean1);
    legend1->Draw("same");
    c3->cd(2);
    strawAstrawB_correlarion->Draw("COLZ");
    c3->SaveAs("img/strawAstrawB_correlarion" + file + ".pdf");

}