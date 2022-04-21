#include <cstdio>
#define tbjinr_cxx

#include "hits.h"

#include <TCanvas.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>

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
    // if (s1 < s2)
    // {
    //     h1->Draw();
    //     h2->Draw("SAME");
    //     h3->Draw("SAME");
    // }
    // else
    // {
    //     h2->Draw();
    //     h3->Draw("SAME");
    //     h1->Draw("SAME");
    // }
    h1->Draw();
    h2->Draw("SAME");
    h3->Draw("SAME");

    h1->Fit("gaus", "", "", 100, 250); // TScint - TGEM3
    h2->Fit("gaus", "", "", -200, -40); // TScint - TStraw
    h3->Fit("gaus", "", "", -150, 50); // TGEM3 - TStraw
    // gStyle->SetOptFit(1111);
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
    sprintf(constant, "Events under peak: %.f", h1->GetFunction("gaus")->Integral(100, 250) / h1->GetBinWidth(1));
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
    sprintf(constant1, "Events under peak: %.f", h2->GetFunction("gaus")->Integral(-200, -40) / h2->GetBinWidth(1));
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
    sprintf(constant2, "Events under peak: %.f", h3->GetFunction("gaus")->Integral(-150, 50) / h3->GetBinWidth(1));
    legend2->AddEntry(h3, constant2);
    legend2->Draw("same");

    three_plots->SaveAs("img/3plots_" + name + "_" + file + ".pdf");
}

struct track
{
    int straw_ch_a;
    int straw_hit_id_a;
    long double strawT_a;
    int straw_ch_b;
    int straw_hit_id_b;
    long double strawT_b;
    vector<array<int, 3> > gemCluster;
    double gem_wm_ch;
    long double gem_av_T;
    long double sciT;
    int sci_hit_id; 
};

void hits::Loop()
{
    if (fChain == 0)
        return;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    vector <track> prelimTrack;

    auto *bcid_sci = new TH1D("bcid_sci", "bcid_sci; BCID", 4200, 0, 4200);
    auto *tdc_sci = new TH1D("tdc_sci", "tdc_sci; TDC", 500, 0, 500);
    auto *adc_sci = new TH1D("adc_sci", "adc_sci; ADC", 1200, 0, 1200);
    auto *bcid_gem = new TH1D("bcid_gem", "bcid_gem; BCID", 4200, 0, 4200);
    auto *tdc_gem = new TH1D("tdc_gem", "tdc_gem; TDC", 500, 0, 500);
    auto *adc_gem = new TH1D("adc_gem", "adc_gem; ADC", 1200, 0, 1200);
    auto *bcid_straw = new TH1D("bcid_straw", "bcid_straw; BCID", 4200, 0, 4200);
    auto *tdc_straw = new TH1D("tdc_straw", "tdc_straw; TDC", 500, 0, 500);
    auto *adc_straw = new TH1D("adc_straw", "adc_straw; ADC", 1200, 0, 1200);
    auto *adc_strawA = new TH1D("adc_strawA", "adc_strawA; ADC", 1200, 0, 1200);
    auto *adc_strawB = new TH1D("adc_strawB", "adc_strawB; ADC", 1200, 0, 1200);
    auto *bcid_straw_ch_96 = new TH1D("bcid_straw_ch_96", "bcid_straw_ch_96; BCID", 4200, 0, 4200);
    auto *tdc_straw_ch_96 = new TH1D("tdc_straw_ch_96", "tdc_straw_ch_96; TDC", 500, 0, 500);
    auto *adc_straw_ch_96 = new TH1D("adc_straw_ch_96", "adc_straw_ch_96; ADC", 1200, 0, 1200);
    auto *spills = new TH1D("spills", "Num of spills; t, sec", 300, 0., 900.);
    auto *strawGemCorr = new TH1D("strawGemCorr", "strawGemCorr; #Deltat, nsec", 400, -200., 200.);
    auto *gem_straw_correlarion_all = new TH2D("gem_straw_correlarion_all", "gem_straw_correlarion_all; Straw, ch; GEM Y-plane, ch", 128, 0, 128, 128, 0, 128);
    auto *strawGem_Cluster_Corr = new TH1D("strawGem_Cluster_Corr", "strawGem_Cluster_Corr; #Deltat, nsec", 1000, -1000., 1000.);
    auto *strawSci_Corr = new TH1D("strawSci_Corr", "strawSci_Corr; #Deltat, nsec", 1000, -1000., 1000.);
    auto *sciGem_Cluster_Corr = new TH1D("sciGem_Cluster_Corr", "sciGem_Cluster_Corr; #Deltat, nsec", 1000, -1000., 1000.);
    auto *gem_Cluster_straw_correlarion_all = new TH2D("gem_Cluster_straw_correlarion_all", "gem_Cluster_straw_correlarion_all; Straw, ch; GEM Y-plane, ch", 128, 0, 128, 128, 0, 128);
    auto *strawGemCorrBcid = new TH1D("strawGemCorrBcid", "strawGemCorrBcid; #Delta, bcid", 10, -5, 5);
    auto *strawLife = new TH2D("strawLife", "strawLife; T, sec; Straw, ch", 400, 0, 400, 128, 0, 128);
    auto *strawAstrawB_correlarion = new TH2D("strawAstrawB_correlarion", "strawAstrawB_correlarion; StrawTypeA, ch; StrawTypeB, ch", 128, 0, 128, 128, 0, 128);
    auto *deltaTstrawAvsB = new TH1D("deltaTstrawAvsB", "deltaTstrawAvsB",1000, -500, 500);

    auto *strawADC_448_Check = new TH1D("strawADC_448_Check", "ADC = 448 check; straw ch",128, 0, 128);
    auto *strawADC_627_Check = new TH1D("strawADC_627_Check", "ADC = 627 check; straw ch",128, 0, 128);
    auto *strawADC_lastBins_Check = new TH1D("strawADC_lastBins_Check", "948 < ADC < 1023 check; straw ch",128, 0, 128);

    auto *strawASci_Corr = new TH1D("strawASci_Corr", "strawASci_Corr; #Deltat, nsec", 1000, -1000., 1000.);
    auto *strawBSci_Corr = new TH1D("strawBSci_Corr", "strawBSci_Corr; #Deltat, nsec", 1000, -1000., 1000.);

    auto *strawBCID_correlarion = new TH2D("strawBCID_correlarion", "t(sci - straw) vs straw BCID; straw BCID; #Delta t, nsec", 4200, 0, 4200, 1000, -500, 500);
    auto *sciBCID_correlarion = new TH2D("sciBCID_correlarion", "t(sci - straw) vs sci BCID; sci BCID; #Delta t, nsec", 4200, 0, 4200, 1000, -500, 500);
    auto *ABBCID_correlarion = new TH2D("ABBCID_correlarion", "straw A vs B BCID; straw A BCID; straw A BCID", 4200, 0, 4200, 4200, 0, 4200);

    auto *deltaTstrawAvsB_2H_correlation = new TH2D("deltaTstrawAvsB_2H_correlation", "#Delta t (scint - straw_{A}) vs #Delta t (scint - straw_{B});#Delta t_{A}, nsec;#Delta t_{B}, nsec",600, -300, 300, 600, -300, 300);
    
    auto *tdc_straw_100 = new TH1D("tdc_straw_100", "tdc_straw ch 100; TDC", 500, 0, 500);
    auto *tdc_straw_99 = new TH1D("tdc_straw_99", "tdc_straw ch 99; TDC", 500, 0, 500);

    auto *deltaTstrawAvsB_hitsTime = new TH1D("deltaTstrawAvsB_hitsTime", "deltaTstrawAvsB_hitsTime",400, -200, 200);
    auto *deltaTstrawAvsB_viaTDC = new TH1D("deltaTstrawAvsB_viaTDC", "deltaTstrawAvsB_viaTDC",400, -200, 200);

    auto *deltaBCIDstrawAvsB = new TH1D("deltaBCIDstrawAvsB", "deltaBCIDstrawAvsB",400, -200, 200);

    auto *deltaBCIDstrawSci = new TH1D("deltaBCIDstrawSci", "deltaBCIDstrawSci",40, -20, 20);


    int strawAcounter = 0;
    int withoutBcounter = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        long double sciT;
        int sciHitId;
        long double gemT;
        int gemCh;
        long double strawT;
        int strawCh;
        int strawHitId;
        int strawTdc;
        long double strawT_B;
        int strawCh_B;
        int strawHitId_B;
        int strawTdc_B;
        Long64_t ientry = LoadTree(jentry);
        
        int strawBCID = 0;
        int strawbBCID = 0;
        int sciBCID = 0;

        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;


        for (int i = 0; i < hits_; i++)
        {
            if (!hits_over_threshold[i])
                continue;
            if (hits_vmm[i] < 2)
            {
                strawCh = hits_vmm[i] * 64 + hits_ch[i];

                if (strawCh == 1)
                {
                    bcid_sci->Fill(hits_bcid[i]);
                    tdc_sci->Fill(hits_tdc[i]);
                    adc_sci->Fill(hits_adc[i]);
                }
                else if (strawCh != 0)
                {
                    bcid_straw->Fill(hits_bcid[i]);
                    tdc_straw->Fill(hits_tdc[i]);
                    adc_straw->Fill(hits_adc[i]);

                    if (strawCh % 4 == 0)
                    {
                        adc_strawA->Fill(hits_adc[i]);
                    }
                    else if (strawCh % 4 == 3)
                    {
                        adc_strawB->Fill(hits_adc[i]);
                    }

                    if (hits_adc[i] == 144)
                    {
                        strawADC_448_Check->Fill(strawCh);
                    }
                    else if (hits_adc[i] == 160)
                    {
                        strawADC_627_Check->Fill(strawCh);
                    }
                    else if (hits_bcid[i] >= 120 && hits_bcid[i] <= 256)
                    {
                        strawADC_lastBins_Check->Fill(strawCh);
                    }
                }

                if (strawCh == 19) 
                {
                    bcid_straw_ch_96->Fill(hits_bcid[i]);
                    tdc_straw_ch_96->Fill(hits_tdc[i]);
                    adc_straw_ch_96->Fill(hits_adc[i]);
                }
                

                strawLife->Fill(hits_time[i]/1e9, strawCh);

                if (strawCh % 4 != 0 || strawCh == 0) 
                {
                    continue;
                }

                strawAcounter++;
                strawT = hits_time[i];
                strawHitId = hits_id[i];
                spills->Fill(strawT / 1e9);
                strawBCID = (int)hits_bcid[i];
                strawTdc = (int)hits_tdc[i];
                


            

                vector<array<int, 3> > GemCluster;
                vector<long double > GemClusterT;
                strawT_B = 0;
                strawCh_B = 0;
                strawHitId_B = 0;
                strawTdc_B = 0;
                double minT = 300;

                for (int j = 0; j < hits_; j++)
                {
                    if (!hits_over_threshold[j])
                        continue;

                    if (i == j)
                        continue;

                    if (abs(strawT - (long double)hits_time[j]) > 300)
                        continue;
                    

                    if (hits_vmm[j] < 2)
                    {
                        int ch = hits_vmm[j] * 64 + hits_ch[j];
                        

                        if (ch % 4 == 3  && abs(strawCh - ch) < 12)
                        // if (ch % 4 == 3)
                        {
                            if (abs(strawT - (long double)hits_time[j]) < minT)
                            {
                                strawT_B = (long double)hits_time[j];
                                strawCh_B = ch;
                                strawHitId_B = hits_id[j];
                                strawbBCID = (int)hits_bcid[j];
                                strawTdc_B = (int)hits_tdc[j];
                                minT = abs(strawT - strawT_B);
                            }
                        }
                    }
                }

                if (minT > 200)
                {
                    strawHitId_B = 0;
                }

                if (strawHitId_B != 0) 
                {
                    strawAstrawB_correlarion->Fill(strawCh, strawCh_B);
                    deltaTstrawAvsB->Fill(strawT - strawT_B);
                    ABBCID_correlarion->Fill(strawBCID, strawbBCID);

                    if(strawCh == 100 && strawCh_B == 99)
                    {
                        tdc_straw_99->Fill(strawTdc_B);
                        tdc_straw_100->Fill(strawTdc);
                        deltaBCIDstrawAvsB->Fill(strawBCID - strawbBCID);

                        deltaTstrawAvsB_hitsTime->Fill(strawT - strawT_B);
                        // deltaTstrawAvsB_viaTDC->Fill(((strawBCID  * 25.0) + ((strawTdc - 48) * 25.0/16.0)) - ((strawbBCID  * 25.0) + ((strawTdc_B - 48) * 25.0/16.0)));
                        deltaTstrawAvsB_viaTDC->Fill((strawTdc - strawTdc_B) * 25.0 / 16);

                    }
                }

                if (strawCh < 75 || strawCh > 120)
                {
                    continue;
                }

                int flag = 0;
                int prevCh = 0;

                // for (int j = 0; j < hits_; j++)
                // {
                //     if (!hits_over_threshold[j])
                //         continue;

                //     if (i == j)
                //         continue;

                //     if (abs(strawT - (long double)hits_time[j]) > 300)
                //         continue;

                //     if (hits_vmm[j] > 1)
                //     {
                //         gemCh = hits_vmm[j] * 64 - 128 + hits_ch[j];
                //         gemT = hits_time[j];
                //         if (abs(strawT - gemT + 35) > 80)
                //         {
                //             continue;
                //         }
                //         if (flag == 0)
                //         {
                //             prevCh = gemCh;
                //             flag++;
                //         }
                //         if (abs(prevCh - gemCh) > 3)
                //         {
                //             continue;
                //         }
                //         else
                //         {
                //             prevCh = gemCh;
                //         }
                //         bcid_gem->Fill(hits_bcid[j]);
                //         tdc_gem->Fill(hits_tdc[j]);
                //         adc_gem->Fill(hits_adc[j]);
                //         strawGemCorr->Fill(strawT - gemT);
                //         gem_straw_correlarion_all->Fill(strawCh, gemCh);
                //         strawGemCorrBcid->Fill(hits_bcid[i] - hits_bcid[j]);

                //         array<int, 3> pair = {{gemCh, (int)hits_adc[j], (int)hits_id[j]}};
                //         GemCluster.push_back(pair);
                //         GemClusterT.push_back(gemT);
                //     }
                // }

                sciT = 0;
                double minSciT = 300;

                for (int j = 0; j < hits_; j++)
                {
                    if (!hits_over_threshold[j])
                        continue;

                    if (i == j)
                        continue;

                    if (hits_vmm[j] == 0 && hits_ch[j] == 1)
                    {
                        // if ((long double)hits_time[j] - strawT  < 0 && (long double)hits_time[j] - strawT > -300)
                        if (abs((long double)hits_time[j] - strawT) < minSciT)
                        {
                            minSciT = abs(strawT - (long double)hits_time[j]);
                            sciT = (long double)hits_time[j];
                            sciHitId = (int)hits_id[j];
                            sciBCID = (int)hits_bcid[j];
                        }
                    }
                }
                if(sciT != 0){
                    strawSci_Corr->Fill(sciT - strawT);
                    strawBCID_correlarion->Fill(strawBCID, sciT - strawT);
                    sciBCID_correlarion->Fill(sciBCID, sciT - strawT);
                    deltaBCIDstrawSci->Fill(strawBCID - sciBCID);
                    if (strawHitId_B != 0)
                    {
                        strawASci_Corr->Fill(sciT - strawT);
                        strawBSci_Corr->Fill(sciT - strawT_B);
                        deltaTstrawAvsB_2H_correlation->Fill(sciT - strawT, sciT - strawT_B);
                        
                    }
                }

                // double gem_wmCh = 0;
                // long double gemAvT = 0;
                // if (GemCluster.size() > 0)
                // {
                //     double sum = 0;
                //     double w_sum = 0;
                //     for (int l = 0; l < GemCluster.size(); l++)
                //     {
                //         gemAvT += GemClusterT[l];
                //         double weight = GemCluster[l][1] / 1024.0;
                //         sum += GemCluster[l][0] * weight;
                //         w_sum += weight;
                //     }
                //     gem_wmCh = sum / w_sum;
                //     gemAvT /= GemCluster.size();
                    
                //     gem_Cluster_straw_correlarion_all->Fill(strawCh, gem_wmCh);
                //     // strawGem_Cluster_Corr->Fill(strawT - gemAvT);
                //     if (sciT != 0)
                //     {
                //         strawSci_Corr->Fill(-strawT + sciT);
                //         sciGem_Cluster_Corr->Fill(gemAvT - sciT);
                //         strawGem_Cluster_Corr->Fill(strawT - gemAvT);
                //     }

                //     if (strawHitId_B != 0)
                //     {
                //         if (sciT != 0)
                //         {
                //             prelimTrack.push_back({strawCh, strawHitId, strawT, strawCh_B, strawHitId_B, strawT_B, GemCluster, gem_wmCh, gemAvT, sciT, sciHitId});
                //         }
                //         else
                //         {
                //             prelimTrack.push_back({strawCh, strawHitId, strawT, strawCh_B, strawHitId_B, strawT_B, GemCluster, gem_wmCh, gemAvT, -1, -1});
                //         }
                        
                //     }
                //     else
                //     {
                //         withoutBcounter++;
                //     }                
                // }

                



                
                

                
            }
        }
    }
    // std::cout << center("Num of hits in Type-A straw:", 40) << " | \t" << strawAcounter << std::endl;
    // std::cout << center("Num of Tracks:", 40) << " | \t" << prelimTrack.size() << std::endl;
    // std::cout << center("Num of Tracks w/o Type-B straw:", 40) << " | \t" << withoutBcounter << std::endl;

    // threePlotDrawF(sciGem_Cluster_Corr, strawSci_Corr, strawGem_Cluster_Corr, "three_hits_correlations");
    // strawSci_Corr->Fit("gaus", "", "", -200, -40); // TScint - TStraw


    TFile *out = new TFile("out/tbJinr" + file + ending, "RECREATE");
    spills->Write("spills");
    bcid_sci->Write("bcid_sci");
    tdc_sci->Write("tdc_sci");
    adc_sci->Write("adc_sci");
    bcid_straw->Write("bcid_straw");
    tdc_straw->Write("tdc_straw");
    adc_straw->Write("adc_straw");
    adc_strawA->Write("adc_strawA");
    adc_strawB->Write("adc_strawB");
    bcid_straw_ch_96->Write("bcid_straw_ch_96");
    tdc_straw_ch_96->Write("tdc_straw_ch_96");
    adc_straw_ch_96->Write("adc_straw_ch_96");
    strawADC_448_Check->Write("strawADC_448_Check");
    strawADC_627_Check->Write("strawADC_627_Check");
    strawADC_lastBins_Check->Write("strawADC_lastBins_Check");
    // bcid_gem->Write("bcid_gem");
    // tdc_gem->Write("tdc_gem");
    // adc_gem->Write("adc_gem");
    // strawGemCorr->Write("strawGemCorr");
    // gem_straw_correlarion_all->Write("gem_straw_correlarion_all");
    // strawGem_Cluster_Corr->Write("strawGem_Cluster_Corr");
    strawSci_Corr->Write("strawSci_Corr");
    strawASci_Corr->Write("strawASci_Corr");
    strawBSci_Corr->Write("strawBSci_Corr");
    // sciGem_Cluster_Corr->Write("sciGem_Cluster_Corr");
    // gem_Cluster_straw_correlarion_all->Write("gem_Cluster_straw_correlarion_all");
    deltaTstrawAvsB->Write("deltaTstrawAvsB");
    strawAstrawB_correlarion->Write("strawAstrawB_correlarion");
    // strawGemCorrBcid->Write("strawGemCorrBcid");
    strawBCID_correlarion->Write("strawBCID_correlarion");
    sciBCID_correlarion->Write("sciBCID_correlarion");
    ABBCID_correlarion->Write("ABBCID_correlarion");
    deltaTstrawAvsB_2H_correlation->Write("deltaTstrawAvsB_2H_correlation");
    strawLife->Write("strawLife");
    tdc_straw_99->Write("tdc_straw_99");
    tdc_straw_100->Write("tdc_straw_100");

    deltaTstrawAvsB_hitsTime->Write("deltaTstrawAvsB_hitsTime");
    deltaTstrawAvsB_viaTDC->Write("deltaTstrawAvsB_viaTDC");

    deltaBCIDstrawAvsB->Write("deltaBCIDstrawAvsB");
    deltaBCIDstrawSci->Write("deltaBCIDstrawSci");

    out->Close();

    // ofstream file_tracks;

    // file_tracks.open("txtFiles/tracks_" + file + ".txt");

    // file_tracks << center("N", 5) << " | "
    //             << center("STRAW A CH", 10) << " | "
    //             << center("STRAW A ID", 10) << " | "
    //             << center("STRAW B CH", 10) << " | "
    //             << center("STRAW B ID", 10) << " | "
    //             << center("SCI ID", 10)     << " | "
    //             << center("GEM3Y claster CHs", 20) << " | "
    //             << center("GEM3Y claster IDs", 20) << " | "
    //             << center("GEM3Y wmCH", 15) << "\n";
    
    // file_tracks << std::string(110 + 8 * 3, '-') << "\n";

    // int straw_ch_a;
    // int straw_hit_id_a;
    // long double strawT_a;
    // int straw_ch_b;
    // int straw_hit_id_b;
    // long double strawT_b;
    // vector<array<int, 3> > gemCluster;
    // double gem_wm_ch;
    // long double gem_av_T;

    // for (int i = 0; i < prelimTrack.size(); i++)
    // {
    //     track tmpTrack = prelimTrack[i];

    //     file_tracks << prd(i + 1, 0, 5) << " | " 
    //                 << prd(tmpTrack.straw_ch_a, 0, 10) << " | "
    //                 << prd(tmpTrack.straw_hit_id_a, 0, 10) << " | "
    //                 << prd(tmpTrack.straw_ch_b, 0, 10) << " | "
    //                 << prd(tmpTrack.straw_hit_id_b, 0, 10) << " | "
    //                 << prd(tmpTrack.sci_hit_id, 0, 10) << " | "
    //                 << prd(tmpTrack.gemCluster[0][0], 0, 20) << " | "
    //                 << prd(tmpTrack.gemCluster[0][2], 0, 20) << " | "
    //                 << prd(tmpTrack.gem_wm_ch, 0, 15) << "\n";

    //     for (int j = 1; j < tmpTrack.gemCluster.size(); j++)
    //     {
    //         file_tracks << std::string(5, ' ')  << " | "
    //                     << std::string(10, ' ') << " | "
    //                     << std::string(10, ' ') << " | "
    //                     << std::string(10, ' ') << " | "
    //                     << std::string(10, ' ') << " | " 
    //                     << std::string(10, ' ') << " | "
    //                     << prd(tmpTrack.gemCluster[j][0], 0, 20) << " | "
    //                     << prd(tmpTrack.gemCluster[j][2], 0, 20) << " | "
    //                     << std::string(15, ' ') << "\n";
    //     }
    //     file_tracks << std::string(110 + 8 * 3, '-') << "\n";
    // }
}
