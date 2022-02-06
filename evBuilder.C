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
    int strawVmm;
    int straw_hit_id;
    vector <array<int, 3>> gemCluster;
    double gem_wm_ch; //the weighted mean ch of cluster in GEM
    // double gem_av_t; //the weighted mean ch of cluster in GEM ---> need to implement
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

void hits::threePlotDrawF(TH1D* h1, TH1D* h2, TH1D* h3, TString name)
{
    h1->SetLineColor(kGreen -2);
	h2->SetLineColor(kMagenta);
	h3->SetLineColor(kBlack);

    TCanvas *three_plots = new TCanvas(name, name);
	three_plots->cd();

	h1 -> SetStats(0);
	h2 -> SetStats(0);
	h3 -> SetStats(0);

	h3->Draw();
	h2->Draw("SAME");
	h1->Draw("SAME");

    double m1 = meanGemScint;
    double s1 = sigmaGemScint;

    double m2 = meanStrawScint;
    double s2 = sigmaStrawScint;

    double m3 = meanStrawGem;
    double s3 = sigmaStrawGem;

	h1->Fit("gaus","","",m1 - 4*s1,m1 + 4*s1); //TScint - TGEM3
	h2->Fit("gaus","","",m2 - 4*s2,m2 + 4*s2); //TScint - TStraw
	h3->Fit("gaus","","",m3 - 4*s3,m3 + 4*s3); //TGEM3 - TStraw
	gStyle->SetOptFit(1111);
	//draw fit parameters as legends:
	Char_t ndf[80];
	Char_t sigma[80];
	Char_t mean[80];
	Char_t constant[80];
	auto legend = new TLegend(0.65, 0.9, 1.0, 0.75, "TScint - TGEM3");
	sprintf(ndf, "#chi^{2}/NDF = %.2f / %.2i", h1->GetFunction("gaus")->GetChisquare(), h1->GetFunction("gaus")->GetNDF());
	legend -> AddEntry(h1, ndf);
	sprintf(sigma, "#sigma = %.1f", h1->GetFunction("gaus")->GetParameter(2));
	legend -> AddEntry(h1, sigma);
	sprintf(mean, "Mean = %.1f", h1->GetFunction("gaus")->GetParameter(1));
	legend -> AddEntry(h1, mean);
	sprintf(constant, "Events under peak: %.f", h1->GetFunction("gaus")->Integral(m1 - 4*s1,m1 + 4*s1) / h1->GetBinWidth(1));
	legend -> AddEntry(h1, constant);
	legend -> Draw("same");

	Char_t ndf1[80];
	Char_t sigma1[80];
	Char_t mean1[80];
	Char_t constant1[80];
	auto legend1 = new TLegend(0.65, 0.75, 1.0, 0.60, "TScint - TStraw");
	sprintf(ndf1, "#chi^{2}/NDF = %.2f / %.2i", h2->GetFunction("gaus")->GetChisquare(), h2->GetFunction("gaus")->GetNDF());
	legend1 -> AddEntry(h2, ndf1);
	sprintf(sigma1, "#sigma = %.1f", h2->GetFunction("gaus")->GetParameter(2));
	legend1 -> AddEntry(h2, sigma1);
	sprintf(mean1, "Mean = %.1f", h2->GetFunction("gaus")->GetParameter(1));
	legend1 -> AddEntry(h2, mean1);
	sprintf(constant1, "Events under peak: %.f", h2->GetFunction("gaus")->Integral(m2 - 4*s2,m2 + 4*s2) / h2->GetBinWidth(1));
	legend1 -> AddEntry(h2, constant1);
	legend1 -> Draw("same");

	Char_t ndf2[80];
	Char_t sigma2[80];
	Char_t mean2[80];
	Char_t constant2[80];
	auto legend2 = new TLegend(0.65, 0.60, 1.0, 0.45, "TGEM3 - TStraw");
	sprintf(ndf2, "#chi^{2}/NDF = %.2f / %.2i", h3->GetFunction("gaus")->GetChisquare(), h3->GetFunction("gaus")->GetNDF());
	legend2 -> AddEntry(h3, ndf2);
	sprintf(sigma2, "#sigma = %.1f", h3->GetFunction("gaus")->GetParameter(2));
	legend2 -> AddEntry(h3, sigma2);
	sprintf(mean2, "Mean = %.1f", h3->GetFunction("gaus")->GetParameter(1));
	legend2 -> AddEntry(h3, mean2);
	sprintf(constant2, "Events under peak: %.f", h3->GetFunction("gaus")->Integral(m3 - 4*s3,m3 + 4*s3) / h3->GetBinWidth(1));
	legend2 -> AddEntry(h3, constant2);
	legend2 -> Draw("same");

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

    auto *h1_before = new TH1D("h1_before", "h1_before", 1000, -1000, 1000);
    auto *h1_2hit = new TH1D("h1_2hit", "h1_2hit", 1000, -1000, 1000);
    auto *h1_3hit = new TH1D("h1_3hit", "h1_3hit", 1000, -1000, 1000);

    auto *h2_before = new TH1D("h2_before", "h2_before", 1000, -1000, 1000);
    auto *h2_2hit = new TH1D("h2_2hit", "h2_2hit", 1000, -1000, 1000);
    auto *h2_3hit = new TH1D("h2_3hit", "h2_3hit", 1000, -1000, 1000);

    auto *h3_before = new TH1D("h3_before", "h3_before", 1000, -1000, 1000);
    auto *h3_2hit = new TH1D("h3_2hit", "h3_2hit", 1000, -1000, 1000);
    auto *h3_3hit = new TH1D("h3_3hit", "h3_3hit", 1000, -1000, 1000);

    auto *gem_strawA_correlarion_all = new TH2D("gem_strawA_correlarion_all", "gem_strawA_correlarion_all; StrawTypeA, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);
    auto *gem_strawB_correlarion_all = new TH2D("gem_strawB_correlarion_all", "gem_strawB_correlarion_all; StrawTypeB, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);
    auto *gem_strawA_correlarion_sci_only = new TH2D("gem_strawA_correlarion_sci_only", "gem_strawA_correlarion_sci_only; StrawTypeA, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);
    auto *gem_strawB_correlarion_sci_only = new TH2D("gem_strawB_correlarion_sci_only", "gem_strawB_correlarion_sci_only; StrawTypeB, ch; GEM3 Y-plane, ch", 128, 0, 128, 256, 0, 256);

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

                gemT = 0;
                gemCh = 0;

                vector <array<int, 3>> GemClusterA;
                vector <array<int, 3>> GemClusterB;

                vector <long double> GemTvector;

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
                            GemTvector.push_back(gemT);
                            h3_before->Fill(gemT - strawT);
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
                    h2_before->Fill(sciT - strawT);

                    for (int k = 0; k < GemTvector.size(); k++)
                    {
                        h1_before->Fill(sciT - GemTvector[k]);
                    }
                        


                    if (TypeA)
                    {
                        tracks.push_back({'A', strawCh, strawVmm, strawHitId, GemClusterA, gem_wmCh, true, strawT, sciT, sciHitId});
                    }
                    else if (TypeB)
                    {
                        tracks.push_back({'B', strawCh, strawVmm, strawHitId, GemClusterB, gem_wmCh, true, strawT, sciT, sciHitId});
                    }
                    else
                        continue;
                }
                else
                {
                    if (TypeA)
                    {
                        tracks.push_back({'A', strawCh, strawVmm, strawHitId, GemClusterA, gem_wmCh, false, strawT, -1, -1});
                    }
                    else if (TypeB)
                    {
                        tracks.push_back({'B', strawCh, strawVmm, strawHitId, GemClusterB, gem_wmCh, false, strawT, -1, -1});
                    }
                    else
                        continue;
                }
            }   
            else
                continue;
        }
    }

    threePlotDrawF(h1_before, h2_before, h3_before, "ALL_possuble_correlations");

    std::cout << "\t -----> " << tracks.size() << "\n";
    int countA = 0;
    int countB = 0;
    for (int i = 0; i < tracks.size(); i++)
    {
        track tmpTrack = tracks[i];

        if (tmpTrack.strawVmm == 11  && tmpTrack.straw_ch == 72)
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
        hitCoord = u * 0.4;
        if (hitCoord < 0 || hitCoord > 6)
            continue;
        if (ChCheck)
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
}