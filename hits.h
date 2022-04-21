//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct 17 15:24:01 2021 by ROOT version 6.22/00
// from TTree hits/hits
// found on file:
// run0055_bc_40p000_tac_60_ccs_7_cs_3_dt_100_mst_0_spc_500_dp_100_cr_0p50-2p00_coin_center-of-mass_smaller-files.root
//////////////////////////////////////////////////////////

#ifndef hits_h
#define hits_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.

class hits
{
public:
    // TString file = "107_at_peak_full";
    // TString file ="1543_at_threshold_full";
    // TString file = "200ns_40MHz_1700V";
    TString file = "gem/GEM_XY_FeSpot";
    TString ending = ".root";
    TTree *fChain;  //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    // 1505 -- 238.4; 19.7
    // 1522 -- 31.1; 21.4
	// 1606 -- 147.3; 55.8
    // 1543 -- 147.2; 55.4
    double meanStrawGem = 147.2;
    double sigmaStrawGem = 55.4;

    // 1505 -- 49.5; 15.7
    // 1522 -- -103.4; 14.2
	// 1606 -- 11.0; 52.5
    // 1543 -- 14.60; 54.5
    double meanStrawScint = 14.60;
    double sigmaStrawScint = 54.5;

    // 1505 -- -193.8; 18.9
    // 1522 -- -136.4; 17.1
    // 1606 -- -136.6; 17.0
    // 1543 -- -136.1; 17.5
    double meanGemScint = -136.1;
    double sigmaGemScint = 17.5;



    bool jinrScint = false; // vmm8 ch63 = true

    // Fixed size dimensions of array or collections stored in the TTree if any.
    static constexpr Int_t kMaxhits = 50000;

    // Declaration of leaf types
    Int_t hits_;
    UInt_t hits_id[kMaxhits];             //[hits_]
    UChar_t hits_det[kMaxhits];           //[hits_]
    UChar_t hits_plane[kMaxhits];         //[hits_]
    UChar_t hits_fec[kMaxhits];           //[hits_]
    UChar_t hits_vmm[kMaxhits];           //[hits_]
    Double_t hits_readout_time[kMaxhits]; //[hits_]
    Double_t hits_time[kMaxhits];         //[hits_]
    UChar_t hits_ch[kMaxhits];            //[hits_]
    UShort_t hits_pos[kMaxhits];          //[hits_]
    UShort_t hits_bcid[kMaxhits];         //[hits_]
    UShort_t hits_tdc[kMaxhits];          //[hits_]
    UShort_t hits_adc[kMaxhits];          //[hits_]
    Bool_t hits_over_threshold[kMaxhits]; //[hits_]
    Double_t hits_chip_time[kMaxhits];    //[hits_]

    // List of branches
    TBranch *b_hits_;               //!
    TBranch *b_hits_id;             //!
    TBranch *b_hits_det;            //!
    TBranch *b_hits_plane;          //!
    TBranch *b_hits_fec;            //!
    TBranch *b_hits_vmm;            //!
    TBranch *b_hits_readout_time;   //!
    TBranch *b_hits_time;           //!
    TBranch *b_hits_ch;             //!
    TBranch *b_hits_pos;            //!
    TBranch *b_hits_bcid;           //!
    TBranch *b_hits_tdc;            //!
    TBranch *b_hits_adc;            //!
    TBranch *b_hits_over_threshold; //!
    TBranch *b_hits_chip_time;      //!

    hits(TTree *tree = 0);
    virtual ~hits();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual void Loop();
    // virtual int gemChConverter(int ch, int vmm);
    // virtual void threePlotDrawF(TH1D* h1, TH1D* h2, TH1D* h3, TString name);
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);
};

#endif

// #ifdef hits_cxx 
hits::hits(TTree *tree) : fChain(0)
{
    
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0)
    {
        TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject(file + ending);
        if (!f || !f->IsOpen())
        {
            f = new TFile(file + ending);
        }
        f->GetObject("hits", tree);
    }
    Init(tree);
}

hits::~hits()
{
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

Int_t hits::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t hits::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (fChain->GetTreeNumber() != fCurrent)
    {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void hits::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("hits", &hits_, &b_hits_);
    fChain->SetBranchAddress("hits.id", hits_id, &b_hits_id);
    fChain->SetBranchAddress("hits.det", hits_det, &b_hits_det);
    fChain->SetBranchAddress("hits.plane", hits_plane, &b_hits_plane);
    fChain->SetBranchAddress("hits.fec", hits_fec, &b_hits_fec);
    fChain->SetBranchAddress("hits.vmm", hits_vmm, &b_hits_vmm);
    fChain->SetBranchAddress("hits.readout_time", hits_readout_time,
                             &b_hits_readout_time);
    fChain->SetBranchAddress("hits.time", hits_time, &b_hits_time);
    fChain->SetBranchAddress("hits.ch", hits_ch, &b_hits_ch);
    fChain->SetBranchAddress("hits.pos", hits_pos, &b_hits_pos);
    fChain->SetBranchAddress("hits.bcid", hits_bcid, &b_hits_bcid);
    fChain->SetBranchAddress("hits.tdc", hits_tdc, &b_hits_tdc);
    fChain->SetBranchAddress("hits.adc", hits_adc, &b_hits_adc);
    fChain->SetBranchAddress("hits.over_threshold", hits_over_threshold,
                             &b_hits_over_threshold);
    fChain->SetBranchAddress("hits.chip_time", hits_chip_time, &b_hits_chip_time);
    Notify();
}

Bool_t hits::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void hits::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain)
        return;
    fChain->Show(entry);
}
Int_t hits::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
// #endif // #ifdef hits_cxx
