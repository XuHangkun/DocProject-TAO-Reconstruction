#include "TAORunData.h"
#include "RadioActiveSource.h"
#include "TChain.h"
#include <string>
#include <vector>
#include "TH1F.h"
#include "TF1.h"
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "globals.h"
#include "SiPM.h"

TAORunData::TAORunData(std::string inSource,std::string filePath,int inNFile,TVector3 inPosition)
{
    source=inSource;
    calibPosition=inPosition;
    radioActiveSource=new RadioActiveSource(source,calibPosition);
    radioActiveSource->ReadInitialFitPars();

    files=filePath;
    NFiles=inNFile;

    isGamma=false;
    gammaNumber=50000;   //defaultly we use 50000 gammas.
    
}

TAORunData::TAORunData(std::string inSource,float kinetic,std::string filePath,int inNFile,TVector3 inPosition)
{
    source=inSource;
    calibPosition=inPosition;
    radioActiveSource=new RadioActiveSource(source,kinetic,calibPosition);
    radioActiveSource->ReadInitialFitPars();

    files=filePath;
    NFiles=inNFile;

    isGamma=false;
    gammaNumber=50000;   //defaultly we use 50000 gammas.
    
}

TAORunData::~TAORunData()
{
    Finalize();
}

void TAORunData::Initialize()
{
    std::cout<<"Initialize "<<std::setw(7)<<source<<" Tree";

    //create the TChain
    mainTree=new TChain("fSiPMTree");
    event_tree=new TChain("event_tree");
    fmuon=new TChain("fmuon");

    //Add All Files to the tree
    std::vector<std::string> fileNames;
    using std::stringstream;
    for(int i=0;i<NFiles;i++){
        std::string tmp_file=files;
        int tmp_pos=tmp_file.find("ProcId");
        if(tmp_pos<0) {
            fileNames.push_back(tmp_file);
            break;
        }
        stringstream tmp_ss;
        tmp_ss<<i;
        tmp_file.replace(tmp_pos,6,tmp_ss.str().c_str());
        std::cout<<tmp_file<<std::endl;
        fileNames.push_back(tmp_file);
    }

    for(int i=0;i<fileNames.size();i++){
        mainTree->Add(fileNames[i].c_str());
        event_tree->Add(fileNames[i].c_str());
        fmuon->Add(fileNames[i].c_str());
    }

    NEntries=mainTree->GetEntries();
    mainTree->SetBranchStatus("*",kFALSE);
    event_tree->SetBranchStatus("*",kFALSE);
    fmuon->SetBranchStatus("*",kFALSE);

    capTime=0;
    mainTree->SetBranchStatus("capTime",kTRUE);
    mainTree->SetBranchAddress("capTime",&capTime);

    
    totalPE=0;
    mainTree->SetBranchStatus("hit_sum",kTRUE);
    mainTree->SetBranchAddress("hit_sum",&totalPE);
    

    NCompton=0;
    mainTree->SetBranchStatus("nCompton",kTRUE);
    mainTree->SetBranchAddress("nCompton",&NCompton);
    
    Edep=0;
    event_tree->SetBranchStatus("miniJUNOGdLSEdep",kTRUE);
    event_tree->SetBranchAddress("miniJUNOGdLSEdep",&Edep);

    NParticles=0;
    mainTree->SetBranchStatus("vertex.n_particles",kTRUE);
    mainTree->SetBranchAddress("vertex.n_particles",&NParticles);


    EDepCenterX=0;
    mainTree->SetBranchStatus("GdLS.Edep_x",kTRUE);
    mainTree->SetBranchAddress("GdLS.Edep_x",&EDepCenterX);

    EDepCenterY=0;
    mainTree->SetBranchStatus("GdLS.Edep_y",kTRUE);
    mainTree->SetBranchAddress("GdLS.Edep_y",&EDepCenterY);

    EDepCenterZ=0;
    mainTree->SetBranchStatus("GdLS.Edep_z",kTRUE);
    mainTree->SetBranchAddress("GdLS.Edep_z",&EDepCenterZ);

    InitT=0;
    mainTree->SetBranchStatus("vertex.t0",kTRUE);
    mainTree->SetBranchAddress("vertex.t0",&InitT);

    InitX=0;
    mainTree->SetBranchStatus("vertex.x0",kTRUE);
    mainTree->SetBranchAddress("vertex.x0",&InitX);

    InitY=0;
    mainTree->SetBranchStatus("vertex.y0",kTRUE);
    mainTree->SetBranchAddress("vertex.y0",&InitY);

    InitZ=0;
    mainTree->SetBranchStatus("vertex.z0",kTRUE);
    mainTree->SetBranchAddress("vertex.z0",&InitZ);


    //Initialize the tree here
    float x_min=radioActiveSource->GetKinetic()*GAMMALY*0.36;
    float x_max=radioActiveSource->GetKinetic()*GAMMALY*1.2;
    int NBins=150;
    std::string histName="FullEnergyPeak_"+radioActiveSource->GetSourceLabel();
    histOfTotalPE=new TH1F(histName.c_str(),histName.c_str(),NBins,x_min,x_max);
    histOfTotalPE->GetXaxis()->SetTitle("Total PE");
    histOfTotalPE->GetYaxis()->SetTitle("Count");

    std::cout<<"-------->Done!"<<std::endl;

}

void TAORunData::GetEntry(int n)
{

    int pdgs[5]={0};
    event_tree->SetBranchStatus("vertex.pdg_code",kTRUE);
    event_tree->SetBranchAddress("vertex.pdg_code",pdgs);

    for(int i=0;i<SiPM::NSIPM;i++){
        firstHitTime[i]=0;
        hitCharge[i]=0;
    }
    const int MaxNum=200000;
    double tmp_hitTime[MaxNum]={0};
    int tmp_hitPMT[MaxNum]={0};
    mainTree->SetBranchStatus("hit_SiPM",kTRUE);
    mainTree->SetBranchAddress("hit_SiPM",tmp_hitPMT);
    mainTree->SetBranchStatus("hit_time",kTRUE);
    mainTree->SetBranchAddress("hit_time",tmp_hitTime);

    mainTree->GetEntry(n);
    event_tree->GetEntry(n);
    currentEntry=n;

    //clean vector
    //judge if the initial contains Gamma or contains positron
    isGamma=false;
    std::vector<int> tmp;
    PDGCode.swap(tmp);
    for(int i=0;i<NParticles;i++){
        PDGCode.push_back(pdgs[i]);
        if(pdgs[i]==22 || pdgs[i]==(-11)){
            isGamma=true;
        }
    }

    for(int i=0;i<totalPE;i++){
        hitCharge[tmp_hitPMT[i]]++;
        if(hitCharge[tmp_hitPMT[i]]==1){
            firstHitTime[tmp_hitPMT[i]]=tmp_hitTime[i];
        }else if(firstHitTime[tmp_hitPMT[i]]>tmp_hitTime[i]){
            firstHitTime[tmp_hitPMT[i]]=tmp_hitTime[i];
        }
    }

    event_tree->SetBranchStatus("vertex.pdg_code",kFALSE);
    mainTree->SetBranchStatus("hit_SiPM",kFALSE);
    mainTree->SetBranchStatus("hit_time",kFALSE);

}

int TAORunData::GetEntries()
{
    return mainTree->GetEntries();
}

std::vector<int>& TAORunData::GetPDGCode()
{
    return PDGCode;
}

std::vector<double>& TAORunData::GetHitTime()
{
    return hitTime;
}

void TAORunData::FillHistOfTotalPE()
{
    //We Should reset hist*/
    histOfTotalPE->Reset();
    int count=0;
    for(int i=0;i<NEntries;i++){
        
        GetEntry(i);
        if(isGamma)
        {
            count++;
            histOfTotalPE->Fill(totalPE);
        }
        if(count>=gammaNumber)
        {
            std::cout<<"Gamma Number : "<<count<<" Reached!"<<std::endl;
            break;
        }
    }
}


void TAORunData::FillHistOfTotalPE_FullEnergy()
{
    //We Should reset hist*/
    histOfTotalPE->Reset();
    int count=0;
    for(int i=0;i<NEntries;i++){
        GetEntry(i);
        if(isGamma)
        {
            count++;
            //std::cout<<Edep<<std::endl;
            if(Edep>(radioActiveSource->GetKinetic()*0.9999) && Edep<(radioActiveSource->GetKinetic()*1.0001))
            {
                histOfTotalPE->Fill(totalPE);
            }
        }
        if(count>=gammaNumber)
        {
            std::cout<<"Gamma Number : "<<count<<" Reached!"<<std::endl;
            break;
        }
    }
}

void TAORunData::FitHistOfTotalPE(std::string func)
{
    if(strcmp(func.c_str(),"MCShape") == 0)
    {
        histOfTotalPE->Fit(radioActiveSource->GetTFMCShape(),"");
    }else if(strcmp(func.c_str(),"Gaus") == 0){
        histOfTotalPE->Fit("gaus","");
    }
}

void TAORunData::AddBkg(float time)
{
    TH1F* signal=histOfTotalPE;

    //Read the file
    TFile* bkgFile=TFile::Open((ANATOP+"/input/Bkg/AssumedBkg.root").c_str());
    TH1F* bkgHist=(TH1F*)bkgFile->Get("Bkg_601.4s");

    //create bkg hist
    TH1F* sampleBkg=new TH1F("SampleBkg","SampleBkg",signal->GetNbinsX(),signal->GetXaxis()->GetXmin(),signal->GetXaxis()->GetXmax());
    float NBkg=bkgHist->GetEntries()*time/601.4;
    for(int i=0;i<NBkg;i++){
        sampleBkg->Fill(bkgHist->GetRandom());
    }
    signal->Sumw2(true);
    signal->Add(sampleBkg,1);
    signal->Sumw2(false);
    bkgFile->Close();
}

void TAORunData::SubBkg(float time)
{
    TH1F* signal=histOfTotalPE;

    //Read the file
    TFile* bkgFile=TFile::Open((ANATOP+"/input/Bkg/AssumedBkg.root").c_str());
    TH1F* bkgHist=(TH1F*)bkgFile->Get("Bkg_601.4s");

    //create bkg hist
    TH1F* sampleBkg=new TH1F("SampleBkg","SampleBkg",signal->GetNbinsX(),signal->GetXaxis()->GetXmin(),signal->GetXaxis()->GetXmax());
    float NBkg=bkgHist->GetEntries()*time/601.4;
    for(int i=0;i<NBkg;i++){
        sampleBkg->Fill(bkgHist->GetRandom());
    }
    signal->Sumw2(true);
    signal->Add(sampleBkg,-1);
    signal->Sumw2(false);
    bkgFile->Close();
}

void TAORunData::ASBkg(float time)
{
    AddBkg(time);
    SubBkg(time);
}

void TAORunData::UpdateTotalPEHist()
{
    delete histOfTotalPE;
    /*initial hist of total PE*/
    float x_min=radioActiveSource->GetKinetic()*GAMMALY*0.36;
    float x_max=radioActiveSource->GetKinetic()*GAMMALY*1.2;
    int NBins=150;
    std::string histName="FullEnergyPeak_"+radioActiveSource->GetSourceLabel();
    histOfTotalPE=new TH1F(histName.c_str(),histName.c_str(),NBins,x_min,x_max);
    histOfTotalPE->GetXaxis()->SetTitle("Total PE");
    histOfTotalPE->GetYaxis()->SetTitle("Count");
}

void TAORunData::UpdateTotalPEHist(float inxmin,float inxmax,int inNBins)
{
    delete histOfTotalPE;
    /*initial hist of total PE*/
    float x_min=inxmin;
    float x_max=inxmax;
    int NBins=inNBins;
    std::string histName="FullEnergyPeak_"+radioActiveSource->GetSourceLabel();
    histOfTotalPE=new TH1F(histName.c_str(),histName.c_str(),NBins,x_min,x_max);
    histOfTotalPE->GetXaxis()->SetTitle("Total PE");
    histOfTotalPE->GetYaxis()->SetTitle("Count");
}

void TAORunData::UpdateTotalPEHist(int inNBins,float* bins)
{
    delete histOfTotalPE;
    /*initial hist of total PE*/
    std::string histName="FullEnergyPeak_"+radioActiveSource->GetSourceLabel();
    histOfTotalPE=new TH1F(histName.c_str(),histName.c_str(),inNBins,bins);
    histOfTotalPE->GetXaxis()->SetTitle("Total PE");
    histOfTotalPE->GetYaxis()->SetTitle("Count");
}

void TAORunData::Finalize()
{

    delete histOfTotalPE;
    delete radioActiveSource;

    /* delete tree*/
    delete mainTree;
    delete event_tree;
    delete fmuon;
}

std::ostream & operator<<(std::ostream & os, const TAORunData TAORun)
{
    using namespace std;
    return os;
}


