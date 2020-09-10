#include "myFunction.h"
#include "TString.h"
#include <string>
#include "globals.h"
#include "VertexReconstruction.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"


void RecUniPositron(float Kinetic,int NFile)
{
        std::string filePattern=DATATOP+Form("/reconstruction/Positron/Positron_%.0fMeV_GdLS_VProcId.root",Kinetic);
        std::string output=ANATOP+Form("/result/PositronStg2/Positron_%.0fMeV_GdLS.root",Kinetic);
        VertexReconstruction* vR=new VertexReconstruction(filePattern,NFile,output);
        vR->Reconstruction();
        vR->Finalize();
}

void PlotPositronRecPerformence()
{
    const int PositronN=9;
    float PositronKinetic[PositronN]={0,1,2,3,4,5,6,7,8};
    TH2F* hist2d[PositronN];
    TGraph* gr[PositronN];
    int color[PositronN]={1,2,3,4,5,6,7,8,14};
    int marker[PositronN]={20,21,22,23,24,25,26,27,32};
    SetGlobalDrawOption();
    TCanvas* c1=new TCanvas("Test","Test");
    TString canvasName=ANATOP+Form("/result/PositronRecPerformance.pdf");
    c1->Print(canvasName+"[");
    c1->SetGrid();
    TLegend* leg=new TLegend(0.15,0.3,0.45,0.85);
    TBox* box=new TBox(0,-150,pow(650,3),300);
    box->SetFillStyle(3003);
    box->SetFillColor(kMagenta);
    for(int i=0;i<PositronN;i++){
        std::string histName=Form("e^{+} kinetic %.1fMeV",PositronKinetic[i]);
        hist2d[i]=new TH2F(histName.c_str(),histName.c_str(),100,0,pow(900,3),100,-150,300);
        hist2d[i]->GetXaxis()->SetTitle("R_{edep}^{3}[mm^{3}]");
        hist2d[i]->GetYaxis()->SetTitle("R_{rec} - R_{edep}[mm]");
        gr[i]=new TGraph();
        gr[i]->GetXaxis()->SetTitle("R_{edep}^{3}[mm^{3}]");
        gr[i]->GetYaxis()->SetTitle("Mean(R_{rec} - R_{edep})[mm]");
        gr[i]->SetLineColor(color[i]);
        gr[i]->SetMarkerColor(color[i]);
        gr[i]->SetMarkerStyle(marker[i]);
        gr[i]->GetYaxis()->SetRangeUser(-100,200);
        gr[i]->GetXaxis()->SetRangeUser(0,pow(900,3));
        leg->AddEntry(gr[i],Form("e^{+} kinetic %.1fMeV",PositronKinetic[i]),"P");
        std::string fileName=ANATOP+Form("/result/Positron/Positron_%.0fMeV_GdLS.root",PositronKinetic[i]);
        std::cout<<fileName<<std::endl;
        TFile* tmp_file=TFile::Open(fileName.c_str());
        TTree* tmp_tree=(TTree*)tmp_file->Get("data");
        tmp_tree->SetBranchStatus("*",kFALSE);
        float RBias=0;
        tmp_tree->SetBranchStatus("RBias",kTRUE);
        tmp_tree->SetBranchAddress("RBias",&RBias);
        float EDepCenterR=0;
        tmp_tree->SetBranchStatus("EDepCenterR",kTRUE);
        tmp_tree->SetBranchAddress("EDepCenterR",&EDepCenterR);
        for(int j=0;j<tmp_tree->GetEntries();j++){
            tmp_tree->GetEntry(j);
            //std::cout<<RBias<<"\t"<<EDepCenterR<<std::endl;
            hist2d[i]->Fill(pow(EDepCenterR,3),-1.0*RBias);
        }
        hist2d[i]->Draw("colz");
        box->Draw("same");
        for(int j=0;j<hist2d[i]->GetNbinsY();j++){
            TH1D* tmp_hist=hist2d[i]->ProjectionY("py",j,j+1);
            gr[i]->SetPoint(j,hist2d[i]->GetXaxis()->GetBinCenter(j),tmp_hist->GetMean());
        }
        delete tmp_tree;
        tmp_file->Close();
        delete tmp_file;
        c1->Print(canvasName);
    }

    for(int i=0;i<PositronN;i++){
        if(i==0) gr[i]->Draw("aP");
        else gr[i]->Draw("p same");
    }
    box->SetY1(0);
    box->SetY2(100);
    box->Draw("same");
    leg->Draw();
    c1->Print(canvasName);
    c1->Print(canvasName+"]");
    
}

void GetReviseFactor()
{
    const int calibPointNum=17;
    int calibX[calibPointNum]={0, 0, 0, 0, 0 ,125 ,225 ,325 ,425 ,217 ,390, 563 ,736 , 250 ,450 ,650 ,850};
    int calibZ[calibPointNum]={0, 250, 450, 650, 850, 217, 390, 563, 736, 125, 225, 325, 425, 0, 0, 0, 0};
    for(int i=0;i<calibPointNum;i++){
        std::string filePattern=DATATOP+Form("/reconstruction/Calib/Positron/Positron_0MeV_X%d_Y0_Z%d_VProcId.root",calibX[i],calibZ[i]);
        std::string output=ANATOP+Form("/result/Rec_Positron_0MeV_X%d_Y0_Z%d_VProcId.root",calibX[i],calibZ[i]);
        VertexReconstruction* vR=new VertexReconstruction(filePattern,5,output);
        vR->Reconstruction();
        vR->Finalize();
    }
}

#include "TGraph2D.h"
void DrawRecToRealFactor()
{

    SetGlobalDrawOption();
    TCanvas* c1=new TCanvas("Test","Test");
    TString canvasName=ANATOP+Form("/result/RecToRealFactor.pdf");
    c1->Print(canvasName+"[");
    c1->SetGrid();
    std::string filePattern=DATATOP+Form("/reconstruction/Calib/Positron_0MeV_X%d_Y0_Z%d_VProcId.root",0,0);
    std::string output=ANATOP+Form("/result/Rec_Positron_0MeV_X%d_Y0_Z%d_VProcId.root",0,0);
    VertexReconstruction* vR=new VertexReconstruction(filePattern,5,output);
    TGraph2D* tmp=new TGraph2D();
    for(int i=0;i<90;i++){
        for(int j=0;j<900;j++){
            tmp->SetPoint(i*900+j,i,j,vR->RecToRealFactor(i,j));
        }
    }
    tmp->Draw("colz");
    c1->Print(canvasName);
    c1->Print(canvasName+"]");

}

void CalculateLinearFactor(const char* source)
{
    const int CalibPointN=25;
    float calibPosZ[CalibPointN]={-850,-800,-750,-700,-650,-600,-550,-500,-400,-300,-200,-100,0,
                                100, 200, 300, 400, 500, 550, 600, 650, 700, 750, 800, 850};
    float recPosZ[CalibPointN]={0};
    float recPosZError[CalibPointN]={0};
    float calibPosZError[CalibPointN]={0};
    SetGlobalDrawOption();
    TCanvas* c1=new TCanvas("Test","Test");
    TString canvasName=ANATOP+Form("/result/%sRecToRealFactor.pdf",source);
    c1->Print(canvasName+"[");
    c1->SetGrid();
    gStyle->SetOptFit(true);
    for(int i=0;i<CalibPointN;i++){
        //reconstruction
        /*
        std::string filePattern=DATATOP+Form("/reconstruction/Calib/%s/%s_%dmm_vProcId.root",source,source,calibPosZ[i]);
        std::string output=ANATOP+Form("/result/Calib/Rec_%s_0.662MeV_X%d_Y0_Z%d.root",source,0,calibPosZ[i]);
        VertexReconstruction* vR=new VertexReconstruction(filePattern,30,output);
        vR->Reconstruction();
        vR->Finalize();
        */
        std::string fileName=ANATOP+Form("/result/Calib/Rec_%s_0.662MeV_X%.0f_Y0_Z%.0f.root",source,0.,calibPosZ[i]);
        std::cout<<fileName<<std::endl;
        TFile* tmp_file=TFile::Open(fileName.c_str());
        TTree* tmp_tree=(TTree*)tmp_file->Get("data");
        tmp_tree->SetBranchStatus("*",kFALSE);
        float recZ=0;
        tmp_tree->SetBranchStatus("recVertexZ",kTRUE);
        tmp_tree->SetBranchAddress("recVertexZ",&recZ);
        std::string histName=Form("%s Calib. Z %.0fmm",source,calibPosZ[i]);
        TH1F* hist=new TH1F(histName.c_str(),histName.c_str(),200,(0.7*calibPosZ[i]-200),(0.7*calibPosZ[i]+200));
        hist->GetXaxis()->SetTitle("Rec Vertex Z[mm]");
        hist->GetYaxis()->SetTitle("Count");
        for(int j=0;j<tmp_tree->GetEntries();j++){
            tmp_tree->GetEntry(j);
            hist->Fill(recZ);
        }
        float maxValue=hist->GetBinCenter(hist->GetMaximumBin());
        hist->Fit("gaus","","",maxValue-40,maxValue+50);
        recPosZ[i]=hist->GetFunction("gaus")->GetParameter(1);
        recPosZError[i]=hist->GetFunction("gaus")->GetParError(1);
        hist->Draw();
        c1->Print(canvasName);
        delete hist;
        delete tmp_tree;
        tmp_file->Close();
        delete tmp_file;
    }
    TGraphErrors* gr_factor=new TGraphErrors(CalibPointN,calibPosZ,recPosZ,calibPosZError,recPosZError);
    //TGraphErrors* gr_factor=new TGraphErrors(CalibPointN,recPosZ,calibPosZ,recPosZError,calibPosZError);
    gr_factor->SetTitle(Form("%s Calibration",source));
    gr_factor->GetXaxis()->SetTitle("Calib Vertex Z[mm]");
    gr_factor->GetYaxis()->SetTitle("Rec Vertex Z[mm]");
    gr_factor->Fit("1++x","","",-660,660);
    gr_factor->Draw("AP");
    gr_factor->SaveAs((ANATOP+"/input/ChargeCenterCorrFactor.root").c_str());
    c1->Print(canvasName);
    c1->Print(canvasName+"]");
    
    
}
