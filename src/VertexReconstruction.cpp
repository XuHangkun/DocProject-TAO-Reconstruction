#include "VertexReconstruction.h"
#include "TTree.h"
#include <string.h>
#include "SiPM.h"
#include "TAORunData.h"
#include "TMath.h"
#include <iostream>
#include "TGraph2D.h"
#include "TFile.h"

VertexReconstruction::VertexReconstruction(std::string inFile,int inNFile,std::string outFile)
{
    currentEventID=0;
    recVertexX=0;
    recVertexY=0;
    recVertexZ=0;
    recVertexR=0;
    initX=0;
    initY=0;
    initZ=0;
    cpuTime=0;
    EDepCenterX=0;
    EDepCenterY=0;
    EDepCenterZ=0;
    EDepCenterR=0;
    RBias=0;
    sipm=new SiPM();
    taoRunData=new TAORunData("Positron",inFile,inNFile);
    outputFileName=outFile;
    grFactor=new TGraph2D();
    Initialize();
    //calibFactorFile=TFile::Open((ANATOP+"/input/Calib/ChargeCenterCorrFactor.root").c_str());
    //calibFactor=(TGraphErrors*)calibFactorFile->Get("Graph");
}

VertexReconstruction::~VertexReconstruction()
{
}

void VertexReconstruction::Initialize()
{
    taoRunData->Initialize();
    outputFile=new TFile(outputFileName.c_str(),"recreate");
    resTree=new TTree("data","Vertex Reconstruction");

    //charge center result
    resTree->Branch("recVertexX",&recVertexX,"recVertexX/F");  
    resTree->Branch("recVertexY",&recVertexY,"recVertexY/F");
    resTree->Branch("recVertexZ",&recVertexZ,"recVertexZ/F");
    resTree->Branch("recVertexR",&recVertexR,"recVertexR/F");
    resTree->Branch("recVertexTheta",&recVertexTheta,"recVertexTheta/F");
    
    //event info
    resTree->Branch("initX",&initX,"initX/F");
    resTree->Branch("initY",&initY,"initY/F");
    resTree->Branch("initZ",&initZ,"initZ/F");
    resTree->Branch("EDepCenterX",&EDepCenterX,"EDepCenterX/F");
    resTree->Branch("EDepCenterY",&EDepCenterY,"EDepCenterY/F");
    resTree->Branch("EDepCenterZ",&EDepCenterZ,"EDepCenterZ/F");
    resTree->Branch("EDepCenterR",&EDepCenterR,"EDepCenterR/F");
    resTree->Branch("EDepCenterTheta",&EDepCenterTheta,"EDepCenterTheta/F");

    //reconstruction comparation
    resTree->Branch("RBias",&RBias,"RBias/F");
    resTree->Branch("ThetaBias",&ThetaBias,"ThetaBias/F");
    

    //calculate info
    resTree->Branch("cpuTime",&cpuTime,"cpuTime/F");


       
}

void VertexReconstruction::Finalize()
{
    resTree->Write();
    outputFile->Close();
    //delete calibFactor;
    //calibFactorFile->Close();
}

bool VertexReconstruction::ChargeCenter(int event_id)
{
    using namespace std;
    taoRunData->GetEntry(event_id);
    if(!taoRunData->GetIsGamma()) {return false; }

    initX=taoRunData->GetInitX();
    initY=taoRunData->GetInitY();
    initZ=taoRunData->GetInitZ();
    cpuTime=0;
    EDepCenterX=taoRunData->GetEdepCenterX();
    EDepCenterY=taoRunData->GetEdepCenterY();
    EDepCenterZ=taoRunData->GetEdepCenterZ();
    EDepCenterR=taoRunData->GetEdepCenterR();
    EDepCenterTheta=taoRunData->GetEdepCenterTheta()*180/TMath::Pi();
    
    float xx=0;
    float yy=0;
    float zz=0;
    float factor=0;
    for(int i=0;i<SiPM::NSIPM;i++){
       //cout<<sipm->GetSiPMPositionX(i)<<endl;
       xx+=sipm->GetSiPMPositionX(i)*sipm->GetGapCorrectFactor(i)*taoRunData->GetHitCharge(i); 
       yy+=sipm->GetSiPMPositionY(i)*sipm->GetGapCorrectFactor(i)*taoRunData->GetHitCharge(i); 
       zz+=sipm->GetSiPMPositionZ(i)*sipm->GetGapCorrectFactor(i)*taoRunData->GetHitCharge(i); 
       factor+=sipm->GetGapCorrectFactor(i)*taoRunData->GetHitCharge(i);

       //xx+=sipm->GetSiPMPositionX(i)*taoRunData->GetHitCharge(i); 
       //yy+=sipm->GetSiPMPositionY(i)*taoRunData->GetHitCharge(i); 
       //zz+=sipm->GetSiPMPositionZ(i)*taoRunData->GetHitCharge(i); 
       //factor+=taoRunData->GetHitCharge(i);
    }
    xx/=factor;
    yy/=factor;
    zz/=factor;
    float rr=sqrt(xx*xx+yy*yy+zz*zz);
    //factor*=RecToRealFactor(acos(xx/rr)*180/TMath::Pi(),rr);
    recVertexX=(xx-2.06)/(0.705);
    recVertexY=(yy-2.06)/(0.705);
    recVertexZ=(zz-2.06)/(0.705);
    /*
    if( acos(zz/rr)<0.123655 && acos(zz/rr)>3.01794 ){
        recVertexX=xx/(calibFactor->Eval(rr)*factor);
        recVertexY=yy/(calibFactor->Eval(rr)*factor);
        recVertexZ=yy/(calibFactor->Eval(rr)*factor);
    }*/
    recVertexR=sqrt(recVertexX*recVertexX+recVertexY*recVertexY+recVertexZ*recVertexZ);
    recVertexTheta=acos(recVertexZ/recVertexR)*180/TMath::Pi();
    RBias=EDepCenterR-recVertexR;
    ThetaBias=EDepCenterTheta-recVertexTheta;
    
    //cout<<"Init :"<<EDepCenterX<<"\t"<<EDepCenterY<<"\t"<<EDepCenterZ<<endl; 
    //cout<<"Rec  :"<<recVertexX<<"\t"<<recVertexY<<"\t"<<recVertexZ<<endl; 
   
    return true;

}


void VertexReconstruction::Reconstruction()
{
    int N=taoRunData->GetEntries();
    int successN=0;
    for(int i=0;i<N;i++){
        bool situ=ChargeCenter(i);
        if(!situ) continue;
        resTree->Fill();

        successN++;
        if((successN)%10000==0){
            std::cout<<(successN)<<" events have been reconstructed!"<<std::endl;
        }
        if(successN>=50000) { break; }
    }
}

float VertexReconstruction::RecToRealFactor(float RecR,float RecTheta)
{
    const int calibThetaNum=4;
    float calibRealTheta[calibThetaNum]={0,30,60,90};
    const int calibRNum=5;
    float calibRealR[calibRNum]={0,250,450,650,850};  //mm

    float calibRecX[calibThetaNum][calibRNum]={ 0.6,    -0.3,   -0.8,   -0.8,   -0.1,
                                                0.6,    88.2,   157.9,  227.2,  298.2,
                                                0.6,    152.3,  274.0,  392.9,  516.0,
                                                0.6,    176.4,  315.3,  454.2,  597.6,
                                                };
    float calibRecXError[calibThetaNum][calibRNum]={0};

    float calibRecY[calibThetaNum][calibRNum]={0};
    float calibRecYError[calibThetaNum][calibRNum]={0};

    float calibRecZ[calibThetaNum][calibRNum]={-0.1,    175.1,  315.7,  452.7,  556.3,
                                               -0.1,    153.0,  274.1,  392.7,  516.9,
                                               -0.1,    87.2,   158.2,  226.7,  297.7,
                                               -0.1,    -0.5,   -0.4,   0.4,    0.4};
    float calibRecZError[calibThetaNum][calibRNum]={0};
    
    float calibRecR[calibThetaNum][calibRNum]={0};
    float factor[calibThetaNum][calibRNum]={1,1,1,1,1,
                                            1,1,1,1,1,
                                            1,1,1,1,1,
                                            1,1,1,1,1};
    for(int i=0;i<calibThetaNum;i++){
        for(int j=1;j<calibRNum;j++){
            calibRecR[i][j]=sqrt(calibRecX[i][j]*calibRecX[i][j]+calibRecY[i][j]*calibRecY[i][j]+calibRecZ[i][j]*calibRecZ[i][j]);
            factor[i][j]=calibRecR[i][j]/calibRealR[j];
        }
        factor[i][0]=factor[i][1];
    }

    TGraph gr_factors[calibThetaNum];
    for(int i=0;i<calibThetaNum;i++){
        gr_factors[i]=TGraph(calibRNum,calibRealTheta,factor[i]);
    }

    //We use linear inter...
    int thetaIndex=0;
    for(int i=0;i<calibThetaNum-1;i++){
        if(calibRealTheta[i]<=RecTheta && calibRealTheta[i+1]>=RecTheta){
            thetaIndex=i;
            break;
        }
    }
    
    return gr_factors[thetaIndex].Eval(RecR);

}
