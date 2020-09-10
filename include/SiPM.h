#ifndef SIPM
#define SIPM

#include "TVector3.h"
#include <vector>
#include <iostream>
#include "TF1.h"
#include "TF2.h"

class SiPM
{
public:
    SiPM();
    ~SiPM();

    static const int NSIPM=4074;      //number of PMT
    //Read SiPM infomation
    void Initialize();
    //Calculate Solid Angle of SiPM
    void CalSiPMSolidAngle();
    void CalGapCorrectFactor();

    //Get and Set function
    float GetSiPMPositionX(int id) { return position[id].X(); }
    float GetSiPMPositionY(int id) { return position[id].Y(); }
    float GetSiPMPositionZ(int id) { return position[id].Z(); }
    float GetGapCorrectFactor(int id) { return gapCorrectFactor[id]; }
    float GetSolidAngle() { return SiPMSolidAngle; }
    float GetCoverRatio() { return CoverRatio; }

    float GetAccuGapCorrFac(TVector3 refer,int id);
    double MicroSiPMSolidAngle(double* x,double* par);
    double MicroSurfSolidAngle(double* x,double* par);

    //Reload operator
    friend std::ostream& operator<<(std::ostream& os,const SiPM sipm);

private:
    float SiPMX;                      //Size of SiPM
    float SiPMY;
    float SiPMPosR;                   //Radius of SiPM center
    float CoverRatio;
    float gapCorrectFactor[NSIPM];    //factor to correct the effect caused by gad becween SiPM
    int ID[NSIPM];                    //SiPM id
    TVector3 position[NSIPM];         //position 
    float darkNoise[NSIPM];           //Dark noise

    float SiPMSolidAngle;             //Solid angle of a SiPM, relative to center
    float eTheta[NSIPM];
    float nTheta[NSIPM];
    float ePhi[NSIPM];
    float nPhi[NSIPM];
    TF2* fMicroSiPMSolidAngle;
    TF2* fMicroSurfSolidAngle;
};
#endif
