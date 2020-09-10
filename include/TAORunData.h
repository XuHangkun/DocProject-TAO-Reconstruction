/*****************************************************************************
*  OpenST Basic tool library                                                 *
*  Copyright (C) 2020 Xu Hangkun  xuhangkun@ihep.ac.cn                       *
*                                                                            *
*  This file is part of OST.                                                 *
*                                                                            *
*  This program is free software; you can redistribute it and/or modify      *
*  it under the terms of the GNU General Public License version 3 as         *
*  published by the Free Software Foundation.                                *
*                                                                            *
*  You should have received a copy of the GNU General Public License         *
*  along with OST. If not, see <http://www.gnu.org/licenses/>.               *
*                                                                            *
*  Unless required by applicable law or agreed to in writing, software       *
*  distributed under the License is distributed on an "AS IS" BASIS,         *
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  *
*  See the License for the specific language governing permissions and       *
*  limitations under the License.                                            *
*                                                                            *
*  @file     RunData.h                                                       *
*  @brief    ABC of physics run data                                         *
*  Details.                                                                  *
*                                                                            *
*  @author   Xu Hangkun                                                      *
*  @email    xuhangkun@ihep.ac.cn                                            *
*  @version  1.0.0.1(版本号)                                                 *
*  @date     2020-08-24                                                      *
*  @license  GNU General Public License (GPL)                                *
*                                                                            *
*----------------------------------------------------------------------------*
*  Remark         : Description                                              *
*----------------------------------------------------------------------------*
*  Change History :                                                          *
*  <Date>     | <Version> | <Author>       | <Description>                   *
*----------------------------------------------------------------------------*
*  2020/08/24 | 1.0.0.1   | Xu Hangkun     | Create file                     *
*----------------------------------------------------------------------------*
*                                                                            *
*****************************************************************************/

#ifndef TAORUNDATA
#define TAORUNDATA

#include "RunData.h"
#include <string>
#include <vector>
#include "TH1F.h"
#include "RadioActiveSource.h"
#include <iostream>
#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TChain.h"
#include "SiPM.h"
#include <cmath>


/**
 * @brief: run data model of TAO
 */
class TAORunData : public RunData
{
private:
    //General information
    TChain* mainTree;         //the tree which include all variables we need
    TChain* event_tree;
    TChain* fmuon;
    int NEntries;

    std::string source;        //source name
    TVector3 calibPosition;
    RadioActiveSource* radioActiveSource;
    std::string files;
    int NFiles;

    //Physical data
    int currentEntry;
    float capTime;              //neutron capture time
    int totalPE;                //total PE of a event
    int NCompton;               //Numbers of compton scattering of a gamma
    float Edep;                 //energy deposit in GdLS
    int NParticles;
    bool isGamma;
    std::vector<int> PDGCode;   //pdg code of primary particles
    std::vector<double> hitTime;//hit time of photon on PMT
    double firstHitTime[SiPM::NSIPM]={0};//first hit time of photon on PMT
    double hitCharge[SiPM::NSIPM]={0};   //hit 
    float EDepCenterX;          //X value of energy deposit position
    float EDepCenterY;          //Y value of energy deposit position
    float EDepCenterZ;          //Z value of energy deposit position
    float InitT;          //X value of energy deposit position
    float InitX;          //X value of energy deposit position
    float InitY;          //Y value of energy deposit position
    float InitZ;          //Z value of energy deposit position
    
    int gammaNumber;            //The number of gamma used in calibration

    TH1F* histOfTotalPE;        //hist Of Total PE

public:
    //filePath is the full path of TAORunData, but it can contain ProcId to include many root file(inNFile)
    TAORunData(std::string inSource,std::string filePath,int inNFile=1,TVector3 inPosition=TVector3(0,0,0));
    TAORunData(std::string inSource,float kinetic,std::string filePath,int inNFile=1,TVector3 inPosition=TVector3(0,0,0));
    ~TAORunData();

    /*Open the file and initialize the tree*/
    virtual void Initialize();

    /*Close the file*/
    virtual void Finalize();

    /*read n'th event in mainTree*/
    virtual void GetEntry(int n);

    /*Get the number of events in the tree*/
    virtual int GetEntries();
    

    /*get function*/
    inline float GetCapTime();
    inline int GetCurrentEntry();
    inline int GetTotalPE();
    bool GetIsGamma() {return isGamma; }
    inline int GetNCompton();
    inline float GetEdep();
    inline int GetNParticles();
    std::vector<int>& GetPDGCode();
    std::vector<double>& GetHitTime();
    double GetFirstHitTime(int pmt_id) { return firstHitTime[pmt_id]; }
    double GetHitCharge(int pmt_id) { return hitCharge[pmt_id]; }
    inline float GetEdepCenterX();
    inline float GetEdepCenterY();
    inline float GetEdepCenterZ();
    float GetEdepCenterR() { return sqrt(EDepCenterX*EDepCenterX+EDepCenterY*EDepCenterY+EDepCenterZ*EDepCenterZ); }
    float GetEdepCenterTheta() { return acos(EDepCenterZ/(GetEdepCenterR())); }
    float GetInitT() { return InitT; }
    float GetInitX() { return InitX; } 
    float GetInitY() { return InitY; } 
    float GetInitZ() { return InitZ; } 
    float GetInitR() { return sqrt(InitX*InitX+InitY*InitY+InitZ*InitZ); }
    float GetInitTheta() { return acos(InitZ/(GetEdepCenterR())); }
    inline int   GetGammaNumber();
    inline RadioActiveSource* GetRadioActiveSource();
    TH1F* GetHistOfTotalPE() { return histOfTotalPE; }
    TF1* GetTFMCShape() { return radioActiveSource->GetTFMCShape(); }

    /* Set function */
    inline void   SetGammaNumber(int inGammaN);

    /*get histgrom of total PE*/
    void FillHistOfTotalPE();

    /*get total PE histogram of full energy peak*/
    void FillHistOfTotalPE_FullEnergy();

    /*Fit Hist of Total PE Spectrum*/
    void FitHistOfTotalPE(std::string func="MCShape");

    /*Add and substract Bkg Effect*/
    void AddBkg( float time=500);  //Add
    void SubBkg(float time=500);  //Substract
    void ASBkg( float time=500);   //Add and Sub

    /*Update hist*/
    void UpdateTotalPEHist();
    void UpdateTotalPEHist(float xmin,float xmax,int NBins);
    void UpdateTotalPEHist(int NBins,float* bins);

    //Reload operator <<
    friend std::ostream & operator<<(std::ostream & os, const TAORunData TAORun);
        
};

inline void TAORunData::SetGammaNumber(int inGammaN)
{
    gammaNumber=inGammaN;
}


inline int TAORunData::GetCurrentEntry()
{
    return currentEntry;
}

inline float TAORunData::GetCapTime()
{
    return capTime;
}

inline int TAORunData::GetTotalPE()
{
    return totalPE;
}

inline int TAORunData::GetNCompton()
{
    return NCompton;
}

inline float TAORunData::GetEdep()
{
    return Edep;
}

inline int TAORunData::GetNParticles()
{
    return NParticles;
}

inline float TAORunData::GetEdepCenterX()
{
    return EDepCenterX;
}

inline float TAORunData::GetEdepCenterY()
{
    return EDepCenterY;
}

inline float TAORunData::GetEdepCenterZ()
{
    return EDepCenterZ;
}

inline int TAORunData::GetGammaNumber()
{
    return gammaNumber;
}

inline RadioActiveSource* TAORunData::GetRadioActiveSource()
{
    return radioActiveSource;
}
#endif
