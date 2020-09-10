#ifndef VERTEXCONSTRUCTION
#define VERTEXCONSTRUCTION
#include "SiPM.h"
#include "TAORunData.h"
#include "TVector3.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include "TGraph2D.h"
#include "TGraphErrors.h"

class VertexReconstruction
{
    public:
        VertexReconstruction(std::string inSubFile,int inNFile,std::string outFile);
        ~VertexReconstruction();

        //Initialize
        void Initialize();
        //Fill grFactor
        float RecToRealFactor(float RecTheta,float RecR);
        //Finalize
        void Finalize();
        
        //Calculate charge center of a event
        bool ChargeCenter(int eventId);

        //Reconstruction
        void Reconstruction();

        //Get and Set
        TGraph2D* GetGrFactor() { return grFactor; }

    private:
        SiPM* sipm;
        TAORunData* taoRunData;
        TGraph2D* grFactor;

        //N'th reconstruction information
        int currentEventID;
        float recVertexX;
        float recVertexY;
        float recVertexZ;
        float recVertexR;
        float recVertexTheta;
        float initX;
        float initY;
        float initZ;
        float cpuTime;
        float EDepCenterX;
        float EDepCenterY;
        float EDepCenterZ;
        float EDepCenterR;
        float EDepCenterTheta;

        float RBias;
        float ThetaBias;

        //Files
        std::string outputFileName;
        TFile* outputFile;
        TTree* resTree;
        
        //calibration res.
        TFile* calibFactorFile;
        TGraphErrors* calibFactor;
    
};


#endif
