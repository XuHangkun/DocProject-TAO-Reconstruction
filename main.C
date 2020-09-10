#include <iostream>
#include "SiPM.h"
#include "VertexReconstruction.h"
#include "globals.h"
#include "myFunction.h"
#include <string>
int main(int argv,char** argc)
{

    using namespace std;
    //DrawRecToRealFactor();
    RecUniPositron(atof(argc[1]),atoi(argc[2]));
    //CalculateLinearFactor(argc[1]);
    //PlotPositronRecPerformence();
}
