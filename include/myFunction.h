#ifndef MYFUNCTION
#define MYFUNCTION
#include <string>

void GetReviseFactor();
void DrawRecToRealFactor();

//Reconstruct positron distribute uniformly in the detector
void RecUniPositron(float Kinetic,int NFile=20);
void PlotPositronRecPerformence();


void CalculateLinearFactor(const char* source);

#endif
