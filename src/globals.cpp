#include "globals.h"
#include "TStyle.h"
#include <string>


void  SetGlobalDrawOption()
{
    gStyle->SetCanvasDefH(560*0.8);
    gStyle->SetCanvasDefY(560*0.8);
    gStyle->SetCanvasDefW(700*0.8);
    gStyle->SetCanvasDefX(700*0.8);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.15);

    gStyle->SetStatX(0.52);
    gStyle->SetStatY(1-gStyle->GetPadTopMargin());


    gStyle->SetLabelFont(132,"X");
    gStyle->SetLabelFont(132,"Y");
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");

    gStyle->SetLegendFont(132);
    gStyle->SetLegendTextSize(0.04);

    gStyle->SetTitleFont(132);
    gStyle->SetTitleSize(0.08);
    gStyle->SetTitleOffset(0.8);

    gStyle->SetTitleXSize(0.06);
    gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleYSize(0.06);
    gStyle->SetTitleYOffset(1.2);

    gStyle->SetOptFit(kFALSE);
    gStyle->SetOptStat(kFALSE);
}
