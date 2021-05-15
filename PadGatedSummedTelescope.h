//
// Created by munk on 17-05-15.
//

#ifndef ANALYZERTEMPLATE_PADGATEDSUMMEDTELESCOPE_H
#define ANALYZERTEMPLATE_PADGATEDSUMMEDTELESCOPE_H


#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <TH1D.h>

class PadGatedSummedTelescope : public AUSA::Sort::AbstractSortedAnalyzer{
public:
    PadGatedSummedTelescope(std::string padName, std::string dssdName);
    virtual void setup(const AUSA::Sort::SortedSetupOutput &output) override;

    virtual void analyze() override;

private:
    std::string padName, dssdName;
    AUSA::Sort::SortedSingleOutput padOut;
    AUSA::Sort::DoubleOutput dOut;

    TH1D hist;
};


#endif //ANALYZERTEMPLATE_PADGATEDSUMMEDTELESCOPE_H
