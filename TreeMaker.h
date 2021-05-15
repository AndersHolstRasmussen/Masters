
#ifndef ANALYZER_EXAMPLE_TREEMAKER_H
#define ANALYZER_EXAMPLE_TREEMAKER_H


#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/util/memory>
#include <ausa/util/DynamicBranchVector.h>

#include <TH1D.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <TCut.h>

class TreeMaker : public AUSA::Sort::AbstractSortedAnalyzer{
public:
    TreeMaker(std::string fname);
    
    virtual void setup(const AUSA::Sort::SortedSetupOutput &output) override;
    virtual void analyze() override;
    virtual void saveToRootFile(AUSA::TFileWrapper& file) override;
    virtual void terminate() override;

private:
    std::string fname;

    std::unique_ptr<TFile> fOut;
    TTree *tOut;

    UInt_t MUL;
    std::unique_ptr<AUSA::DynamicBranchVector<Double_t>> bENERGY;
};


#endif //ANALYZERTEMPLATE_PADGATEDSUMMEDTELESCOPE_H
