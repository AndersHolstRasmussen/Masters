#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/geometry/Box.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/setup/SingleSidedSiliconDetector.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Ion.h>
#include <ausa/eloss/Default.h>
#include <ausa/util/memory>
#include <ausa/constants/Mass.h>
#include <ausa/output/OutputConvenience.h>
#include <Math/Vector3D.h>
#include "Hit.h"
#include <TROOT.h>
#include <ctime>

using namespace std;
using namespace AUSA;
using namespace ROOT::Math;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace AUSA::Constants;

class Li8Analysis : public AbstractSortedAnalyzer {
public:
    Li8Analysis(Target &target, TFile *output) : target(target) {

        IonRange = 166 * 1e-6;
        targetThickness = 226 * 1e-6;
        ALPHA = new ParticleType("He4");
        Be8 = new ParticleType("Be8");
        numberOfAlphas = 0;
        numberOfBetas = 0;

        tree = new TTree("tree", "analyzed tree");
        tree->Branch("mul", &mul);
        tree->Branch("mulAlpha", &mulAlpha);
        tree->Branch("mulBeta", &mulBeta);
        tree->Branch("num", &NUM);
        tree->Branch("clock", &CLOCK);
        // tree->Branch("cosang", &cosang);

        
        v_p0 = make_unique<DynamicBranchVector<TVector3>>(*tree, "p0");
        v_p1 = make_unique<DynamicBranchVector<TVector3>>(*tree, "p1");
        v_dir0 = make_unique<DynamicBranchVector<TVector3>>(*tree, "dir0");
        v_pos0 = make_unique<DynamicBranchVector<TVector3>>(*tree, "pos0");
        v_dir1 = make_unique<DynamicBranchVector<TVector3>>(*tree, "dir1");
        v_pos1 = make_unique<DynamicBranchVector<TVector3>>(*tree, "pos1");
        v_Ex = make_unique<DynamicBranchVector<double>>(*tree, "Ex", "1");

        v_Edep = make_unique<DynamicBranchVector<double>>(*tree, "Edep", "mulAlpha");
        v_E = make_unique<DynamicBranchVector<double>>(*tree, "E", "mulAlpha");
        v_angToBeam = make_unique<DynamicBranchVector<double>>(*tree, "angToBeam", "mulAlpha");
        v_ang = make_unique<DynamicBranchVector<double>>(*tree, "angle", "mulAlpha");
        cosA = make_unique<DynamicBranchVector<double>>(*tree, "cosA", "1");
        v_cosAngAll = make_unique<DynamicBranchVector<double>>(*tree, "cosAngAll", "mul");
        v_dE = make_unique<DynamicBranchVector<double>>(*tree, "dE", "mulAlpha");
        v_F = make_unique<DynamicBranchVector<short>>(*tree, "FI", "mulAlpha");
        v_B = make_unique<DynamicBranchVector<short>>(*tree, "BI", "mulAlpha");
        v_i = make_unique<DynamicBranchVector<short>>(*tree, "i", "mulAlpha");
        v_ibeta = make_unique<DynamicBranchVector<int>>(*tree, "ibeta", "mulBeta");
        v_Ebeta = make_unique<DynamicBranchVector<double>>(*tree, "Ebeta", "mulBeta");
        v_EPbeta = make_unique<DynamicBranchVector<double>>(*tree, "EPbeta", "mulBeta");
        v_betaAlpaAngle0 = make_unique<DynamicBranchVector<double>>(*tree, "betaAlphaAngle0", "1");
        v_betaAlpaAngle1 = make_unique<DynamicBranchVector<double>>(*tree, "betaAlphaAngle1", "1");
        v_bTheta = make_unique<DynamicBranchVector<double>>(*tree, "bTheta", "mulBeta");
        v_bPhi = make_unique<DynamicBranchVector<double>>(*tree, "bPhi", "mulBeta");
        v_ptot = make_unique<DynamicBranchVector<double>>(*tree, "ptot", "1");
        v_Etot = make_unique<DynamicBranchVector<double>>(*tree, "Etot", "1");
        SiCalc = defaultRangeInverter("He4", "Silicon");
        for (auto &layer : target.getLayers()) {
            targetCalcs.push_back(defaultRangeInverter(Ion::predefined("He4"), layer.getMaterial()));
        }


    }

    void setup(const SortedSetupOutput &output) override {
        AbstractSortedAnalyzer::setup(output);
        for(size_t i = 0; i < output.dssdCount(); i++){
            auto dlF = getFrontDeadLayer(output.getDssdOutput(i).detector());
            auto dlB = getBackDeadLayer(output.getDssdOutput(i).detector());
            deadlayerF.push_back(dlF);
            deadlayerB.push_back(dlB);
        }
    }

    void analyze() override {
        clear();
        findHits();
        determineTypes();  
        if (numberOfAlphas < 2) return;
        doAnalysis();
        if (moment > 40000) return;
        if (mulBeta < 1) return;
        if(cosang > -0.95) return;
        CLOCK = output.getScalerOutput("CLOCK").getValue();
        tree->Fill();
        NUM++;
        

    }

    void findHits() {
        for(int i = 0; i < output.dssdCount(); i++){
            auto &out = output.getDssdOutput(i);
            auto &pad = output.getSingleOutput(i);
            auto &det = out.detector();
            auto m = AUSA::mul(out);


            for(UInt_t j = 0; j < m; j++){
                Hit hit;
                auto dE = fEnergy(out, j) - bEnergy(out, j);
                hit.dE = dE;
                
                auto eDssd = bEnergy(out, j);
                auto ePad = pad.energy(j);
                if(i == 0) ePad = 0;

                hit.deposited = eDssd;
                hit.paddeposited = ePad;

                
                auto BI = bSeg(out, j);
                auto FI = fSeg(out, j);
                hit.fseg = short(FI);
                hit.bseg = short(BI);

                auto position = out.detector().getUniformPixelPosition(FI, BI);
                auto origin = target.getCenter() - TVector3(0, 0, targetThickness / 2.) + TVector3(0, 0, IonRange);
                auto direction = (position - origin).Unit();
                auto angle = direction.Angle(-det.getNormal());

                hit.position = position;
                hit.direction = direction;
                hit.origin = origin;
                hit.angle = angle;
                hit.theta = hit.direction.Theta();
                hit.phi = hit.direction.Phi();


                hit.TF = fTime(out, j);
                hit.TB = bTime(out, j);
                hit.TPad = pad.time(0);

                hit.detectorMul = m;
                hit.numInDet = j;
                hit.index = i;
                hits.emplace_back(move(hit));
            }
        }
    }

    void determineTypes() {
        for(int i = 0; i < hits.size(); i++){
            Hit *hit = &hits[i];
            hit->canBeAlpha = false;
            hit->canBeBeta = false;
            bool thickDetector = (hit->index == 1 || hit->index == 5);      // Thick detectors are the ones with id 1 or 5 (det2 and detD)
            if(hit->paddeposited != 0 || thickDetector){                    // Can be a beta if it comes from a thick detector or a pad  
                hit->canBeBeta = true;
                numberOfBetas++;
                hit->EBeta = hit->deposited + hit->paddeposited;
                if(hit->detectorMul == 1){
                    // If there only one particle has hit, and it is a beta, it cannot be an alpha
                    
                    // hit->lVector = {sqrt(2 * hit->paddeposited * ELECTRON_MASS) * hit->direction, hit->paddeposited + ELECTRON_MASS};
                    // continue;
                }
            }

            hit->canBeAlpha = true;
            numberOfAlphas++;
            auto tF = deadlayerF[hit->index] / abs(cos(hit->angle));
            double E = hit->deposited;
            E += SiCalc->getTotalEnergyCorrection(E, tF);
            for (auto &intersection : target.getIntersections(hit->position, hit->origin)) {
                auto &calc = targetCalcs[intersection.index];
                E += calc->getTotalEnergyCorrection(E, intersection.transversed);
                hit->Ectarget += calc->getTotalEnergyCorrection(E, intersection.transversed);
                hit->tarTrav += intersection.transversed;
            }
            hit->E = E;
            hit->lVector = {sqrt(2 * E * ALPHA_MASS) * hit->direction, E + ALPHA_MASS}; 

        }
    }

    vector<int> determineDoubleAlpha() {
        vector<int> index = {0, 0};
        double angle = 100;
        double momentumSum = std::numeric_limits<double>::infinity();
        for(int i = 0; i < hits.size(); i++){
            for(int j = i +1; j < hits.size(); j++){
                if(!hits[i].canBeAlpha && !hits[j].canBeAlpha) continue;

                Hit *hit1 = &hits[i];
                Hit *hit2 = &hits[j];
                auto TempMomentumSum = (hits[i].lVector.Vect() + hits[j].lVector.Vect()).Mag();
                if (TempMomentumSum < momentumSum) {
                    index[0] = i;
                    index[1] = j;
                    momentumSum = TempMomentumSum;
                }
                
                moment = momentumSum;

                // hits[index[0]].canBeBeta = false;
                // hits[index[1]].canBeBeta = false;

            }
        }
        return index;
        
    }
    
    double angleBetweenTwoHits(Hit *hit1, Hit *hit2){
        double x1 = hit1->direction.X();
        double x2 = hit2->direction.X();

        double y1 = hit1->direction.Y();
        double y2 = hit2->direction.Y();

        double z1 = hit1->direction.Z();
        double z2 = hit2->direction.Z();

        double dotprod = x1*x2 + y1*y2 + z1*z2;
        double len1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
        double len2 = sqrt(pow(x2, 2) + pow(y2, 2) + pow(z2, 2));
        double cosa = dotprod / (len1 * len2);
        
        return cosa;
    }

    vector<double> getXYZ(double theta, double phi){
        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);
        vector<double> v = {x,y,z};
        return v;
    }

    double angleToBeam(Hit *hit){
        auto beam  = TVector3(0, 0, IonRange).Unit();
        double x1 = hit->direction.X();
        double x2 = beam.X();

        double y1 = hit->direction.Y();
        double y2 = beam.Y();

        double z1 = hit->direction.Z();
        double z2 = beam.Z();

        double dotprod = x1*x2 + y1*y2 + z1*z2;
        double len1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
        double len2 = sqrt(pow(x2, 2) + pow(y2, 2) + pow(z2, 2));
        double cosa = dotprod / (len1 * len2);

        return cosa;
    }

    void doAnalysis() { 
        vector<int> alphaIndex = determineDoubleAlpha();
        int i = alphaIndex[0];
        int j = alphaIndex[1];
        Hit *a0, *a1;
        a0 = &hits[i];
        a1 = &hits[j];
        double cosAlphaAngle = angleBetweenTwoHits(a0, a1);
        cosang = cosAlphaAngle;
        cosA->add(cosang);

        double momentumSum = (a0->lVector.Vect() + a1->lVector.Vect()).Mag();
        double energySum = a0->E + a1->E;
        v_ptot->add(momentumSum);
        v_Etot->add(energySum);
        double Ex = energySum - Be8->mass + 2 * ALPHA_MASS;
        v_Ex->add(Ex);



        vector<Hit *> orderedHits{a0, a1};
        v_pos0->add(a0->position);
        v_dir0->add(a0->direction);
        v_p0->add(a0->lVector.Vect());
        v_pos1->add(a1->position);
        v_dir1->add(a1->direction);
        v_p1->add(a1->lVector.Vect());

        for (auto &hit : orderedHits){
            v_E->add(hit->E);
            v_Edep->add(hit->deposited);
            v_F->add(hit->fseg);
            v_B->add(hit->bseg);
            v_i->add(static_cast<short>(hit->index));
            v_dE->add(hit->dE);
            v_ang->add(hit->angle * TMath::RadToDeg());
            v_angToBeam->add(angleToBeam(hit));

        }
        double betaAng;
        for (int in = 0; in < hits.size(); in++) {
            double kk = angleBetweenTwoHits(&hits[in], &hits[in+1]);
            v_cosAngAll->add(kk);
            Hit *hit = &hits[in];
            
            if (!hit->canBeBeta) continue;
            if (hit == a0) continue;
            if (hit == a1) continue;
            v_betaAlpaAngle0->add(angleBetweenTwoHits(hit, a0));
            v_betaAlpaAngle1->add(angleBetweenTwoHits(hit, a1));
            v_ibeta->add((hit->index));
            v_Ebeta->add(hit->EBeta);
            v_bPhi->add(hit->phi);
            v_bTheta->add(hit->theta);
            mulBeta++;
            v_EPbeta->add(hit->paddeposited);
        }
        if(mulBeta > 0){
            
        }
        mul = hits.size();
        mulAlpha = 2;

    }

    void terminate() override {
        AbstractSortedAnalyzer::terminate();
        gDirectory->WriteTObject(tree);
    }

    void clear() {
        mul = 0;
        mulBeta = 0;
        mulAlpha = 0;
        numberOfBetas = 0;
        numberOfAlphas = 0;
        cosang = 0;
        moment = 0;
        hits.clear();
        AUSA::clear(
            *v_E, *v_Edep,
            *v_F, *v_B,
            *v_Ebeta, *v_ibeta, *v_EPbeta, 
            *v_dE, *v_angToBeam, *v_ang,
            *cosA,
            *v_ptot, *v_Etot,
            *v_betaAlpaAngle1, *v_betaAlpaAngle0,
            *v_bPhi, *v_bTheta, *v_i,
            *v_cosAngAll,
            *v_p0, *v_p1, *v_Ex, *v_dir0, *v_dir1, *v_pos0, *v_pos1
        );



    }



    Target &target;
    ParticleType *ALPHA, *Be8;
    TTree *tree;
    UInt_t mul{}, CLOCK{}, mulAlpha{}, mulBeta{};
    int NUM;
    unique_ptr<DynamicBranchVector<TVector3>> v_dir0, v_pos0, v_dir1, v_pos1, v_p0, v_p1;
    unique_ptr<DynamicBranchVector<double>> v_E, v_Edep, v_dE, v_Ebeta, v_EPbeta, cosA, v_angToBeam, v_ang, v_cosAngAll, v_Ex;
    unique_ptr<DynamicBranchVector<short>> v_F, v_B, v_i;
    unique_ptr<DynamicBranchVector<int>> v_ibeta;
    unique_ptr<DynamicBranchVector<double>> v_ptot, v_Etot, v_betaAlpaAngle0, v_betaAlpaAngle1, v_bTheta, v_bPhi;
    unique_ptr<EnergyLossRangeInverter> SiCalc;
    vector<unique_ptr<EnergyLossRangeInverter>> targetCalcs;
    double cosang, moment;
    vector<double> deadlayerF, deadlayerB, deadlayerP;
    double IonRange, targetThickness;
    vector<Hit> hits;
    int numberOfAlphas, numberOfBetas;


};

void print(vector<string> const &input) {
 
}

int main(int argc, char *argv[]) {
    clock_t start = clock();
    vector<string> input;
    int run = 0;

    auto setup = JSON::readSetupFromJSON("/home/anders/i257/setup/setup.json");
    auto target = JSON::readTargetFromJSON("/home/anders/i257/setup/targets/target.json");
    double targetThickness = target.getThickness() *1e6;

    string file = "/home/anders/i257/data/sorted/225_N102m.root";

    SortedReader reader(*setup);
    reader.add(file);
    auto base = stripFileExtension(extractFileName(file));

    TString outfile = ("/home/anders/i257/data/Li8/" + base + "lio.root").c_str();
    TFile output(outfile, "RECREATE");

    auto analysis = make_shared<Li8Analysis>(target, &output);
    reader.attach(analysis);
    reader.run();


    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("\n\x1b[96mTime elapsed: %.2f min\n", elapsed/60);
}
