#include <TVector3.h>
#include <TLorentzVector.h>
#include "ParticleType.h"

struct Hit {
    double deposited;
    double paddeposited;

    double E, FE, BE, Edssd, dE, EBeta, fbdiff;

    unsigned int ADCF, ADCB;
    double deltaPeak;
    double TF, TB, TPad, T;
    double targetAngle;
    double tarTrav, detTrav;
    double Ectarget;
    TVector3 direction, position, origin;
    double theta;
    double angle;
    double phi;
    int numInDet, detectorMul;
    short fseg, bseg;
    std::vector<ParticleType> type;
    bool canBeAlpha;
    bool canBeBeta;
    size_t index;

    TLorentzVector lVector, lVectorBeta, lVectorAlpha;
};
