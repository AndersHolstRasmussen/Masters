// For handling command line arguments (getopt)
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>

// Parses setup.json file
#include "ausa/json/IO.h"
#include "ausa/util/FileUtil.h"

// Reads Sorted output
#include "ausa/sort/SortedReader.h"

// Library analyzer
#include "ausa/sort/analyzer/SummedSpectrumPlotter.h"
#include "PadGatedSummedTelescope.h"
#include "TreeMaker.h"

// ROOT classes
#include <TH2.h>
#include <TMath.h>

using namespace std;
using namespace AUSA;
using namespace AUSA::Sort;

const int len1 = 256;
const int len2 = 512;
const int len6 = 1536;
const int len6squared = 2359296;
const int len2times6 = len2 * len6;

double xs[len6];
double ys[len6];
double zs[len6];
TVector3 cords[len6];
TVector3 norms[len6];
double angs[len6squared];
double weight[len6squared];
double alphaAngs[len6];
double singleAng[256], singleW[256];
double det1w[256];
ofstream outputFile;
ofstream outputFile2;
ofstream outputFile3;
std::string filename = "efficiencyOutputAllDet.csv";
std::string filename2 = "effDet2.csv";
std::string filename3 = "angEffDet2.csv";


double angFromXYZ(double x1, double y1, double z1, double x2, double y2, double z2){
    double dotprod = x1*x2 + y1*y2 + z1*z2;
    double len1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
    double len2 = sqrt(pow(x2, 2) + pow(y2, 2) + pow(z2, 2));
    double cosa = dotprod / (len1 * len2);
    return cosa;
}

double ang(TVector3 v1, TVector3 v2){
    // auto dot = v1.Dot(v2);
    v1.Angle(v2);
    return v1.Dot(v2) / (v1.Mag() * v2.Mag());
}

double efficiency(double x, double y, double z, TVector3 norm){
    // TVector3 direction = TVector3(x, y, z);
    // auto mag = direction.Mag();
    // auto magnorm = norm.Mag();
    // double theta = direction.Dot(norm) / mag*magnorm;
    // return cos(theta) / mag*mag;

    double len = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double xn = -norm.X();
    double yn = -norm.Y();
    double zn = -norm.Z();
    double costheta = angFromXYZ(x, y, z, xn, yn, zn);
    return (costheta) / (pow(len, 2));
}


int main(int argc, char **argv) {
    // Read in setup configuration
    auto setup = JSON::readSetupFromJSON("/home/anders/i257/setup/setup.json");
    auto target = JSON::readTargetFromJSON("/home/anders/i257/setup/targets/target.json");
    double IonRange = 166 * 1e-6;
    double targetThickness = 226 * 1e-6;
    auto beam  = TVector3(0, 0, IonRange).Unit();
    auto tarPos = target.getCenter() - TVector3(0, 0, targetThickness / 2.) + TVector3(0, 0, IonRange);

 
    int k = 0;
    int dets[6] = {1, 5, 0, 2, 3, 4};
    for(int d : dets){
        for(int i = 1; i < 17; i++){
            for(int j = 1; j < 17; j++){
                auto pixPos = setup->getDSSD(d)->getPixelPosition(i, j);
                auto direction = (pixPos - tarPos);
                auto norm = setup->getDSSD(d)->getNormal();
                auto theta = direction.Theta();
                auto phi = direction.Phi(); 
                norms[k] = norm;          
                xs[k] = direction.x();
                ys[k] = direction.y();
                zs[k] = direction.z();
                k++;
            }
        }
    }

    for(int i = 0; i < 256; i++){
        det1w[i] = efficiency(xs[i], ys[i], zs[i], norms[i]);
    }

    k = 0;
    // for (int i = 0; i < 256; i++){
    //     for (int j = 0; j < 256; j++){
    //         singleAng[k] = angFromXYZ(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j]);
    //         singleW[k] = efficiency(xs[i], ys[i], zs[i], norms[i]) * efficiency(xs[j], ys[j], zs[j], norms[j]);
    //         k++;
    //     }
    // }
    

    k = 0;
    for(int i = 0; i < len6; i++){
        for(int j = 0; j < len6; j++){
            angs[k] = angFromXYZ(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j]);
            weight[k] = efficiency(xs[i], ys[i], zs[i], norms[i]) * efficiency(xs[j], ys[j], zs[j], norms[j]);
            k++;
        }
    }

    // k = 0;
    // for (int i = 0; i < len6; i++){
    //     alphaAngs[k] = angFromXYZ(xs[i], ys[i], zs[i], beam.X(), beam.Y(), beam.Z());
    //     k++;
    // }
    
    // outputFile3.open(filename3);
    // for(int i = 0; i < 256*256; i++){
    //     outputFile3 << singleAng[i] << "\t" << singleW[i] << endl;
    // }
    // outputFile3.close();

    outputFile2.open(filename2);
    for(int i = 0; i < 256; i++){
        outputFile2 << det1w[i] << endl;
    }
    outputFile2.close();

    outputFile.open(filename);
    for(int i = 0; i < len6squared; i++){
        outputFile << angs[i] << "\t" << weight[i] << endl; // 
    }
    outputFile.close();

    return 0;
}


