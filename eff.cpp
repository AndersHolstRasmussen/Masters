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
double angs[len6squared];

ofstream outputFile;
std::string filename = "efficiencyOutput.csv";


double angFromXYZ(double x1, double y1, double z1, double x2, double y2, double z2){
    double dotprod = x1*x2 + y1*y2 + z1*z2;
    double len1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
    double len2 = sqrt(pow(x2, 2) + pow(y2, 2) + pow(z2, 2));
    double cosa = dotprod / (len1 * len2);
    return cosa;
}


int main(int argc, char **argv) {
    // Read in setup configuration
    auto setup = JSON::readSetupFromJSON("/home/anders/i257/setup/setup.json");
    auto target = JSON::readTargetFromJSON("/home/anders/i257/setup/targets/target.json");
    auto tarPos = target.getCenter();

 
    int k = 0;
    int dets[6] = {1, 5, 0, 2, 3, 4};
    for(int d : dets){
        for(int i = 1; i <= 16; i++){
            for(int j = 1; j <= 16; j++){
                auto pixPos = setup->getDSSD(d)->getPixelPosition(i, j);
                auto direction = (pixPos - tarPos).Unit();
                xs[k] = direction.X();
                ys[k] = direction.Y();
                zs[k] = direction.Z();
                k++;
            }
        }
    }

    k = 0;
    for(int i = 0; i < len2; i++){
        for(int j = 0; j < len6; j++){
            angs[k] = angFromXYZ(xs[i], ys[i], zs[i], xs[j], ys[j], zs[j]);
            k++;
        }
    }

    outputFile.open(filename);
    for(int i = 0; i < len2times6; i++){
        outputFile << angs[i] << endl;
    }
    outputFile.close();

    return 0;
}


