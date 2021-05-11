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
#include <fstream>
#include <string>


using namespace std;
using namespace ROOT::Math;


// double ang(double theta, double phi){
//     vector<double> v1 = getXYZ(theta, phi);
//     vector<double> v2 = getXYZ(theta, phi);
//     double dotprod = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
//     double len1 = sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
//     double len2 = sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
//     double cosa = dotprod / (len1 * len2);
//     return cosa;
// }

// vector<double> getXYZ(double theta, double phi){
//     double x = sin(theta) * cos(phi);
//     double y = sin(theta) * sin(phi);
//     double z = cos(theta);
//     vector<double> v = {x,y,z};
//     return v;
// }

double angFromXYZ(double x1, double y1, double z1, double x2, double y2, double z2){
    double dotprod = x1*x2 + y1*y2 + z1*z2;
    double len1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
    double len2 = sqrt(pow(x2, 2) + pow(y2, 2) + pow(z2, 2));
    double cosa = dotprod / (len1 * len2);
    return cosa;
}




double pixel_x[1536];
double pixel_y[1536];
double pixel_z[1536];
double angs[2359296];


void createDetector1And2(double distance, double x, double y, int start){
    double y_safe = y;
    int k = start;
    for(int i = 0; i < 16; i++){
        for(int j = 0; j < 16; j++){
            pixel_x[k] = x;
            pixel_y[k] = y;
            pixel_z[k] = distance;
            y += 9;
            k++;
        }
        x += 9;
        y = y_safe;
    }
}

void createDetector3And4(double distance, double z, double y, int start){
    double y_safe = y;
    int k = start;
    for(int i = 0; i < 16; i++){
        for(int j = 0; j < 16; j++){
            pixel_x[k] = distance;
            pixel_y[k] = y;
            pixel_z[k] = z;
            k++;
            y += 9;
        }
        z += 9;
        y = y_safe;
    }
}

void createDetector5And6(double distance, double x, double z, int start){
    double z_safe = z;
    int k = start;
    for(int i = 0; i < 16; i++){
        for(int j = 0; j < 16; j++){
            pixel_x[k] = x;
            pixel_y[k] = distance;
            pixel_z[k] = z;
            z += 9;
            k++;
        }
        x += 9;
        z = z_safe;
    }
}

ofstream outputFile;
ofstream fs;
std::string filename = "efficiencyOutput.csv";

int main(int argc, char const *argv[])
{   
    double init = sqrt(pow(28.7, 2) + pow(28.7, 2) + pow(0.4, 2)); // milimeter
    double x = -67.5;
    double y = -67.5;
    createDetector1And2(init, x, y, 0);
    createDetector1And2(-init, x, y, 256);
    createDetector3And4(init, x, y, 512);
    createDetector3And4(-init, x, y, 768);
    createDetector5And6(init, x, y, 1024);
    createDetector5And6(-init, x, y, 1280);
    int k = 0;
    for(int i = 0; i < 1536; i++){
        for(int j = 0; j < 1536; j++){
            angs[k] = angFromXYZ(pixel_x[i], pixel_y[i], pixel_z[i], pixel_x[j], pixel_y[j], pixel_z[j]);
            // cout << angs[k] << endl;
            k++;
        }
    }


    outputFile.open(filename);
    // outputFile << "data" << endl;
    for(int i = 0; i < sizeof(angs)/sizeof(*angs); i++){
        outputFile << angs[i] << endl;
    }
    outputFile.close();


    return 0;
}

