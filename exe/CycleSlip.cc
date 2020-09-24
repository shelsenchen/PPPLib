//
// Created by cc on 9/15/20.
//

#include "CmnFunc.h"
#include "ReadFiles.h"
#include "Solver.h"
INITIALIZE_EASYLOGGINGPP

namespace plt=matplotlibcpp;
using namespace PPPLib;

int main(int argc,char** argv)
{
    tPPPLibConf C;

//    C.fileC.rover="/home/cc/dataset/data_mgex/sgoc3350.19o";
    vector<int> t;
    vector<double> b;

    for(int i=0;i<10;i++){
        t.push_back(i);
        b.push_back(1.0+i*0.5);
    }

    plt::plot(t,b);
    plt::show()
}
