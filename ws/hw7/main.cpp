#include "AMPCore.h"
#include "hw/HW2.h"
#include "hw/HW5.h"
#include "hw/HW7.h"
#include "MyPointMotionPlanner.h"
#include "HelpfulClass.h"

using namespace amp;

int main(int argc, char** argv) {
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // Visualizer::showFigures();
    // amp::HW6::grade(wf1g, wf2g, Astar, "collin.hudson@colorado.edu", argc, argv);
    // HW7::hint();
    amp::Problem2D w1 = HW5::getWorkspace1();
    w1.y_min = -3;
    w1.y_max = 3;
    w1.x_min = -1;
    w1.x_max = 11;
    amp::Problem2D w2 = HW2::getWorkspace1();
    amp::Problem2D w3 = HW2::getWorkspace2();
    amp::Path2D path;
    std::vector<Eigen::Vector2d> collision_points;
    {LOG("PRM");

        // std::list<std::vector<double>> pathTimes;
        // std::list<std::vector<double>> pathLens;
        // std::vector<double> tempT;
        // std::vector<double> tempL;
        // // std::vector<std::pair<int,double>> NR = {{200,0.5},{200,1},{200,1.5},{200,2},{500,0.5},{500,1},{500,1.5},{500,2}};
        // // std::vector<std::string> NRLabels = {"{200,0.5}","{200,1}","{200,1.5}","{200,2}","{500,0.5}","{500,1}","{500,1.5}","{500,2}"};
        // std::vector<std::pair<int,double>> NR = {{200,1},{200,2},{500,1},{500,2},{1000,1},{1000,2}};
        // std::vector<std::string> NRLabels = {"{200,1}","{200,2}","{500,1}","{500,2}","{1000,1}","{1000,2}"};
        // std::vector<double> solns(NR.size(),0);
        // for(int n = 0; n < NR.size(); n++){
        //     LOG("n: " << NR[n].first << " r: " << NR[n].second);
        //     MyPRM prm;
        //     tempT.clear();
        //     tempL.clear();
        //     prm.getS() = true;
        //     prm.getN() = NR[n].first;
        //     prm.getR() = NR[n].second;
        //     for(int i = 0; i < 100; i++){
        //         path = prm.plan(w3);
        //         if(HW7::check(path,w3,false)){
        //             tempT.push_back(prm.getT());
        //             tempL.push_back(path.length());
        //             solns[n] ++;
        //         }
                
        //     }
        //     pathTimes.push_back(tempT);
        //     pathLens.push_back(tempL);
        // }
        // Visualizer::makeBoxPlot(pathTimes, NRLabels, std::string("W3 runtimes"),std::string("{n,r}"),std::string("Runtime (ms)"));
        // Visualizer::makeBoxPlot(pathLens, NRLabels,std::string("W3 path lengths"), std::string("{n,r}"),std::string("Path Length"));
        // Visualizer::makeBarGraph(solns, NRLabels,std::string("W3 Num Solutions"), std::string("{n,r}"),std::string("#solutions"));

        // MyPRM prm;
        // prm.getN() = 200;
        // prm.getR() = 0.5;
        // prm.getW() = true;
        // path = prm.plan(w1);
        // HW7::check(path,w1,collision_points);
        // Visualizer::makeFigure(w1, path, collision_points);
        // prm.pathSmoother(w1,path);
        // Visualizer::makeFigure(w1, path, collision_points);
        // // do{
        //     collision_points.clear();
        //     path = prm.plan(w2);
        // // }while(!HW7::check(path,w2,collision_points,false));
        // // HW7::check(path,w2,collision_points);
        // // Visualizer::makeFigure(w2, path, collision_points);
        // // prm.pathSmoother(w2,path);
        // // Visualizer::makeFigure(w2, path, collision_points);
        
        // // do{
        //     collision_points.clear();
        //     path = prm.plan(w3);
        // // }while(!HW7::check(path,w3,collision_points,false));
        // // HW7::check(path,w3,collision_points);
        // // Visualizer::makeFigure(w3, path, collision_points);
        // // prm.pathSmoother(w3,path);
        // // Visualizer::makeFigure(w3, path, collision_points);
    }
    {LOG("RRT");
        MyGoalBiasRRT RRT;
        RRT.getN() = 5000;
        RRT.getG() = 0.05;
        RRT.getS() = 0.5;
        RRT.getE() = 0.25;
        RRT.getW() = true;
        collision_points.clear();
        path = RRT.plan(w1);
        HW7::check(path,w1,collision_points);
        // Visualizer::makeFigure(w1, path, collision_points);
        collision_points.clear();
        path = RRT.plan(w2);
        HW7::check(path,w2,collision_points);
        // Visualizer::makeFigure(w2, path, collision_points);
        collision_points.clear();
        path = RRT.plan(w3);
        HW7::check(path,w3,collision_points);
        // Visualizer::makeFigure(w3, path, collision_points);
        
        // std::list<std::vector<double>> pathTimes;
        // std::list<std::vector<double>> pathLens;
        // std::vector<double> tempT;
        // std::vector<double> tempL;
        // std::vector<amp::Problem2D> NR = {w1,w2,w3};
        // std::vector<std::string> NRLabels = {"W1","W2","W3"};
        // std::vector<double> solns(NR.size(),0);
        // for(int n = 0; n < NR.size(); n++){
        //     MyGoalBiasRRT RRT;
        //     tempT.clear();
        //     tempL.clear();
        //     for(int i = 0; i < 100; i++){
        //         path = RRT.plan(NR[n]);
        //         if(HW7::check(path,NR[n],false)){
        //             tempT.push_back(RRT.getT());
        //             tempL.push_back(path.length());
        //             solns[n] ++;
        //         }
                
        //     }
        //     pathTimes.push_back(tempT);
        //     pathLens.push_back(tempL);
        // }
        // Visualizer::makeBoxPlot(pathTimes, NRLabels, std::string("RRT Runtimes"),std::string("Workspace"),std::string("Runtime (ms)"));
        // Visualizer::makeBoxPlot(pathLens, NRLabels,std::string("RRT Path lengths"), std::string("Workspace"),std::string("Path Length"));
        // Visualizer::makeBarGraph(solns, NRLabels,std::string("RRT Num Solutions"), std::string("Workspace"),std::string("#solutions"));
    }
    // Visualizer::showFigures();
    // MyPRM prm;
    // MyGoalBiasRRT RRT;
    // amp::HW7::grade(prm, RRT, "collin.hudson@colorado.edu", argc, argv);

    return 0;
}