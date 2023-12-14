#include "AMPCore.h"
#include "hw/HW8.h"
#include "MyFlightPlanner.h"
#include "HelpfulClass.h"

using namespace amp;

int main(int argc, char** argv) {
    amp::RNG::seed(amp::RNG::randiUnbounded());
    UASProblem prob(3,1,30,0.5,4.0);
    

    MyFlightPlanner fPlanner;
    MultiAgentPath2D soln =  fPlanner.plan(prob);
    // LOG("Now in main, len = " << soln.agent_paths[0].waypoints.size());
    Visualizer::makeFigure(prob,soln);
    LOG("Figure made");
    Visualizer::showFigures();
    LOG("Figure shown");
    // return 0;
}