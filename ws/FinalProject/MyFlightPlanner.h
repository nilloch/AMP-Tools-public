#pragma once

#include "AMPCore.h"
#include "hw/HW7.h"
#include "hw/HW8.h"
#include <Eigen/LU>
#include "HelpfulClass.h"

using Node = uint32_t;

struct UASProblem : public amp::MultiAgentProblem2D {

    amp::MultiAgentPath2D GApaths; //xy position of ground agents
    std::vector<int> endGAt; //vector of times (number of steps) for each GA to reach their final position
    int maxTime = 1; //Largest element in endGAt
    int numGA = 1;
    int numUAV = 1;
    UASProblem(uint32_t n_GA = 3, uint32_t n_UAV = 2, uint32_t n_Obs = 10, double min_Obs = 1.0, double max_Obs = 2.0, double size_UAV = 0.2);
};

class MyFlightPlanner : public MyGoalBiasRRTND{
    public:
        /// @brief Solve a motion planning problem. Derive class and override this method
        /// @param problem Multi-agent motion planning problem
        /// @return Array of paths that are ordered corresponding to the `agent_properties` field in `problem`.
        amp::MultiAgentPath2D plan(const UASProblem& problem);
};

class FlightChecker : public checkPath{
    public:
        bool inLOS(const Eigen::Vector2d state0, const Eigen::Vector2d state1, const amp::Environment2D& obs){
            return !lineCollision2D(state0, state1, obs);
        };

        void makeLOS(const UASProblem& problem);

        void updateLOS(Eigen::VectorXd state, const UASProblem& problem, int time);

        bool checkLOS(int numGA);

    private:
        std::vector<std::set<int>> losGraph;

};

