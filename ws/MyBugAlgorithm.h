#pragma once

#include "AMPCore.h"
#include "hw/HW2.h"

/// @brief Declare your bug algorithm class here. Note this class derives the bug algorithm class declared in HW2.h
class MyBugAlgorithm : public amp::BugAlgorithm {
    public:
        Eigen::Vector2d bugXY;

        Eigen::Vector2d bugNext;

        double bugStep = 0.05;

        amp::Polygon currentOb;

        int vertIdx;

        Eigen::Vector2d wallPoint;

        Eigen::Vector2d qLeave;
        
        double qLScore;
        
        bool hit = false;

        bool followMode;

        
        // Override and implement the bug algorithm in the plan method. The methods are declared here in the `.h` file
        virtual amp::Path2D plan(const amp::Problem2D& problem) override;

        // Add any other methods here...
        
        Eigen::Vector2d getNext(Eigen::Vector2d bugXY, Eigen::Vector2d goal, double stepSize){
            Eigen::Vector2d nextXY;
            double ySign = (1 - 2*(goal(1) - bugXY(1) < 0));
            double xSign = (1 - 2*(goal(0) - bugXY(0) < 0));
            //std::cout << "ySign: " << ySign << "xSign: " << xSign << std::endl;
            //If line to goal is vertical, move vertically by step size
            if((goal(0) - bugXY(0)) == 0){
                nextXY << bugXY(0), bugXY(1) + stepSize*ySign;
            }
            //If line to goal is horizontal, move horizontally by step size
            else if((goal(1) - bugXY(1)) == 0){
                nextXY << bugXY(0) + stepSize*xSign, bugXY(1);
            }
            else{
                //double sign = (1 - 2*(goal(0) - bugXY(0) < 0));
                //nextXY << bugXY(0) + stepSize*(xSign), (goal(1) - bugXY(1))/(goal(0) - bugXY(0))*(bugXY(0) + stepSize*(xSign));
                double theta = atan2((goal(1) - bugXY(1)),(goal(0) - bugXY(0)));
                //std::cout << "theta: " << theta << " cos(theta) = " << cos(theta) << " sin(theta) = "<< sin(theta) << std::endl;
                nextXY << bugXY(0) + cos(theta)*stepSize, bugXY(1) + sin(theta)*stepSize;
            }
            return nextXY;
        }
        Eigen::Vector2d getNextFollow(Eigen::Vector2d bugXY, Eigen::Vector2d obV1, Eigen::Vector2d obV2, double stepSize){
            // obV2 must be clockwise next from obV1!
            Eigen::Vector2d nextXY;
            double ySign = (1 - 2*(obV2(1) - obV1(1) < 0));
            double xSign = (1 - 2*(obV2(0) - obV1(0) < 0));
            //If wall to follow is vertical, move vertically by step size
            if((obV2(0) - obV1(0)) == 0){
                nextXY << bugXY(0), bugXY(1) + stepSize*ySign;
            }
            //If wall to follow is horizontal, move horizontally by step size
            else if((obV2(1) - obV1(1)) == 0){
                nextXY << bugXY(0) + stepSize*xSign, bugXY(1);
            }
            else{
                //nextXY << bugXY(0) + stepSize*(xSign), (obV2(1) - obV1(1))/(obV2(0) - obV1(0))*(bugXY(0) + stepSize*(xSign));
                double theta = atan2((obV2(1) - obV1(1)),(obV2(0) - obV1(0)));
                //std::cout << "theta: " << theta << std::endl;
                nextXY << bugXY(0) + cos(theta)*stepSize, bugXY(1) + sin(theta)*stepSize;
            }
            return nextXY;
        }
        Eigen::Vector2d getNextCorner(Eigen::Vector2d bugXY, Eigen::Vector2d corner, double theta){
            // nextPoint must be clockwise next from corner!
            Eigen::Vector2d nextXY;
            //nextXY << corner(0) + bugStep*sin(theta), corner(1) + bugStep*cos(theta);
            nextXY << corner(0) + bugStep*cos(theta), corner(1) + bugStep*sin(theta);
            return nextXY;
        }
        double getT(Eigen::Vector2d bugXY, Eigen::Vector2d bugNext, Eigen::Vector2d v1, Eigen::Vector2d v2){
            Eigen::Matrix2d Tmat1;
            Eigen::Matrix2d Tmat2;
            Tmat1.col(0) = bugXY-v1;
            Tmat1.col(1) = v1-v2;
            Tmat2.col(0) = bugXY-bugNext;
            Tmat2.col(1) = v1-v2;
            return Tmat1.determinant()/Tmat2.determinant();
        }
        double getU(Eigen::Vector2d bugXY, Eigen::Vector2d bugNext, Eigen::Vector2d v1, Eigen::Vector2d v2){
            Eigen::Matrix2d Tmat1;
            Eigen::Matrix2d Tmat2;
            Tmat1.col(0) = bugXY-v1;
            Tmat1.col(1) = bugXY-bugNext;
            Tmat2.col(0) = bugXY-bugNext;
            Tmat2.col(1) = v1-v2;
            return Tmat1.determinant()/Tmat2.determinant();
        }
        bool evalTU( double t, double u){
            return t >= 0 && t <= 1 && u >= 0 && u <= 1;
        }
        Eigen::Vector2d getHit(Eigen::Vector2d bugXY, Eigen::Vector2d bugNext, double t){
            Eigen::Vector2d interP;
            interP << bugXY(0) + t*(bugNext(0) - bugXY(0)), bugXY(1) + t*(bugNext(1) - bugXY(1));
            return interP;
        }
        void stepCheck(const amp::Problem2D& problem){
            hit = false;
            for(auto Ob : problem.obstacles){
                if(!hit){
                    for(int j = 0; j < Ob.verticesCCW().size(); j++){
                        if(!hit){
                            hit = checkHit(Ob,j); //hitPoint set in checkHit function
                        }
                    }
                }
            }
        }
        bool checkHit(amp::Polygon Ob, int startVert){
            int endVert;
            std::vector<Eigen::Vector2d> vertCW = Ob.verticesCW();
            int numV = vertCW.size() - 1;
            (startVert == 0) ? (endVert = numV) : endVert = startVert - 1;
            double t = getT(bugXY,bugNext,vertCW[startVert],vertCW[endVert]);
            double u = getU(bugXY,bugNext,vertCW[startVert],vertCW[endVert]);
            if((bugNext - vertCW[startVert]).norm() < bugStep/2){
                //std::cout << "hitpoint: " << getHit(bugXY,bugNext,t) << std::endl;
                wallPoint = vertCW[startVert];
                /*if(!followMode){
                    //qHit = vertCW[startVert];
                }*/
                currentOb = Ob;
                vertIdx = abs(startVert - numV);
                /*
                std::cout << "Hit corner, current obstacle 1st ccw vertex: " << currentOb.verticesCCW()[0] << std::endl;
                std::cout << "vertIdx: "<< vertIdx << std::endl;
                std::cout << "current obstacle vertIdx vertex: " << currentOb.verticesCCW()[vertIdx] << std::endl;
                std::cout << "endVert: "<< endVert <<std::endl;*/
                return true;
            }
            else if(evalTU(t,u)){
                //std::cout << "t: " << t << " u: " << u << " startVert: " << startVert<< " endVert: " << endVert << " hitpoint: " << getHit(bugXY,bugNext,t) << std::endl;
                //std::cout << " startVert Vertex: " << vertCW[startVert] << " endVert Vertex: " << vertCW[endVert] << std::endl;
                wallPoint = getHit(bugXY,bugNext,t);
                /*if(!followMode){
                    //qHit = getHit(bugXY,bugNext,t);
                    qHit = getNext(bugXY,problem.q_goal,getAlignDist(false));
                }*/
                currentOb = Ob;
                //std::cout << "current obstacle 1st ccw vertex: " << currentOb.verticesCCW()[0] << std::endl;
                if((startVert == numV && endVert == 0) || (startVert == 0 && endVert == numV)){
                    vertIdx = 0;
                }
                else{
                    startVert = abs(startVert - numV);
                    endVert = abs(endVert - numV);
                    (startVert > endVert) ? (vertIdx = startVert) : (vertIdx = endVert); 
                }
                //std::cout << "vertIdx: "<< vertIdx << std::endl;
                //std::cout << "current obstacle vertIdx vertex: " << currentOb.verticesCCW()[vertIdx] << std::endl;
                //std::cout << "endVert: "<< endVert <<std::endl;
                return true;
            }
            return false;
        }
        bool checkMLine(Eigen::Vector2d init, Eigen::Vector2d goal){
            // Checks for intersection of next step and a closer M-line point
            double t = getT(bugXY,bugNext,init,goal);
            double u = getU(bugXY,bugNext,init,goal);
            if(evalTU(t,u)){
                Eigen::Vector2d testLeave = getHit(bugXY,bugNext,t);
                if(checkLeavePoint(testLeave,goal)){
                    // Set next step to be on M-line and to stop following the boundary
                    bugNext = testLeave;
                    followMode = false;
                    hit = true;
                    return true;
                }

            }
            return false;
        }
        double getQLScore(Eigen::Vector2d bugXY, Eigen::Vector2d goal){
            Eigen::Vector2d diff =  (goal - bugXY);
            return diff.norm();
        }
        bool checkLeavePoint(Eigen::Vector2d position, Eigen::Vector2d goal){
            // returns true if leavePoint changed
            double testScore = getQLScore(position,goal);
            if(testScore <= qLScore || qLScore == 0 ){
                qLeave = position;
                qLScore = testScore;
                return true;
            }
            else{
                return false;
            }
        }
        double getAlignDist(bool alignToCorner){
            if(alignToCorner){
                //assumes bug is bugStep away from the wall and in wall follow mode
                int endVert;
                (vertIdx == 0) ? (endVert = currentOb.verticesCCW().size() - 1) : endVert = vertIdx - 1;
                Eigen::Vector2d corner = currentOb.verticesCCW()[endVert];
                double norm = ((bugXY - corner).norm());
                double theta = acos(bugStep/norm);
                return norm*(sin(theta));                
            }
            else{
                Eigen::Vector2d pointD = currentOb.verticesCCW()[vertIdx];
                double theta1 = atan2((wallPoint(1) - pointD(1)),(wallPoint(0) - pointD(0)));
                double theta5 = atan2((wallPoint(1) - bugXY(1)),(wallPoint(0) - bugXY(0)));
                double theta2;
                //std::cout << "WALLPOINT: " << wallPoint << " pointD: " << pointD << " bugXY: " << bugXY << std::endl;
                if((wallPoint(1) - pointD(1)) == 0 && (wallPoint(0) - pointD(0)) == 0){
                    //moving towards corner
                    int endVert;
                    (vertIdx == 0) ? (endVert = currentOb.verticesCCW().size() - 1) : endVert = vertIdx - 1;
                    pointD = currentOb.verticesCCW()[endVert];
                    theta1 = atan2((wallPoint(1) - pointD(1)),(wallPoint(0) - pointD(0)));
                    theta2 = theta1 - theta5;
                }
                if((wallPoint(1) - pointD(1)) == 0){
                    //wall parallel to x-axis
                    theta2 = theta5;
                }
                else if((wallPoint(0) - pointD(0)) == 0){
                    //wall parallel to y-axis
                    theta2 = atan2((wallPoint(0) - bugXY(0)),(wallPoint(1) - bugXY(1)));
                }
                else{
                    theta2 = theta1 - theta5;
                }
                //std::cout << "theta1: " << theta1 << " theta5: " << theta5 << " theta2: " << theta2 << " bugStep/sin(theta2): " << bugStep/sin(theta2) << " norm: " << (wallPoint - bugXY).norm() << " output: " << (wallPoint - bugXY).norm() - bugStep/sin(theta2) << std::endl;
                return (wallPoint - bugXY).norm() - abs(bugStep/sin(theta2));
            }
            
        }
        amp::Path2D followBug1(const amp::Problem2D& problem, amp::Path2D path){
            int startVert;
            int endVert;
            int stepsAway = 0;
            //move once more toward obstacle edge if able
            bugNext = getNext(bugXY,problem.q_goal,getAlignDist(false));
            stepCheck(problem);//hit set to false at start of stepCheck
            if(!hit){
                path.waypoints.push_back(bugNext);
                bugXY = bugNext;
                hit = true;
            }
            // Set distance from bug to goal at its probable max (initial hit point)
            checkLeavePoint(bugXY, problem.q_goal);
            //Switch to boundary following mode
            followMode = true;
            // set exit point to initial hit point
            Eigen::Vector2d exitPoint = bugXY;
            //std::cout << "exitPoint: " << exitPoint << std::endl;
            do{
                    startVert = vertIdx;
                    (startVert == 0) ? (endVert = currentOb.verticesCCW().size() - 1) : endVert = startVert - 1;
                    //std::cout << " Headed to vertex " << endVert << " With start vert: "<< startVert << std::endl;
                    //step toward nearest clockwise vertex on obstacle
                    bugNext = getNextFollow(bugXY,currentOb.verticesCCW()[startVert],currentOb.verticesCCW()[endVert],2*bugStep);
                    stepCheck(problem); //hit set to false at start of stepCheck
                    if(!hit)
                    {
                        bugNext = getNextFollow(bugXY,currentOb.verticesCCW()[startVert],currentOb.verticesCCW()[endVert],bugStep);
                        path.waypoints.push_back(bugNext);
                        bugXY = bugNext;
                        checkLeavePoint(bugXY, problem.q_goal);
                        stepsAway++;
                        //std::cout << "Norm to end: " << (bugXY - currentOb.verticesCCW()[endVert]).norm() << std::endl;

                        //Check for corner (end of wall line segment)
                        if((bugXY - currentOb.verticesCCW()[endVert]).norm() <= 1.25*bugStep){
                            int nextVert;
                            (endVert == 0) ? (nextVert = currentOb.verticesCCW().size() - 1) : nextVert = endVert - 1;
                            double phi = atan2(currentOb.verticesCCW()[nextVert](1) - currentOb.verticesCCW()[endVert](1), 
                            currentOb.verticesCCW()[nextVert](0) - currentOb.verticesCCW()[endVert](0)) - 
                            atan2(currentOb.verticesCCW()[vertIdx](1) - currentOb.verticesCCW()[endVert](1),
                            currentOb.verticesCCW()[vertIdx](0) - currentOb.verticesCCW()[endVert](0));
                            double theta0 = atan2(bugXY(1) - currentOb.verticesCCW()[endVert](1), bugXY(0) - currentOb.verticesCCW()[endVert](0)) - M_PI/2;
                            double thetaEnd = theta0 + M_PI/2;
                            //Check that taking the corner is safe
                            bugNext = getNextFollow(bugXY,currentOb.verticesCCW()[startVert],currentOb.verticesCCW()[endVert], 2*bugStep);
                            stepCheck(problem); //hit set to false at start of stepCheck
                            int endVertTemp;
                            (endVert == 0) ? (endVertTemp = currentOb.verticesCCW().size() - 1) : endVertTemp = endVert - 1;
                            bugNext = getNextFollow(bugNext,currentOb.verticesCCW()[endVert],currentOb.verticesCCW()[endVertTemp], 2*bugStep);
                            stepCheck(problem);
                            //Take smooth corner
                            if(!hit){
                                for(double  theta = thetaEnd; theta >= theta0; theta = theta - abs(M_PI - abs(phi))/15){
                                    if(!hit){
                                        bugNext = getNextCorner(bugXY, currentOb.verticesCCW()[endVert], theta);
                                        stepCheck(problem);
                                        path.waypoints.push_back(bugNext);
                                        bugXY = bugNext;
                                        checkLeavePoint(bugXY, problem.q_goal);
                                        if(!followMode && (bugXY - exitPoint).norm() < 1.05*bugStep && stepsAway > 2){
                                            //std::cout << "hit leave point on corner, bugXY: "<< bugXY << std::endl;
                                            hit = true; //Indicates we've hit the leave point
                                        }
                                    }
                                    else{
                                        // hit wall while cornering
                                        if(followMode){
                                            std::cout << "Error: hit wall while cornering!" << std::endl;
                                            path.waypoints.push_back(problem.q_goal);
                                            return path;
                                        }
                                    }
                                }
                                // update 'wall starts from vertex' to wall segment end vertex
                                checkLeavePoint(bugXY, problem.q_goal);
                                vertIdx = endVert;
                            }
                        }
                        // Failsafe
                        /*
                        if(abs(bugXY(1)) > 30 || abs(bugXY(0)) > 30 || stepsAway > 5000){
                            std::cout << "Bug has left the working area!" << " bugxy: " << bugXY << std::endl;
                            path.waypoints.push_back(problem.q_goal);
                            return path;
                        }*/
                    }
                    else{
                        //move once more toward obstacle wall if able
                        bugNext = getNext(bugXY,wallPoint,getAlignDist(false));
                        stepCheck(problem);//hit set to false at start of stepCheck
                        if(!hit){
                            path.waypoints.push_back(bugNext);
                            bugXY = bugNext;
                            hit = true;
                        }
                    }
                    // Check if bug has hit original hit point
                    if(followMode && (bugXY - exitPoint).norm() < 1.05*bugStep && stepsAway > 2){
                        // Update point of loop exit to leave point
                        std::cout << "Now retracing to leave point because a) " << (bugXY - exitPoint).norm() << " and b) " << stepsAway << std::endl;
                        std::cout << "Leave point score:  " << qLScore << " and leavePoint: "<< qLeave << std::endl;
                        exitPoint = qLeave;
                        followMode = false;
                    }
                    else if(!followMode && (bugXY - problem.q_goal).norm() <= 1.05*qLScore && stepsAway > 2){
                        //std::cout << "DEBUG: Hit leave point! qLScore: " << qLScore << " norm: " << (bugXY - problem.q_goal).norm() << " bugXY: " << bugXY << " qLeave: " << qLeave << std::endl;
                    }
                }while((bugXY - exitPoint).norm() >= 1.05*bugStep || stepsAway < 2);
                std::cout << "Now leaving object because a) " << (bugXY - exitPoint).norm() << " and b) " << stepsAway << std::endl;
                std::cout << "Leave point score:  " << qLScore << " and leavePoint: "<< qLeave << std::endl;
                // reset leave point score
                qLScore = 0;
                return path;
        }
        amp::Path2D followBug2(const amp::Problem2D& problem, amp::Path2D path){
            int startVert;
            int endVert;
            int stepsAway = 0;
            //move once more toward obstacle edge if able
            bugNext = getNext(bugXY,problem.q_goal,getAlignDist(false));
            stepCheck(problem);//hit set to false at start of stepCheck
            if(!hit){
                path.waypoints.push_back(bugNext);
                bugXY = bugNext;
                hit = true;
            }
            // Set distance from bug to goal at its probable max (initial hit point)
            checkLeavePoint(bugXY, problem.q_goal);
            //Switch to boundary following mode
            followMode = true;
            // set exit point to initial hit point
            Eigen::Vector2d exitPoint = bugXY;
            do{
                    startVert = vertIdx;
                    (startVert == 0) ? (endVert = currentOb.verticesCCW().size() - 1) : endVert = startVert - 1;
                    //std::cout << " Headed to vertex " << endVert << " With start vert: "<< startVert << std::endl;
                    //step toward nearest clockwise vertex on obstacle
                    bugNext = getNextFollow(bugXY,currentOb.verticesCCW()[startVert],currentOb.verticesCCW()[endVert],2*bugStep);
                    stepCheck(problem); //hit set to false at start of stepCheck
                    if(!hit)
                    {
                        bugNext = getNextFollow(bugXY,currentOb.verticesCCW()[startVert],currentOb.verticesCCW()[endVert],bugStep);
                        // Check that bug is sufficiently far away from initial hit point before checking for m-line intersections
                        if(followMode && (bugXY - exitPoint).norm() > 1.05*bugStep){
                            checkMLine(problem.q_init, problem.q_goal);
                        }
                        path.waypoints.push_back(bugNext);
                        bugXY = bugNext;
                        stepsAway++;

                        //Check for corner (end of wall line segment)
                        if((bugXY - currentOb.verticesCCW()[endVert]).norm() <= 1.25*bugStep){
                            int nextVert;
                            (endVert == 0) ? (nextVert = currentOb.verticesCCW().size() - 1) : nextVert = endVert - 1;
                            double phi = atan2(currentOb.verticesCCW()[nextVert](1) - currentOb.verticesCCW()[endVert](1), 
                            currentOb.verticesCCW()[nextVert](0) - currentOb.verticesCCW()[endVert](0)) - 
                            atan2(currentOb.verticesCCW()[vertIdx](1) - currentOb.verticesCCW()[endVert](1),
                            currentOb.verticesCCW()[vertIdx](0) - currentOb.verticesCCW()[endVert](0));
                            double theta0 = atan2(bugXY(1) - currentOb.verticesCCW()[endVert](1), bugXY(0) - currentOb.verticesCCW()[endVert](0)) - M_PI/2;
                            double thetaEnd = theta0 + M_PI/2;
                            //Check that taking the corner is safe
                            bugNext = getNextFollow(bugXY,currentOb.verticesCCW()[startVert],currentOb.verticesCCW()[endVert], 2*bugStep);
                            stepCheck(problem); //hit set to false at start of stepCheck
                            int endVertTemp;
                            (endVert == 0) ? (endVertTemp = currentOb.verticesCCW().size() - 1) : endVertTemp = endVert - 1;
                            bugNext = getNextFollow(bugNext,currentOb.verticesCCW()[endVert],currentOb.verticesCCW()[endVertTemp], 2*bugStep);
                            stepCheck(problem);
                            //Take smooth corner
                            if(!hit){
                                for(double  theta = thetaEnd; theta >= theta0; theta = theta - abs(M_PI - abs(phi))/15){
                                    if(!hit){
                                        bugNext = getNextCorner(bugXY, currentOb.verticesCCW()[endVert], theta);
                                        stepCheck(problem);
                                        checkMLine(problem.q_init, problem.q_goal);
                                        path.waypoints.push_back(bugNext);
                                        bugXY = bugNext;
                                    }
                                    else{
                                        // hit wall while cornering
                                        if(followMode){
                                            std::cout << "Error: hit wall while cornering!" << std::endl;
                                            path.waypoints.push_back(problem.q_goal);
                                            return path;
                                        }
                                    }
                                }
                                // update 'wall starts from vertex' to wall segment end vertex
                                vertIdx = endVert;
                            }
                        }
                        // Failsafe
                        /*
                        if(abs(bugXY(1)) > 30 || abs(bugXY(0)) > 30 || stepsAway > 10000){
                            std::cout << "Bug has left the working area!" << " bugxy: " << bugXY << std::endl;
                            path.waypoints.push_back(problem.q_goal);
                            return path;
                        }*/
                    }
                    else{
                        //move once more toward obstacle wall if able
                        bugNext = getNext(bugXY,wallPoint,getAlignDist(false));
                        stepCheck(problem);//hit set to false at start of stepCheck
                        if(!hit){
                            //checkMLine(problem.q_init, problem.q_goal);
                            path.waypoints.push_back(bugNext);
                            bugXY = bugNext;
                            hit = true;
                        }
                    }
                    // Check if bug has hit original hit point
                    /*if(followMode && (bugXY - exitPoint).norm() < 1.05*bugStep && stepsAway > 2){
                        // Update point of loop exit to leave point
                        std::cout << "Never found closer point on M-Line, solution not possible :((((" << std::endl;
                        std::cout << "Leave point score:  " << qLScore << " and leavePoint: "<< qLeave << std::endl;
                        path.waypoints.push_back(problem.q_goal);
                        followMode = false;
                    }*/
                }while(followMode);
                // reset leave point score
                qLScore = 0;
                return path;
        }
    private:
        // Add any member variables here...
};