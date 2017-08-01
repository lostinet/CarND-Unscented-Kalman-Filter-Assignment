#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    // initialize rmse
    rmse << 0,0,0,0;
    
    // check validate input
    if (estimations.size() != ground_truth.size() || estimations.size() == 0)
    {
        std::cout << "invalid data" << std::endl;
        return rmse;
    }
    
    // culculate rmse thru interations
    for(int i = 0; i < estimations.size(); ++i)
    {
        VectorXd residual = estimations[i] - ground_truth[i];
        
        residual = residual.array() * residual.array();
        
        rmse += residual;
        
    }
    
    rmse = rmse/estimations.size();
    
    rmse = rmse.array().sqrt();
    
    return rmse;
    
}
