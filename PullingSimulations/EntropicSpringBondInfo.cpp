/**
 * Copyright (C) 2022 Lehigh University.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 * Author: David Rutkowski (dmr518@lehigh.edu)
 */

#include <iostream>
#include "EntropicSpringBondInfo.h"

using namespace std;

EntropicSpringBondInfo::EntropicSpringBondInfo()
{
    kT = 300.0;
    N = 1;
    b = 1.0;
}

EntropicSpringBondInfo::EntropicSpringBondInfo(double newkT, int newN, double newB)
{
    kT = newkT;
    N = newN;
    b = newB;
}

double EntropicSpringBondInfo::calcForce(double dist)
{
    /*if(dist <= eqDist)
    {
        //entropic spring (rather then expansive force in this region it is a compressive force)
        return -3.0 * kT / (double)N / (b*b) * dist;
    }
    else
    {
        // normal spring when every fictious segment is stretched to maximum extent added to entropic spring force at eqDist
        return -k * (dist - eqDist) -3.0 * kT / (double)N / (b*b) * eqDist;
    }*/
    
    double parenthesisValue = dist/(b*N);
    double invLValue = invLangevin(parenthesisValue);
    
    return -kT/b*invLValue;
}

double EntropicSpringBondInfo::invLangevin(double x)
{
    double x2 = x*x;
    double x4 = x2*x2;
    double x6 = x2*x4;
    double x8 = x4*x4;
    
    double val = x*(3.0 - 1.00651*x2 - 0.962251*x4 + 1.4353*x6 - 0.48953*x8) / (1.0-x) / (1.0+1.01524*x);
    
    return val;
}