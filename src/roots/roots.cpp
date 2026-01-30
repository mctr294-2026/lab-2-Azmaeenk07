#include "roots.hpp"
#include <cmath> // needed for math tools, like abs()

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
 
    if (f(a) == 0) { *root = a; return true; } // if a is already a root
    if (f(b) == 0) { *root = b; return true; } // if b is already a root 

    // Step 1: ensure f(a) & f(b) have opposite signs
    if (f(a) * f(b) > 0){
        return false;
    }

    // Step 4: repeat process until f(c) is small enough (1e6 iteration limit)
    for (int count = 0; count < 1e6; count++){
                      
        double c = (a + b)/2; // Step 2: continually calulates the midpoint of a and b 

        // Step 3: if f(c) is small enough, it is the root. Otherwise replace a or b with c until c is small enough
        if (std::abs(f(c)) < 1e-6 || std::abs(b - a) < 1e-6) {
            *root = c;
            return true;
        }
        else if (f(a) * f(c) > 0) { 
            a = c; // moves a up to c
        }        
        else {
            b = c;
        }

    }

    return false;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {

    // Step 1: ensure f(a) & f(b) have opposite signs
    if (f(a) * f(b) >= 0){
        return false;
    }

    // Step 5: repeat process until f(c) is small enough (1e6 iteration limit)
    for (int count = 0; count < 1e6; count++){
                      
        double c = a - (f(a) * (b - a)) / (f(b) - f(a)); // Step 2: midpoint of [a,b] calculation

        // Step 3 & 4: if f(c) is small enough, it is the root. Otherwise replace a or b with c until c is small enough
        if (std::abs(f(c)) < 1e-6 || std::abs(b - a) < 1e-6) {
            *root = c;
            return true;
        }
        else if (f(a) * f(c) > 0) { 
            a = c; 
        }        
        else {
            b = c;
        }

    }

    return false;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
    
    // Step 1: start with the initial guess c (from the header args)
    double xn = c;

    // Step 3: repeat process until we find the answer 
    for (int count = 1; count < 1e6; count++) {

        // we need the function and its derivative value (g is the derivative)
        double fx = f(xn);
        double dfx = g(xn);

        // if derivative is zero, the program needs to fail; we can't divide by sero 
        if (std::abs(dfx) < 1e-6){
            return false;
        }

        // Step 2: formula for subsequent approximations
        double next_xn = xn - (fx/dfx);

        // header says fail if we leave the interval
        if (next_xn < a || next_xn > b) {
            return false;
        }

        // Step 3: we have the answer if the change is lower than our tolerance
        if (std::abs(next_xn-xn) < 1e-6) {
            *root = next_xn;
            return true;
        }

        xn = next_xn; // updates x to new value for every iteration until the root is found

    }

    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {

    // using a and b as our two starting points since secant needs 2
    double prev_xn = a;
    double xn = b;

    // Step 3: loop the process until the change becomes smaller than the tolerance
    for (int count = 0; count < 1e6; count++) {

        double fxn = f(xn);
        double fprev_xn = f(prev_xn);

        // need to stop the program from crashing when we divide by sero
        if (std::abs(fxn - fprev_xn) < 1e-6) {
            return false;
        }

        // Step 2: secant approximation formula
        double next_xn = xn - fxn * (xn - prev_xn) / (fxn - fprev_xn);

        // header says fail if we leave the interval
        if (next_xn < a || next_xn > b) {
            return false;
        }

        // Step 3: change between new and current guess is negligible; we have the root
        if (std::abs(next_xn - xn) < 1e-6) {
            *root = next_xn;
            return true;
        }

        prev_xn = xn; // old becomes current
        xn = next_xn; // current becomes new

    }

    return false;
}