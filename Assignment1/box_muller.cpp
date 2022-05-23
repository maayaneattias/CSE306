#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>    // std::min


static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0.,1.);
static std::normal_distribution<double> normal(0.,1.);

void boxMuller(double stdev, double &x, double &y){
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    
    x = sqrt(-2 * log(r1)) * cos(2*M_PI*r2)*stdev;
    y = sqrt(-2 * log(r1)) * sin(2*M_PI*r2)*stdev;

}

double gaussian(double stdev, double x){
    double param = 1 / (stdev*sqrt(2*M_PI));
    double param2 = exp(-(pow(x,2.))/(2*pow(stdev,2.)));
    return param * param2;
}

int main(){
    double x,y,z,w,p;
    double stdev = 1;
    int N = 100000;
    double total  = 0;
    int i = 0;
    for (int i=0;i<N;i++){
        boxMuller(stdev,x,y); 
        boxMuller(stdev,z,w); 
        if(x< -M_PI/2 | x > M_PI/2) continue;
        if(y< -M_PI/2 | y > M_PI/2) continue;
        if(z< -M_PI/2 | z > M_PI/2) continue;
        p = gaussian(stdev,x) *gaussian(stdev,y)  *gaussian(stdev,z); 
        total += cos(x*y*z)/p;
    }
    total /= N;
    std::cout << total << std::endl;
}