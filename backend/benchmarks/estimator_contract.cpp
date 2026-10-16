// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#include "../include/estimators.hpp"
#include <iostream>
#include <iomanip>

int main(){
    double worst=0.;
    for(int n=1;n<=37;++n){
        std::vector<std::vector<double>> states(n,std::vector<double>(2));
        for(int p=0;p<n;++p){states[p][0]=((p/3)%2?-.5:.5);states[p][1]=((p/5)%3-1)*.25;}
        std::vector<SSEPathSegment> path;
        for(int p=0;p<n;++p){
            if(!path.empty() && path.back().values==states[p])path.back().end=p+1;
            else path.push_back({p,p+1,states[p]});
        }
        auto estimated=sse_path_correlation(path,n);
        double fractions[5]={0,.125,.25,.375,.5};
        for(int t=0;t<5;++t){
            double f=fractions[t],brute=0.,mass=std::pow(1-f,n);
            for(int m=0;m<=n;++m){
                double cyclic=0.;
                for(int p=0;p<n;++p)for(int c=0;c<2;++c)cyclic+=states[p][c]*states[(p+m)%n][c]/(2*n);
                brute+=mass*cyclic;
                if(m<n)mass*=double(n-m)*f/((m+1)*(1-f));
            }
            worst=std::max(worst,std::abs(brute-estimated[t]));
        }
    }
    std::cout<<std::setprecision(17)<<worst<<'\n';
    return worst<2e-12?0:1;
}
