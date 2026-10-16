// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#include "../include/variables.hpp"
#include "../include/estimators.hpp"
#include <iostream>
#include <iomanip>

int main(){
    RNG rng(938102);std::vector<int> ranks{2,7,15};
    double sum[3]={},product[3][3]={};const int samples=100000;
    bool ordered=true;
    for(int i=0;i<samples;++i){
        auto times=sse_selected_event_times(23,ranks,rng);
        ordered=ordered && times[0]>0 && times[0]<times[1] && times[1]<times[2] && times[2]<1;
        for(int a=0;a<3;++a){sum[a]+=times[a];for(int b=0;b<3;++b)product[a][b]+=times[a]*times[b];}
    }
    std::cout<<std::setprecision(17)<<"{\"ordered\":"<<(ordered?"true":"false")<<",\"samples\":"<<samples<<",\"means\":[";
    for(int a=0;a<3;++a)std::cout<<(a?",":"")<<sum[a]/samples;
    std::cout<<"],\"covariance\":[";
    for(int a=0;a<3;++a){std::cout<<(a?",":"")<<"[";for(int b=0;b<3;++b)std::cout<<(b?",":"")<<(product[a][b]/samples-sum[a]*sum[b]/(double(samples)*samples));std::cout<<"]";}
    std::cout<<"]}\n";
}
