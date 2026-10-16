// Copyright (c) 2026 Vikas Vijigiri. All rights reserved.
#ifndef SSE_ESTIMATORS_HPP
#define SSE_ESTIMATORS_HPP
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>
#include <stdexcept>

// Selected order statistics of n uniform operator times. Grouping their
// Dirichlet(1,...,1) spacings gives integer-shape gamma gaps, so only the
// off-diagonal event ranks need to be sampled for the boundary audit.
template<class Random>
std::vector<double> sse_selected_event_times(int n,const std::vector<int>& ranks,Random& rng){
    if(ranks.empty())return {};
    std::vector<double> gaps;gaps.reserve(ranks.size()+1);
    int previous=0;double total=0.;
    for(int rank:ranks){
        if(rank<=previous || rank>n)throw std::runtime_error("Invalid SSE event rank.");
        double value=rng.gamma(rank-previous);gaps.push_back(value);total+=value;previous=rank;
    }
    total+=rng.gamma(n+1-previous);
    double cumulative=0.;for(auto& value:gaps){cumulative+=value;value=cumulative/total;}
    return gaps;
}

// Exact conditional binomial average of cyclic propagation-index dipole
// correlations. Integrating it over tau gives Sandvik's Eq. (27).
// Piecewise constant paths avoid an O(n^2) correlation calculation.
struct SSEPathSegment {
    int begin, end;
    std::vector<double> values;
};

struct SSEBinomialTails {
    int n;
    double fraction;
    std::vector<double> probability, moment;
    SSEBinomialTails(int count,double f):n(count),fraction(f),probability(n+2),moment(n+2) {
        std::vector<double> mass(n+1);
        int mode=std::min(n,int((n+1)*f));
        mass[mode]=1.;
        for(int m=mode;m<n;++m)mass[m+1]=mass[m]*(n-m)*f/((m+1)*(1-f));
        for(int m=mode;m>0;--m)mass[m-1]=mass[m]*m*(1-f)/((n-m+1)*f);
        double total=0;for(double value:mass)total+=value;
        for(int m=n;m>=0;--m){probability[m]=probability[m+1]+mass[m]/total;moment[m]=moment[m+1]+m*mass[m]/total;}
    }
    double positive_part(int threshold) const {
        if(threshold<0)return n*fraction-threshold;
        if(threshold>=n)return 0.;
        return moment[threshold+1]-threshold*probability[threshold+1];
    }
    double interval_overlap(const SSEPathSegment& a,const SSEPathSegment& b) const {
        int length_a=a.end-a.begin,length_b=b.end-b.begin,offset=b.begin-a.begin;
        double overlap=0.;
        for(int winding=-1;winding<=1;++winding){
            int t=offset-winding*n;
            overlap+=positive_part(t-length_a)-positive_part(t)-positive_part(t-length_a+length_b)+positive_part(t+length_b);
        }
        return overlap;
    }
};

inline std::array<double,5> sse_path_correlation(const std::vector<SSEPathSegment>& path,int n) {
    std::array<double,5> result{};
    if(path.empty())return result;
    if(n==0){for(double v:path[0].values)result[0]+=v*v/path[0].values.size();result.fill(result[0]);return result;}
    for(const auto& segment:path)for(double v:segment.values)result[0]+=(segment.end-segment.begin)*v*v/(n*segment.values.size());
    // Cache depends only on expansion order, never on a spin configuration.
    static std::map<int,std::vector<SSEBinomialTails>> cache;
    auto found=cache.find(n);
    if(found==cache.end()){
        if(cache.size()>=256)cache.clear();
        std::vector<SSEBinomialTails> tails;
        for(double f:{.125,.25,.375,.5})tails.emplace_back(n,f);
        found=cache.emplace(n,std::move(tails)).first;
    }
    for(const auto& a:path)for(const auto& b:path){
        double dot=0.;for(size_t i=0;i<a.values.size();++i)dot+=a.values[i]*b.values[i]/a.values.size();
        for(int t=0;t<4;++t)result[t+1]+=dot*found->second[t].interval_overlap(a,b)/n;
    }
    return result;
}
#endif
