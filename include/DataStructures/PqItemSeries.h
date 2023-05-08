//
// Created by Zeyu Wang on 2021/8/7.
//

#ifndef DUMPY_PQITEMSERIES_H
#define DUMPY_PQITEMSERIES_H

class PqItemSeries {
public:
    float* ts{};
    double dist{};
    bool needFree= false, needDeepCopy = false;

    PqItemSeries()= default;

    PqItemSeries(float*t, float*q, bool free=false, bool ndc = false);

    PqItemSeries(float*t, double d, bool free= false, bool ndc = false){
        ts = t;
        dist = d;
        needDeepCopy = ndc;
        needFree = free;
    }

    ~PqItemSeries(){
        if(needFree)
            delete[] ts;
    }

    void copyData();

};

struct PqItemSeriesMaxHeap{
    bool operator()(const PqItemSeries* x, const PqItemSeries* y){
        return x->dist < y->dist;
    }
};

#endif //DUMPY_PQITEMSERIES_H
