//
// Created by Zeyu Wang on 2021/8/9.
//
#include <cassert>
#include "../include/DataStructures/PqItemSeries.h"
#include "../include/Utils/TimeSeriesUtil.h"
#include "../include/Const.h"

PqItemSeries::PqItemSeries(float *t, float *q, bool free, bool ndc) {
    needDeepCopy = ndc;
    needFree = free;
    ts = t;
    dist = TimeSeriesUtil::euclideanDist(ts, q, Const::tsLength);

}

void PqItemSeries::copyData() {
    if(!needDeepCopy)   return;
    assert(ts != nullptr);
    auto *tmp = new float[Const::tsLength];
    copy(ts, ts + Const::tsLength, tmp);
    ts = tmp;
    needFree = true;
    needDeepCopy = false;
}
