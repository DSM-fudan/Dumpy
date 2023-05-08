//
// Created by Zeyu Wang on 2021/8/7.
//

#ifndef DUMPY_TIMESERIESUTIL_H
#define DUMPY_TIMESERIESUTIL_H
#include "../DataStructures/PqItemSeries.h"
#include "../Const.h"
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <bitset>


class TimeSeriesUtil {
public:
    static string timeSeries2Line(float* timeSeries) {
        string sb;
        for (int i=0; i < Const::tsLength; ++i) {
            sb += formatFloatValue(timeSeries[i], 4) + " , ";
        }
        sb.append("\n");
        return sb;
    }

    static string timeSeries2Line(vector<float>* timeSeries) {
        string sb;
        for (int i=0; i < Const::tsLength; ++i) {
            sb += formatFloatValue((*timeSeries)[i], 4) + " , ";
        }
        sb.append("\n");
        return sb;
    }

    static string formatFloatValue(float val, int fixed) {
        ostringstream oss;
        oss << setprecision(fixed) << val;
        return oss.str();
    }

    static double euclideanDist(const float* ts_1, const float* ts_2, int len);

    /**
     * using avg and std
     *
     * @param node
     * @param queryTs
     * @return
     * @throws IOException
     */

    template<typename ... Args>
    static string str_format(const string &format, Args... args);

    static void heap_data_copy(vector<PqItemSeries *> &heap);

    static int intersectionTsSets(const vector<PqItemSeries *> *tsSet1, vector<float *> *tsSet2);

    static double euclideanDist(const float *ts_1, const float *ts_2, int len, double bound);

    static int intersectionTsSetsCardinality(const vector<PqItemSeries *> &tsSet1, const vector<PqItemSeries *> &tsSet2);

    static  double dtw(const float *A, const float *B, int len, int r, double bsf);

    static double dtwsimd(const float* A, const float* B, float* cb, int len, int r, float bsf,
                                          float* tSum, float* pCost, float* rDist);

    static bool isSame(const PqItemSeries *ts1, const float *ts2);

    static float lb_keogh_data_bound(float* qo, float* upperLemire, float* lowerLemire, float* cb, int len, float bsf);


};


#endif //DUMPY_TIMESERIESUTIL_H
