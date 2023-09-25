//
// Created by Zeyu Wang on 2021/8/7.
//

#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/FileUtil.h"
#include <bitset>
#include <immintrin.h>
#define dist(x,y) (((x)-(y)*((x)-(y))))


// Mihalis Tsoukalos and Nikos Platis
#include <memory>

using namespace std;

void TimeSeriesUtil::heap_data_copy(vector<PqItemSeries *> &heap){
    for(auto *pis:heap)
        pis->copyData();
}

bool TimeSeriesUtil::isSame(const PqItemSeries *ts1, const float *ts2)
{

    float *t1 = (ts1->ts);
    for(int i = 0; i < Const::tsLength; ++i) {
        if(abs(t1[i] - ts2[i])>1e-5)
            return false;
    }
    return true;;
}

int TimeSeriesUtil::intersectionTsSets(const vector<PqItemSeries *> *tsSet1, vector<float *> *tsSet2){
    int intersectionNum = 0, size= tsSet1->size();
    for(const PqItemSeries* currentTs : *tsSet1)
        for(int i=0;i<size;++i){
            if(isSame(currentTs, (*tsSet2)[i])) {
                intersectionNum++;
                break;
            }
        }

    return intersectionNum;
}

int TimeSeriesUtil::intersectionTsSetsCardinality(const vector<PqItemSeries *> &tsSet1, const vector<PqItemSeries *> &tsSet2){
    int intersectionNum = 0;
    for(PqItemSeries* currentTs : tsSet1)
        for (PqItemSeries* targetTs : tsSet2)
            if (isSame(currentTs, targetTs->ts)) {
                intersectionNum += 1;
                break;
            }
    return intersectionNum;
}

float TimeSeriesUtil::lb_keogh_data_bound(float* qo, float* upperLemire, float* lowerLemire, float* cb, int len,
                                          float bsf)
{
    float lb = 0;
    float uu=0,ll=0,d=0;
    int i=0;

    int len1 = (len/8)*8;
    __m256 tu256, tl256, cb256, Q, calc1, calc2;
    __m128 temp1, temp2;
    float *cbtmp = static_cast<float *>(malloc(sizeof(float) * 8));

    for(i=0; i<len1&&lb<bsf; i+=8)
    {
        Q = _mm256_loadu_ps(&qo[i]);
        tu256 = _mm256_loadu_ps(&upperLemire[i]);
        tl256 = _mm256_loadu_ps(&lowerLemire[i]);
        //tu256 = _mm_setr_ps(tu[order[i]],tu[order[i+1]],tu[order[i+2]],tu[order[i+3]]);
        //tl256 = _mm_setr_ps(tl[order[i]],tl[order[i+1]],tl[order[i+2]],tl[order[i+3]]);
        calc1 = _mm256_min_ps(Q,tu256);
        calc1 = _mm256_sub_ps(Q,calc1);

        calc2 = _mm256_max_ps(Q,tl256);
        calc2 = _mm256_sub_ps(calc2,Q);
        calc1 = _mm256_add_ps(calc1,calc2);

        calc1 = _mm256_mul_ps(calc1,calc1);

        _mm256_storeu_ps(cbtmp,calc1);

        calc1 = _mm256_hadd_ps(calc1,calc1);
        calc1 = _mm256_hadd_ps(calc1,calc1);
        temp1 = _mm256_extractf128_ps(calc1,1);
        temp2 = _mm_add_ss(_mm256_castps256_ps128(calc1),temp1);
        lb += _mm_cvtss_f32(temp2);

        cb[i]=cbtmp[0];cb[i+1] = cbtmp[1];cb[i+2]=cbtmp[2];cb[i+3]=cbtmp[3];
        cb[i+4]=cbtmp[4];cb[i+5] = cbtmp[5];cb[i+6]=cbtmp[6];cb[i+7]=cbtmp[7];
    }

    for(;i<len&&lb<bsf;i++)
    {
        uu = upperLemire[i];
        ll = lowerLemire[i];
        d= 0;
        if(qo[i] > uu)
        {
            d = dist(qo[i],uu);
        }
        else if(qo[i] < ll)
        {
            d = dist(qo[i],ll);
        }
        lb += d;
        cb[i] = d;
    }

    free( cbtmp);
    return lb;
}

double TimeSeriesUtil::euclideanDist(const float* ts_1, const float* ts_2, int len) {
    double sum = 0, dp;
    for (int i = 0; i < len; i++) {
        dp = ts_1[i] - ts_2[i];
        sum += dp * dp;
    }
    return sum;
}

double TimeSeriesUtil::dtw(const float* A, const float* B, int len, int r, double bsf)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=numeric_limits<double>::max();

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=numeric_limits<double>::max();

    for (i=0; i<len; i++)
    {
        k = max(0,r-i);
        min_cost = numeric_limits<double>::max();

        for(j=max(0,i-r); j<=min(len-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=dist(A[0], B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = numeric_limits<double>::max();
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = numeric_limits<double>::max();
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = numeric_limits<double>::max();
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(A[i], B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (i+r < len-1 && min_cost  >= bsf)
        {   free(cost);
            free(cost_prev);
            return min_cost ;
        }
        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}

double TimeSeriesUtil::dtwsimd(const float* A, const float* B, float* cb, int len, int r, float bsf, float* tSum,
                               float* pCost, float* rDist)
{
    int length = 2*r + 1;
    // SIMD register
    //__m256 a256, b256;

    int start, end;
    float minCost = 0.0f;
    // the first line
    for(int k=0; k<=r; k++)
    {
        rDist[k] = dist(A[0], B[k]);
    }

    tSum[0] = rDist[0];
    for(int ij=1; ij<=r; ij++)
        tSum[ij] = tSum[ij-1] + rDist[ij];


    pCost[0] = tSum[0];
    for(int ij= 1; ij<=r; ij++)
    {
        pCost[ij] = min(tSum[ij-1],tSum[ij]);
    }
    pCost[r+1] = tSum[r];

    for(int i=1; i < len - 1; i++)
    {
        start = max(0,i-r);
        end = min(len - 1, i + r);

        for(int k=start; k<=end; k++)
        {
            rDist[k-start] = dist(A[i], B[k]);
        }

        for(int k=start; k<=end; k++)
        {
            tSum[k-start] =  pCost[k-start] + rDist[k-start];
        }

        minCost = tSum[0];
        for(int k=start+1; k<=end; k++)
        {
            if(tSum[k-1-start]<pCost[k-start])
            {
                tSum[k-start] = tSum[k-1-start] + rDist[k-start];
            }
            minCost = min(minCost, tSum[k - start]);
        }
        if(i+r < len - 1 && minCost + cb[i + r + 1] >= bsf)
        {
            return minCost + cb[i+r+1];
        }

        if((end-start+1) < length && start == 0)
        {
            pCost[start-start] = tSum[start-start];
            for(int ij= start+1; ij<=end; ij++)
            {
                pCost[ij-start] = min(tSum[ij-1-start],tSum[ij-start]);
            }
            pCost[end+1-start] = tSum[end-start];
        }
        else
        {
            for(int ij= start+1; ij<=end; ij++)
            {
                pCost[ij-1-start] = min(tSum[ij-1-start],tSum[ij-start]);
            }
            pCost[end-start] = tSum[end-start];
        }
    }

    // the last line
    start = len - 1 - r;
    end = len - 1;

    for(int k=start; k<=end; k++)
    {
        rDist[k-start] = (A[len - 1] - B[k]) * (A[len - 1] - B[k]);
    }

    for(int k=start; k<=end; k++)
    {
        tSum[k-start] =  pCost[k-start] + rDist[k-start];
    }
    for(int k=start+1; k<=end; k++)
    {
        if(tSum[k-1-start]<pCost[k-start])
        {
            tSum[k-start] = tSum[k-1-start] + rDist[k-start];
        }
    }
    double ret = tSum[r];
    return ret;
}

double TimeSeriesUtil::euclideanDist(const float* ts_1, const float* ts_2, int len, double bound) {
    double sum = 0, dp;
    for (int i = 0; i < len && sum < bound; i++) {
        dp = ts_1[i] - ts_2[i];
        sum += dp * dp;
    }
    return sum;
}

template<typename ... Args>
string TimeSeriesUtil::str_format(const string &format, Args ... args)
{
    auto size_buf = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1;
    std::unique_ptr<char[]> buf(new(std::nothrow) char[size_buf]);

    if (!buf)
        return string{};

    std::snprintf(buf.get(), size_buf, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size_buf - 1);
}
