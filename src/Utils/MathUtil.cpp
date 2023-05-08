//
// Created by Zeyu Wang on 2021/8/7.
//

#include "../../include/Utils/MathUtil.h"
#include <cmath>
#include <cassert>

//标准差σ=sqrt(s^2)
double MathUtil::deviation(const float* ts, int start, int end, double dAve) {
    int len = end - start;
    double dVar=0;
    for(int i = start; i < end; ++i)
        dVar += (ts[i] - dAve) * (ts[i] - dAve);
    return sqrt(dVar / len);
}

double MathUtil::avg(const double* timeSeries, int start, int end) {
    double sum = 0;
    for(int i=start;i<end;++i)
        sum += timeSeries[i];
    return sum / (end - start);
}

double MathUtil::deviation(double* timeSeries, int start, int end) {
    double mean = avg(timeSeries, start, end), standardDeviation = 0;

    for(int i=start; i<end; ++i)
        standardDeviation += ((timeSeries[i] - mean) * (timeSeries[i] - mean));

    return sqrt(standardDeviation/(end - start));
}

double MathUtil::deviation(double* timeSeries, int len) {
    return deviation(timeSeries, 0, len);
}

int MathUtil::bitCount(int input){
    int res = 0;
    while (input!=0) {
        res += (input % 2);
        input /= 2;
    }
    return res;
}

// k must be less than 17, n must less than 31
int MathUtil::nChooseK(int n, int k) {
    k = min(k, (n - k));
    if (k <= 1) {  // C(n, 0) = 1, C(n, 1) = n
        return k == 0 ? 1 : n;
    }
    int limit = numeric_limits<int>::max() >> (31 - n);
    int cnk = 0;
    for (int i = 3; i < limit; i++) {
        if (bitCount(i) == k) {
            cnk++;
        }
    }
    return cnk;
}

//将points分割成两段，生成两倍长的数组，要求每段长度至少为minLength，否则不分割
int MathUtil::split(const vector<int> &points, int minLength, int len, vector<int> &ret) {
    auto* newPoints = new int[len + len];
    int c = 0;
    for (int i = 0; i < len; i++) {
        int length = points[i]; //id==0
        if (i > 0)
            length = (int) (points[i] - points[i - 1]);

        if (length >= minLength * 2) {
            int start = 0;
            if (i > 0)
                start = points[i - 1];
            newPoints[c++] = (int) (start + length / 2);
        }
        newPoints[c++] = points[i];
    }
    ret.resize(c);
    copy(newPoints, newPoints + c, ret.begin());
    return c;
}

void MathUtil::get1BitChangedNums(int x, int n, int *array, int start)
{
    int temp = 1;
    array[start++] = temp ^ x;
    for(int i=2; i<=n; ++i)
    {
        temp = temp << 1;
        array[start++] = temp ^ x;
    }
}

void MathUtil::get2BitsChangedNums(int x, int n, int *array, int start)
{
    assert(n >= 2);
    int temp, mask;
    for(int i = 1; i < n; i++)
    {
        mask = 1 << i;
        temp = 1 << (i - 1);
        for(int j = i + 1; j <= n; ++j, mask <<= 1)
        {
            int realMask = mask + temp;
            array[start++] = x ^ realMask;
        }
    }
}

void MathUtil::get3BitsChangedNums(int x, int n, int *array, int start)
{
    assert(n>=3);
    int temp, mask;
    for(int i = 1; i < n - 1; i++)
    {
        for(int j = i + 1; j < n; ++j)
        {
            temp = (1 << (i - 1)) + (1 << (j - 1));
            mask = 1 << j;
            for(int k = j + 1; k <= n; ++k, mask <<= 1)
            {
                int realMask = mask + temp;
                array[start++] = x ^ realMask;
            }
        }
    }
}

void MathUtil::get4BitsChangedNums(int x, int n, int *array, int start)
{
    assert(n >= 4);
    int temp, mask;
    for(int i = 1; i < n - 2; i++)
    {
        for(int j = i + 1; j < n - 1; ++j)
        {
            for(int k = j + 1; k < n; ++k)
            {
                temp = (1 << (i - 1)) + (1 << (j - 1)) + (1 << (k - 1));
                mask = 1 << k;
                for(int m = k + 1; m <= n; ++m, mask <<= 1)
                {
                    int realMask = mask + temp;
                    array[start++] = x ^ realMask;
                }
            }
        }
    }
}

int MathUtil::bitDiffNum(int i, int j, int n){
    int differentBits = i ^ j;
    int res = 0;
    for(int k=0; k<n; ++k){
        if(differentBits == 0)
            break;
        int tmp = differentBits >> 1;
        if(tmp + tmp + 1 == differentBits)
            res++;
        differentBits = tmp;
    }
    return res;
}

double MathUtil::errorRatio(vector<PqItemSeries*>&approx, vector<PqItemSeries*>&exact, int k){
    double sum = 0;
    for(int i=0;i<approx.size();++i){
        if(approx[i]->dist == 0 && exact[i]->dist == 0) sum +=1;
        else sum += (sqrt(approx[i]->dist) / sqrt(exact[i]->dist));
    }
    return sum / k;
}

double MathUtil::tightness(vector<PqItemSeries *> &approx, vector<PqItemSeries *> &exact, int k) {
    double sum = 0;
    for(int i=0;i<approx.size();++i){
        if(approx[i]->dist == 0 && exact[i]->dist == 0) sum +=1;
        else sum += (sqrt(exact[i]->dist) / sqrt(approx[i]->dist));
    }
    return sum / k;
}

int * MathUtil::generateMask(int segno){
    int *mask = new int[segno];
    mask[0] = 1 << (segno - 1);
    for(int i=1;i<segno;++i)
        mask[i] = mask[i-1] >> 1;
    return mask;
}

int MathUtil::generateMaskSettingKbits(const int* pos, int k, int len) {
    int res = 0, cur = 0;
    for(int i = 0;i<len;++i){
        res <<= 1;
        if(cur < k && pos[cur] == i){
            res += 1;
            ++cur;
        }
    }
    return res;
}