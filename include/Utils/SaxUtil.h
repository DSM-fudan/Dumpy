//
// Created by Zeyu Wang on 2021/8/7.
//

#ifndef DUMPY_SAXUTIL_H
#define DUMPY_SAXUTIL_H
#include <string>
#include <vector>

using namespace std;

typedef struct dequeue
{   int *dq;
    int size,capacity;
    int f,r;
} dequeue;

class SaxUtil {
public:
    static double* breakpoints;

    static double *bp8;

    static double* readFromstring(string str);

    static double *readDoubleFromFileAtOnce();

    static float * paaFromTsFloat(const float* ts, int tsLengthPerSegment, int segmentNum);

    static vector<int> * saxFromTs(float*ts, int tsLengthPerSegment, int segmentNum, int cardinality);

    static int findFirstGE(const double* array, int start, int length, double target);// satisfy condition: array[?] >= target  and the first one

    /**
     This function prints a sax record.
     */
    static void saxPrint(int* sax, int bits_cardinality, int segment_num);

    static void printBinary(long n, int size);

    static vector<unsigned short> * saxFromPaa(float *paa, int segmentNum, int cardinality);

    static int invSaxHeadFromSax(vector<int> *sax, int bitsCardinality, int segmentNum);

    static void saxFromTs(const float *ts, unsigned short *sax, int tsLengthPerSegment, int segmentNum, int cardinality);

    static double LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, int bits_cardinality);


    static void generateSaxFile(const string &fn, const string &output);

    static void id2Sax2(int id, unsigned short *sax, int segment_num);

    static int invSaxHeadkFromSax(const unsigned short *sax, int bitsCardinality, int segmentNum, int k);

    static void paaFromTs(const float *ts, float *paa, int tsLengthPerSegment, int segmentNum);

    static void generatePaaFile(const string &fn, const string &output);

    static int findFirstGE(const int *array, int start, int length, int target);

    static void getValueRange(int sax_single, int bits_cardinality, double *lb, double *ub);

    static int extendSax(float *paa, const int *bits_cardinality, vector<int> &segments);

    static int invSaxHeadFromPaa(const float *paa, int tsLengthPerSegment, int segmentNum);

    static int extendSax(const unsigned short *sax, const int *bits_cardinality);

    static int getNewId(const float *paa, const float *split_line);

    static int getNewId(const float *paa, const float *split_line, vector<int> &segments);

    static int extendSax(const unsigned short *sax, const int *bits_cardinality, vector<int> &segments);

    static int invSaxHeadFromSax(const unsigned short *sax, int bitsCardinality, int segmentNum);

    static double
    LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, const int *bits_cardinality,
                        vector<int> &chosen_segs,
                        int new_id);

    static double LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, const int *bits_cardinality);

    static int extendSax(const unsigned short *sax, const int *bits_cardinality, vector<int> &segments,
                  const unsigned short *parent_sax);

    static void paaAndSaxFromTs(const float *ts, float *paa, unsigned short *sax, int tsLengthPerSegment, int segmentNum,
                         int cardinality);

    static double LowerBound_Paa_iSax(const float *paa, const unsigned short *sax);

    static double getMidLineFromSaxSymbolbc8(unsigned short symbol);

    static void lower_upper_lemire(const float *t, int len, int r, float *l, float *u);

    static double *paaFromTs(const float *ts, int tsLengthPerSegment, int segmentNum);

    static double minidist_paa_to_isax_DTW(const double *paaU, const double *paaL, const unsigned short *sax,
                                    const int *sax_cardinalities);

    static double minidist_paa_to_isax_DTW(const double *paaU, const double *paaL, const unsigned short *sax,
                                                    const int* bits_cardinality, vector<int>&chosen_segs, int new_id);

    static double getMinDist1stLayer(const float *paa, int id);

    static double getMinDist1stLayerDTW(const double *paaU, const double *paaL, int id);
};



#endif //DUMPY_SAXUTIL_H
