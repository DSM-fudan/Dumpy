//
// Created by Zeyu Wang on 2021/8/25.
//

#ifndef DUMPY_CONST_H
#define DUMPY_CONST_H
#include <string>
#include <iostream>
#include <sys/time.h>
#include "../include/Utils/INIReader.h"

// Mihalis Tsoukalos and Nikos Platis
#include <limits>


using namespace std;
class Const {

public:
   // sec:expr
    static string dataset;
    static int index, ops, query_num, series_num, k, visited_node_num, dtw_window_size, nprobes;
    static double search_ratio, dtw_window_percent;

    //sec: parameter
    // segment number is special, it needs to be specified here
    const static int segmentNum = 16;
    static int th, bitsCardinality, delta, r_max, fbl_size, max_diff, fbl_series_num, pre_read;
    static double fuzzy_f_1st, fuzzy_f, alpha_1st, alpha, small_perc, max_mask_bit_percentage, f_low, f_high;

    //sec: others
    static string graphfn, breakpointsfn;
    static int bitsReserve;

    //sec: dataset
    static string paafn, saxfn, idxfn, fuzzyidxfn, memoryidxfn, datafn, queryfn;
    static int maxK, tsLength;

    // 2-nd level parameter

    static int tsLengthPerSegment, cardinality, tsLengthBytes, vertexNum;
    static long offset;

    static void readConfig(){
        INIReader reader("../config.ini");

        if (reader.ParseError() < 0) {
            cout << "Can't load '.ini'\n";
            return;
        }
        dataset = reader.Get("expr", "dataset","");
        cout << "dataset: " << dataset<<endl;
        if(dataset.empty())   exit(-1);

        index = reader.GetInteger("expr", "index",-1);
        cout << "index: " << index<<endl;

        ops = reader.GetInteger("expr", "ops",-1);
        cout << "ops: " << ops<<endl;

        query_num = reader.GetInteger("expr", "query_num", -1);
        cout << "query_num: " << query_num << endl;

        series_num = reader.GetInteger("expr", "series_num", -1);
        cout << "series_num: " << series_num << endl;

        k = reader.GetInteger("expr", "k", -1);
        cout << "k: " << k << endl;

        dtw_window_percent = reader.GetReal("expr", "dtw_window_percent", -1);
        cout << "dtw_window_percent: " << dtw_window_percent << endl;

        visited_node_num = reader.GetInteger("expr", "visit_node_num", -1);
        cout << "visited_node_num: " << visited_node_num << endl;

        nprobes = reader.GetInteger("expr", "nprobes", -1);
        cout << "nprobes: " << nprobes << endl;

        th = reader.GetInteger("parameter", "th", 10000);
        cout << "th: " << th << endl;

        fbl_size = reader.GetInteger("parameter", "fbl_size", 1024 * 4);
        cout << "fbl_size: " << fbl_size << endl;

        max_diff = reader.GetInteger("parameter", "max_diff", 2);
        cout << "max_diff: " << max_diff << endl;

        bitsCardinality = reader.GetInteger("parameter", "bitsCardinality", 8);
        cout << "bitsCardinality: " << bitsCardinality << endl;

        fuzzy_f_1st = reader.GetReal("parameter", "fuzzy_f_1st", 0.2);
        cout << "boundary_1st: " << fuzzy_f_1st << endl;

        fuzzy_f = reader.GetReal("parameter", "fuzzy_f", 0.3);
        cout << "boundary: " << fuzzy_f << endl;

        delta = reader.GetInteger("parameter", "delta", 10000);
        cout << "max_replica: " << delta << endl;

        alpha = reader.GetReal("parameter", "alpha", 0.5);
        cout << "weighting factor: " << alpha << endl;

        f_low = reader.GetReal("parameter", "f_low", 0.5);
        cout << "f_low: " << f_low << endl;

        f_high = reader.GetReal("parameter", "f_high", 0.5);
        cout << "weighting factor: " << f_high << endl;

        max_mask_bit_percentage = reader.GetReal("parameter", "max_mask_bit_percentage", 0.5);
        cout << "filling_factor: " << max_mask_bit_percentage << endl;

        small_perc = reader.GetReal("parameter", "small_perc", 0.5);
        cout << "filling_factor: " << small_perc << endl;

        r_max = reader.GetInteger("parameter", "max_radius", 6);
        cout << "max_radius: " << r_max << endl;

        graphfn = reader.Get("other", "graphfn","");
        cout << "graphfn: " << graphfn <<endl;

        breakpointsfn = reader.Get("other", "breakpointsfn","");
        cout << "breakpointsfn: " << breakpointsfn <<endl;

        bitsReserve = reader.GetInteger("other", "bitsReserve", 6);
        cout << "bitsReserve: " << bitsReserve << endl;

        tsLength = reader.GetInteger(dataset, "tsLength", -1);
        cout << "tsLength: " << tsLength << endl;
        if(tsLength == -1)  exit(-1);

        maxK = reader.GetInteger(dataset, "maxK", -1);
        cout << "maxK: " << maxK << endl;

        paafn = reader.Get(dataset, "paafn","");
        cout << "paafn: " << paafn <<endl;

        saxfn = reader.Get(dataset, "saxfn","");
        cout << "saxfn: " << saxfn <<endl;

        idxfn = reader.Get(dataset, "idxfn","");
        cout << "idxfn: " << idxfn <<endl;

        fuzzyidxfn = reader.Get(dataset, "fuzzyidxfn","");
        cout << "fuzzyidxfn: " << fuzzyidxfn <<endl;

        memoryidxfn = reader.Get(dataset, "memoryidxfn","");
        cout << "memoryidxfn: " << memoryidxfn <<endl;

        datafn = reader.Get(dataset, "datafn","");
        cout << "datafn: " << datafn <<endl;

        queryfn = reader.Get(dataset, "queryfn","");
        cout << "queryfn: " << queryfn <<endl;

        tsLengthPerSegment = tsLength / segmentNum;
        cardinality = 1<<bitsCardinality;
        dtw_window_size = dtw_window_percent * Const::tsLength;
        tsLengthBytes = tsLength * 4;
        offset = ((long )(cardinality - 1) * (cardinality - 2)) / 2;
        vertexNum = 1 << segmentNum;
        fbl_series_num = (fbl_size/ tsLengthBytes) * 1024 * 1024 ;
        int max_series_num = numeric_limits<unsigned>::max() / Const::tsLength;
        fbl_series_num = fbl_series_num > max_series_num ? max_series_num: fbl_series_num;

        cout << "buffer_series_num: " << fbl_series_num << endl;
    }

    static string now() {
        time_t t = time(nullptr);
        char buffer[9] = {0};
        strftime(buffer, 9, "%H:%M:%S", localtime(&t));
        return (buffer);
    }

    static void logPrint(const string& content){
        cout << now() << ": " << content <<endl;
    }

    static void timer_start(struct timeval *timer)
    {
        (void)gettimeofday(timer, nullptr);
    }

    static double timer_end(const struct timeval *start)
    {
        struct timeval now{};
        (void)gettimeofday(&now, NULL);
        return (now.tv_sec - start->tv_sec)*1000000 + ((double)now.tv_usec - start->tv_usec) ;
    }

};


#endif //DUMPY_CONST_H
