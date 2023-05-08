//
// Created by Zeyu Wang on 2021/12/12.
//
#include "../include/Const.h"

// sec:expr
string Const::dataset = "";
int Const::index = -1, Const::ops = -1, Const::query_num = -1, Const::series_num = -1, Const::k = -1,
Const::visited_node_num = -1, Const::dtw_window_size = -1, Const::nprobes=-1;
double Const::dtw_window_percent = -1;

//sec: parameter
// segment number is special, it needs to be specified here
int Const::th = -1, Const::bitsCardinality = -1, Const::delta = -1, Const::r_max = -1, Const::fbl_size = -1,
Const::max_diff = -1, Const::fbl_series_num = -1, Const::pre_read = 1;
double Const::fuzzy_f_1st = -1, Const::fuzzy_f = -1, Const::alpha = -1, Const::small_perc = -1, Const::max_mask_bit_percentage=-1, Const::f_high = -1, Const::f_low=-1;

//sec: others
string Const::graphfn = "", Const::breakpointsfn = "";
int Const::bitsReserve = -1;

//sec: dataset
string Const::paafn = "", Const::saxfn = "", Const::idxfn = "", Const::memoryidxfn = "",
Const::fuzzyidxfn, Const::datafn = "", Const::queryfn = "";
int Const::tsLength = -1, Const::maxK = 0;
// 2-nd level parameter

int Const::tsLengthPerSegment = -1, Const::cardinality = -1, Const::tsLengthBytes = -1, Const::vertexNum = -1;
long Const::offset = -1;