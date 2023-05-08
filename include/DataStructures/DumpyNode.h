//
// Created by Zeyu Wang on 2021/11/8.
//

#ifndef DUMPY_DUMPYNODE_H
#define DUMPY_DUMPYNODE_H
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <unordered_set>
#include "../Const.h"
#include "TimeSeries.h"
#include "PqItemSeries.h"

using namespace std;

struct partUnit{
    int id;
    int size;
    int pid;

    static bool comp_size(partUnit*x, partUnit*y){
        return x->size > y->size;
    }
};

struct SAX_BUF_UNIT{
    int size{0};
    unsigned short *buffer{nullptr};
};

struct PAA_INFO{
    double paa_variance[Const::segmentNum]{};
    int paa_up_size[Const::segmentNum]{},
    paa_below_size[Const::segmentNum]{};
};

struct FBL_UNIT{
    int size{0};
    float *buffer{nullptr};
    int pos{0};
};

struct LBL_UNIT{
    vector<float*>buffer;
};

struct NODE_RECORDER;

class DumpyNode {
    DumpyNode(){;}
    explicit DumpyNode(int _layer) {
        layer = _layer;
    }
    DumpyNode(int _layer, int pid) {
        layer = _layer;
        partition_id = pid;
        bits_cardinality[0] = -1;
    }
    // only for node in 1st layer
    DumpyNode(int _layer, int _size, int _id) {
        layer = _layer;
        size = _size;
        id = _id;
        file_id = to_string(id);
        chosenSegments.clear();
        for(auto &i:sax)    i=0;
        for(auto &i:bits_cardinality)   i = 1;
    }
    DumpyNode(const DumpyNode* parent, int _size, int _id){
        layer = parent->layer + 1;
        id = _id;
        size = _size;
        file_id  = parent->file_id + "_" + to_string(id);
        chosenSegments.clear();
    }
    DumpyNode(const DumpyNode* parent, int pid){
        layer = parent->layer + 1;
        partition_id = pid;
        file_id  = parent->file_id + "_" + to_string(partition_id);
        bits_cardinality[0] = -1;
    }

    static void loadPaa(const string & paafn);
    static int loadSax(const string & saxfn);
    static long generateSaxAndPaaTbl();
    PAA_INFO* statPaa();
    void chooseSegment(PAA_INFO *paa, int chosen_num);
    static int chooseOneSegment(PAA_INFO *node);
    void generateSaxAndCardIn1stLayer(int new_id);
    void generateSaxAndCardinality(DumpyNode *node, int new_id);
    void generateSaxAndCardIn1stLayer4LeafNode(int new_id);
    void generateSaxAndCardinality4LeafNode(DumpyNode *node, int new_id);
    static int partition1stLayer(partUnit *nodes_map, vector<vector<int>> *g,double filling_factor);
    static int partition(partUnit *nodes_map);
    void growIndex();
    void growIndexFuzzy(unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl, vector<vector<int>> *g);

    void fuzzySeriesInPartUnit(partUnit *part_units, int actual_size, int chosen_num, vector<int> &node_offsets,
                               vector<int> &series_index_list,
                               unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl, int _id) const;
    void fuzzy(partUnit *part_units, vector<int> &actual_sizes, vector<vector<int>> &node_offsets,
               vector<vector<int>> &series_index_list, int chosen_num,
               unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl) const;
    void fuzzySeriesInPartUnitInFirstLayer(partUnit *part_units, vector<int> &node_offsets, int _id,
                                           unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl,
                                           vector<vector<double >> &paa_mu_part_units) const;
    void fuzzyFirstLayer(partUnit *part_units, const int *nav_ids,
                         unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl,
                         vector<vector<double>> &paa_mu_part_units) const;

    int getMaxHeight();
    int getNodeNum();
    int getTotalSize();
    int getSumHeight();
    int get1stLayerInterNodesNo();
    int get1stLayerInterNodeSeriesNo();
    int get1stLayerNodesNo();


public:
    const static int power_2[];
    static int* mask;
    static unsigned short *saxes;
    static float *paas;
    static int a,b,c;
    static int ***combines;
    static int * combine_num;

    static void loadCombines();

    unsigned short sax[Const::segmentNum]{};
    vector<int> chosenSegments{};
    vector<DumpyNode*>children;
    int size = 0;
    int id = -1;
    string file_id{};   // to identify a particular leaf node
    int layer = 0;
    int leaf_num = 0;
    int bits_cardinality[Const::segmentNum]{};
    int partition_id = -1;
    vector<int> offsets{};

    DumpyNode *route(const unsigned short *_sax);
    static DumpyNode *BuildIndex(string &datafn, string &saxfn);

    [[nodiscard]] string getFileName() const{
        if(layer == 1)  return "1_" + to_string(partition_id);
        return to_string(layer) + "-" + file_id;
    }
    [[nodiscard]] bool isLeafNode() const   {return children.empty();}
    [[nodiscard]] bool isInternalNode() const {return  size > Const::th;}
    void search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const;
    void searchDTWSIMD(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
                                  float* upperLemire, float* lowerLemire) const;

    static DumpyNode*  BuildIndexFuzzy(const string & datafn, const string & saxfn, const string &paafn, vector<vector<int>>* g);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & partition_id; ar & layer;
        ar & file_id;   ar & id;
        ar & size;
        ar & sax;   ar & bits_cardinality;
        ar & children;
        ar & chosenSegments;
    }
    void save2Disk(const string &output) {
        ofstream ofs(output.c_str(), ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << (*this);
        ofs.close();
    }
    static DumpyNode *loadFromDisk(const string &saxfn, const string &idxfn, bool need_sax);

    void getIndexStats();
    int getLeafNodeNum();
    int assignLeafNum();
    int getBiasLeafNodeNum();

    DumpyNode *route1step(const unsigned short *_sax);

    void search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
                std::unordered_set<float *, createhash, isEqual> *hash_set) const;

    static long generateSaxTbl();

    static int partition(partUnit *nodes_map, int chosen_segment_number);


    void determineFanout(int *lambda_min, int *lambda_max) const;

    void determineSegments();

    double compute_score(vector<int> &node_sizes, int *plan, int lambda, vector<double> &data_seg_stdev) const;

    void visitPlanFromBaseTable(unordered_set<int> &visited, int cur_lambda, const int *plan, vector<int> &base_tbl,
                                double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                vector<double> &data_seg_stdev, double base_score);

    void searchDTW(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const;
};

struct NODE_RECORDER{
    int actual_size{0};
    vector<int>series_index_list{};

    explicit NODE_RECORDER(int act_size, DumpyNode *node) {
        actual_size = act_size;
        // internal node stores all series index list
        series_index_list.resize(node->size);
        for(int i=0;i<node->size;++i)
            series_index_list[i] = i;
    }

    NODE_RECORDER(int size, vector<int>&index_list){
        actual_size = size;
        series_index_list.resize(index_list.size());
        copy(index_list.begin(),  index_list.end(), series_index_list.begin());
    }

    NODE_RECORDER(){actual_size = 0;}
};


#endif //DUMPY_DUMPYNODE_H
