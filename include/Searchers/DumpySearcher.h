//
// Created by Zeyu Wang on 2022/1/16.
//

#ifndef DUMPY_DUMPYSEARCHER_H
#define DUMPY_DUMPYSEARCHER_H
#include <vector>
#include "../DataStructures/DumpyNode.h"

class DumpySearcher {
public:
    static void approxSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                      vector<PqItemSeries *> *heap, const string &index_dir);

    static vector<PqItemSeries *> *
    approxSearch(DumpyNode *root, float *query, int k, vector<vector<int>> *g, const string &index_dir);

    static vector<PqItemSeries *> *exactSearch(DumpyNode *root, float *query, int k, vector<vector<int>> *g);

    static vector<PqItemSeries*>* exactSearchDTW(DumpyNode* root, float *query, int k, vector<vector<int>> *g);

    static vector<PqItemSeries *> *ngSearch(DumpyNode *root, float *query, int k, int nprobes);

    static vector<PqItemSeries *> *ngSearchFuzzy(DumpyNode *root, float *query, int k, int nprobes);

    static vector<PqItemSeries *> *
    approxIncSearch(DumpyNode *root, float *query, int k, const string &index_dir, int node_num);

    static void
    approxIncSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                             vector<PqItemSeries *> *heap, const string &index_dir,
                             int &node_num);

    static void approxIncSearchInterNodeFuzzy(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                              vector<PqItemSeries *> *heap, const string &index_dir, int &node_num,
                                              unordered_set<float *, createhash, isEqual> *hash_set);

    static void approxIncSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                 vector<PqItemSeries *> *heap, const string &index_dir,
                                                 int &node_num, unordered_set<DumpyNode*>&visit,
                                                 unordered_set<float*, createhash, isEqual>*hash_set);

    static void approxIncSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                        vector<PqItemSeries *> *heap, const string &index_dir,
                                                        int &node_num, unordered_set<DumpyNode*>&visit);

    static vector<PqItemSeries *> *
    approxIncSearchFuzzy(DumpyNode *root, float *query, int k, const string &index_dir, int node_num);

    static vector<PqItemSeries *> *
    approxSearchDTW(DumpyNode *root, float *query, int k, vector<vector<int>> *g, const string &index_dir);

    static void
    approxSearchInterNodeDTW(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                             vector<PqItemSeries *> *heap,
                             const string &index_dir);

    static void approxIncSearchInterNodeDTW(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                           vector<PqItemSeries *> *heap, const string &index_dir,int &node_num);

    static vector<PqItemSeries *> * approxIncSearchDTW(DumpyNode *root, float *query, int k, const string &index_dir,
                                                                      int node_num);
};


#endif //DUMPY_DUMPYSEARCHER_H
