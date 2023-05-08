#include <iostream>
#include <cstdio>
#include <vector>
#include <thread>
#include <chrono>
#include "../include/DataStructures/DumpyNode.h"
#include "../include/DataStructures/GraphConstruction.h"
#include "../include/Searchers/DumpySearcher.h"
#include "../include/Utils/FileUtil.h"
#include "../include/Utils/MathUtil.h"
#include "../include/Utils/TimeSeriesUtil.h"
using namespace std;

vector<vector<int>>* loadGraphSkeleton(){
    int vd = 0;
    for(int i=1; i<=Const::bitsReserve;++i)
        vd += MathUtil::nChooseK(Const::segmentNum, i);
    auto nnList = new vector<vector<int>>(Const::vertexNum, vector<int>(vd, -1));

    if(!FileUtil::checkFileExists(Const::graphfn.c_str())){
        cout << "File not exists!" << Const::graphfn << endl;
        exit(-1);
    }
    FILE *f = fopen(Const::graphfn.c_str(), "rb");

    for(int i=1;i<Const::vertexNum;++i)
        fread(&((*nnList)[i][0]), sizeof(int), vd, f);

    return nnList;

}

void constructGraph(){
    GraphConstruction::buildAndSave2Disk();
}

void buildDumpy(){
    auto g = loadGraphSkeleton();
    DumpyNode* root = DumpyNode::BuildIndex(Const::datafn, Const::saxfn);
    root->save2Disk(Const::idxfn + "root.idx");
}

void approxSearchOneNode() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::idxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void approxSearchMoreNode() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearch(root, queries + i * Const::tsLength,
                                                                           Const::k, Const::idxfn,
                                                                           Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void approxSearchOneNodeDTW() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearchDTW(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::idxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void approxSearchMoreNodeDTW() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearchDTW(root, queries + i * Const::tsLength,
                                                                           Const::k, Const::idxfn,
                                                                           Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void buildDumpyFuzzy(){
    auto g = loadGraphSkeleton();
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode* root = DumpyNode::BuildIndexFuzzy(Const::datafn, Const::saxfn, Const::paafn, g);
    root->save2Disk(Const::fuzzyidxfn + "root.idx");
}

void approxSearchOneNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::fuzzyidxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void approxSearchMoreNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        auto start = chrono::system_clock::now();
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearchFuzzy(root, queries + i * Const::tsLength,
                                                                                Const::k, Const::fuzzyidxfn,
                                                                                Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void exactSearchDumpy() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *exactKnn = DumpySearcher::exactSearch(root, queries + i * Const::tsLength, Const::k, g);
        Const::logPrint("Results:");
        for (int j = 0; j < exactKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*exactKnn)[j]->ts) << endl;
    }
}

void exactSearchDumpyDTW() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *exactKnn = DumpySearcher::exactSearchDTW(root, queries + i * Const::tsLength, Const::k, g);
        Const::logPrint("Results:");
        for (int j = 0; j < exactKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*exactKnn)[j]->ts) << endl;
    }
}

void ngSearchDumpy() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->assignLeafNum();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::ngSearch(root, queries + i * Const::tsLength,
                                                                    Const::k, Const::nprobes);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void ngSearchDumpyFuzzy() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->assignLeafNum();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::ngSearchFuzzy(root, queries + i * Const::tsLength,
                                                                    Const::k, Const::nprobes);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void statIndexDumpy(){
    DumpyNode* root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->getIndexStats();
}

void statIndexDumpyFuzzy(){
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode* root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    root->getIndexStats();
}

int main() {
    Const::readConfig();

    switch (Const::index) {
        case 0:
            constructGraph();
            break;
        case 1:
            switch (Const::ops) {
                case 0:
                    buildDumpy();
                    break;
                case 1:
                    approxSearchOneNode();
                    break;
                case 2:
                    exactSearchDumpy();
                    break;
                case 3:
                    statIndexDumpy();
                    break;
                case 4:
                    approxSearchMoreNode();
                    break;
                case 5:
                    approxSearchOneNodeDTW();
                    break;
                case 6:
                    approxSearchMoreNodeDTW();
                    break;
                case 7:
                    ngSearchDumpy();
                    break;
                case 8:
                    exactSearchDumpyDTW();
                    break;
                default:
                    break;
            }
            break;
        case 2:
            if(Const::ops == 0){
                buildDumpyFuzzy();
                break;
            }
            else if(Const::ops == 1){
                approxSearchOneNodeFuzzy();
                break;
            }else if(Const::ops == 3){
                statIndexDumpyFuzzy();
                break;
            }else if(Const::ops == 4){
                approxSearchMoreNodeFuzzy();
                break;
            }else if(Const::ops == 7){
                ngSearchDumpyFuzzy();
                break;
            }
            break;
        default:    break;
    }
}
