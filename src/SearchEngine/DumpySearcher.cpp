//
// Created by pengwang5 on 2022/1/16.
//

#include <set>
#include <unordered_set>
#include <chrono>
#include <cmath>
#include <queue>
#include "../../include/Searchers/DumpySearcher.h"
#include "../../include/DataStructures/PqItemSeries.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Const.h"

static DumpyNode* targetNode;
vector<PqItemSeries *> * DumpySearcher::approxSearch(DumpyNode *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    DumpyNode *cur = (root->children)[head];
    if(cur == nullptr){
        DumpyNode *node = nullptr;
        for(int i=0;i<DumpyNode::a + DumpyNode::b + DumpyNode::c;++i){
            if(root->children[(*g)[head][i]] != nullptr){
                node = root->children[(*g)[head][i]];
                break;
            }
        }
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        }else { node->search(k, queryTs, *heap, index_dir); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->search(k, queryTs, *heap, index_dir); targetNode = cur;}
    }else approxSearchInterNode(cur, queryTs, sax, k, heap, index_dir);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void DumpySearcher::approxSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir) {
    DumpyNode *cur = root->route(sax);
    if(!cur->isInternalNode()){
        cur->search(k, queryTs, *heap, index_dir);
        targetNode = cur;
        return;
    }

    // below is only for a nullptr target leaf node, then we search the nearest sibling
    double min_dist = numeric_limits<double>::max(), max_size = 0;
    DumpyNode *node;
    for(int i=0;i<cur->children.size();++i){
        if(cur->children[i] == nullptr)  continue;
        double dist;
        if(!cur->children[i]->isInternalNode())
            dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->sax, cur->bits_cardinality, cur->children[i]->chosenSegments, i);
//        if(cur->children[i]->isLeafNode())  dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        else dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }else if(dist == min_dist && cur->children[i]->size > max_size){
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }
    }

    // we only concern whether the nearest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->search(k, queryTs, *heap, index_dir); targetNode = node;}
}

vector<PqItemSeries *> * DumpySearcher::approxSearchDTW(DumpyNode *root, float *query, int k, vector<vector<int>> *g,
                                                        const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    DumpyNode *cur = (root->children)[head];
    if(cur == nullptr){
        DumpyNode *node = nullptr;
        for(int i=0;i<DumpyNode::a + DumpyNode::b + DumpyNode::c;++i){
            if(root->children[(*g)[head][i]] != nullptr){
                node = root->children[(*g)[head][i]];
                break;
            }
        }
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNodeDTW(node, queryTs, sax, k, heap, index_dir);
        }else { node->searchDTW(k, queryTs, *heap, index_dir); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->searchDTW(k, queryTs, *heap, index_dir); targetNode = cur;}
    }else approxSearchInterNodeDTW(cur, queryTs, sax, k, heap, index_dir);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void DumpySearcher::approxSearchInterNodeDTW(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir) {
    DumpyNode *cur = root->route(sax);
    if(!cur->isInternalNode()){
        cur->searchDTW(k, queryTs, *heap, index_dir);
        targetNode = cur;
        return;
    }

    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    SaxUtil::lower_upper_lemire(queryTs->ts, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    double *lowerPaa = SaxUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    double *upperPaa = SaxUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);


    // below is only for a nullptr target leaf node, then we search the nearest sibling
    double min_dist = numeric_limits<double>::max(), max_size = 0;
    DumpyNode *node;
    for(auto & i : cur->children){
        if(i == nullptr)  continue;
        double dist = SaxUtil::minidist_paa_to_isax_DTW(upperPaa, lowerPaa, i->sax, i->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = i->size;
            node = i;
        }else if(dist == min_dist && i->size > max_size){
            max_size = i->size;
            node = i;
        }
    }

    // we only concern whether the nearest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNodeDTW(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->searchDTW(k, queryTs, *heap, index_dir); targetNode = node;}
}
static double *low_paa, *up_paa;
static float *t_paa;
bool comp_Dumpy_dtw(const DumpyNode* x, const DumpyNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::minidist_paa_to_isax_DTW(up_paa, low_paa, x->sax, x->bits_cardinality) < SaxUtil::minidist_paa_to_isax_DTW(up_paa, low_paa, y->sax, y->bits_cardinality);
}
vector<PqItemSeries *> * DumpySearcher::approxIncSearchDTW(DumpyNode *root, float *query, int k, const string &index_dir,
                                                           int node_num) {
    auto* queryTs = new TimeSeries(query);
    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    SaxUtil::lower_upper_lemire(query, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    low_paa= SaxUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    up_paa = SaxUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);

    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNodeDTW(root, queryTs, sax, k, heap, index_dir, node_num);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void DumpySearcher::approxIncSearchInterNodeDTW(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                vector<PqItemSeries *> *heap, const string &index_dir,int &node_num) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    DumpyNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(!cur->isInternalNode()){
            cur->searchDTW(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNodeDTW(cur, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<DumpyNode*>candidates;
    unordered_set<DumpyNode*>cands;
    for(DumpyNode *node: parent->children)
        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();
    sort(candidates.begin(), candidates.end(), comp_Dumpy_dtw);



    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(!candidates[i]->isInternalNode()) {
            candidates[i]->searchDTW(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            approxIncSearchInterNodeDTW(candidates[i], queryTs, sax, k, heap, index_dir, node_num);
        }
    }

}

struct PqItem{
    DumpyNode* parent;
    int id;
    double dist{};

    PqItem(DumpyNode* _parent, int _id, double d){ parent = _parent;id= _id; dist = d;}
    PqItem(){ id=-1; dist = 0;parent = nullptr;}

    bool operator <(const PqItem & pt) const{
        if(dist != pt.dist)
            return dist < pt.dist;
        if(parent->layer != pt.parent->layer)
            return parent->layer > pt.parent->layer;
        return parent < pt.parent;
    }

    bool operator >(const PqItem& pt) const{
        if(dist != pt.dist)
            return dist > pt.dist;
        if(parent->layer != pt.parent->layer)
            return parent->layer < pt.parent->layer;
        return parent > pt.parent;
    }
};

struct cmp_PqItem{
    bool operator()(const PqItem& a, const PqItem& pt) const{
        if(a.dist != pt.dist)
            return a.dist > pt.dist;
        if(a.parent->layer != pt.parent->layer)
            return a.parent->layer < pt.parent->layer;
        return a.parent > pt.parent;
    }
};

vector<PqItemSeries*>*DumpySearcher::exactSearch(DumpyNode* root, float *query, int k, vector<vector<int>> *g){

    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::idxfn);
    unordered_set<DumpyNode*>visited;
    visited.insert(targetNode);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItem>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.count(root->children[i]) > 0)    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItem(root, i, dist));
    }

    double top_dist;
    DumpyNode* node;
    int len;
    while(!pq.empty()){
        top_dist = pq.begin()->dist;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->children[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);
        if(node->isInternalNode()){
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);
//                double dist_simd = SaxUtil::LowerBound_Paa_iSax_SIMD(queryTs->paa, node->sax,
//                                                                 node->bits_cardinality, node->chosenSegments, i);

                if(dist < bsf){
                    pq.insert(PqItem(node, i, dist));
                }
            }
        }else{
//            node->search_SIMD_reordered(k, queryTs, *heap, Const::idxfn);
            node->search(k,queryTs,*heap,Const::idxfn);
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries*>*DumpySearcher::exactSearchDTW(DumpyNode* root, float *query, int k, vector<vector<int>> *g){
    vector<PqItemSeries*>* heap = approxSearchDTW(root, query, k, g, Const::idxfn);
    unordered_set<DumpyNode*>visited;
    char _ = 0;
    visited.insert(targetNode);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);
    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    SaxUtil::lower_upper_lemire(query, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    double *lowerPaa = SaxUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    double *upperPaa = SaxUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);


    priority_queue<PqItem, vector<PqItem>, cmp_PqItem> pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.count(root->children[i]) > 0)    continue;
        double dist = SaxUtil::getMinDist1stLayerDTW(upperPaa, lowerPaa, i);
        pq.emplace(root, i, dist);
    }

    int len;
    PqItem cur; double top_dist, local_bsf; DumpyNode* node;

    while(!pq.empty()){
        cur = pq.top();
        pq.pop();
        top_dist = cur.dist;
        if(top_dist >= bsf)  break;
        node = cur.parent->children[cur.id];
        if(visited.count(node) > 0) continue;
        visited.insert(node);
        if(node->isInternalNode()){
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::minidist_paa_to_isax_DTW( upperPaa, lowerPaa, node->sax,
                                             node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf)
                    pq.emplace(node, i, dist);
            }
        }else{
            node->searchDTWSIMD(k, queryTs, *heap, Const::idxfn, upperLemire, lowerLemire);
//            node->searchDTW(k, queryTs, *heap, Const::idxfn);
            bsf = (*heap)[0]->dist;
        }
    }

    delete queryTs;
    sort_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}


void searchSubTree(DumpyNode *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap, const string &index_dir,
                   int &node_num){
    unordered_set<DumpyNode*>visited;
    if(!root->isInternalNode()){
        if(root != targetNode) {
            root->search(k, queryTs, *heap, index_dir);
            --node_num;
        }
        return;
    }
    for(auto child:root->children){
        if(child == nullptr || child == targetNode || visited.count(child) > 0)    continue;
        visited.insert(child);
        if(!child->isInternalNode()){
            child->search(k,queryTs,*heap, index_dir);
            --node_num;
        }else{
            searchSubTree(child, queryTs, k, heap, index_dir, node_num);
        }
    }
}

void searchSubTree(DumpyNode *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap, const string &index_dir,
                   int &node_num, unordered_set<float*, createhash, isEqual>*hash_set){
    unordered_set<DumpyNode*>visited;
    if(!root->isInternalNode()){
        if(root != targetNode) {
            root->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        return;
    }
    for(auto child:root->children){
        if(child == nullptr || child == targetNode || visited.count(child) > 0)    continue;
        visited.insert(child);
        if(!child->isInternalNode()){
            child->search(k,queryTs,*heap, index_dir, hash_set);
            --node_num;
        }else{
            searchSubTree(child, queryTs, k, heap, index_dir, node_num, hash_set);
        }
    }
}


vector<PqItemSeries *> *DumpySearcher::ngSearch(DumpyNode *root, float *query, int k, int nprobes){
    auto *queryTs = new TimeSeries(query);
    auto heap = new vector<PqItemSeries*>();
    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    DumpyNode* root_subtree = root->route1step(sax);
    unordered_set<DumpyNode*>visited;
    if(root_subtree){
        if(!root_subtree->isInternalNode()){
            visited.insert(root_subtree);
            root_subtree->search(k, queryTs, *heap, Const::idxfn);
            --nprobes;
        }else if(root_subtree->leaf_num <= Const::pre_read){
            // a small subtree
            if(nprobes >= root_subtree->leaf_num){
                int _ =root_subtree->leaf_num;
                visited.insert(root_subtree);
                searchSubTree(root_subtree, queryTs, k, heap, Const::idxfn, _);
                nprobes -= root_subtree->leaf_num;
            }else{
                int rest = nprobes;
                approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::idxfn, rest,  visited);
                nprobes = rest;
            }
        }else{
            // a big subtree
            int to_search = min(nprobes, Const::pre_read);
            int _ = to_search;
            approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::idxfn, to_search,  visited);
            nprobes = nprobes - _ + to_search;
        }
    }

    if(nprobes <= 0){
        delete queryTs;
        sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
        return heap;
    }

    double bsf = heap->size() < k ? numeric_limits<double>::max(): (*heap)[0]->dist;

    set<PqItem>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.count(root->children[i]) > 0)    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItem(root,i, dist));
    }
    int cur_probe = 0;
    while(!pq.empty() && cur_probe < nprobes){
        double top_dist;
        DumpyNode* node;
        top_dist = pq.begin()->dist;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->children[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);

        if(node->isInternalNode()){
            int len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    pq.insert(PqItem(node, i, dist));
                }
            }
        }else{
            node->search(k, queryTs, *heap, Const::idxfn);
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

vector<PqItemSeries *> *DumpySearcher::ngSearchFuzzy(DumpyNode *root, float *query, int k, int nprobes){
    auto *queryTs = new TimeSeries(query);
    auto heap = new vector<PqItemSeries*>();
    auto*hash_set = new unordered_set<float*, createhash, isEqual>();
    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];
//    approxSearchInner(root,queryTs, k,heap,g,sax, Const::fidxfn);
//    --nprobes;

    DumpyNode* root_subtree = root->route1step(sax);
    unordered_set<DumpyNode*>visited;
    if(root_subtree){
        if(!root_subtree->isInternalNode()){
            visited.insert(root_subtree);
            root_subtree->search(k, queryTs, *heap, Const::fuzzyidxfn, hash_set);
            --nprobes;
        }else if(root_subtree->leaf_num <= Const::pre_read){
            // a small subtree
            if(nprobes >= root_subtree->leaf_num) {
                int _ = root_subtree->leaf_num;
                visited.insert(root_subtree);
                searchSubTree(root_subtree, queryTs, k, heap, Const::fuzzyidxfn, _, hash_set);
                nprobes -= root_subtree->leaf_num;
            }else{
                int rest = nprobes;
                approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fuzzyidxfn, rest,  visited,
                                         hash_set);
                nprobes = rest;
            }
        }else{
            // a big subtree
            int to_search = min(nprobes, Const::pre_read);
            int _ = to_search;
            approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fuzzyidxfn, to_search,  visited,
                                     hash_set);
            nprobes = nprobes - _ + to_search;
        }
    }

    if(nprobes <= 0){
        delete queryTs;
        sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
        return heap;
    }

    double bsf = heap->size() < k ? numeric_limits<double>::max(): (*heap)[0]->dist;

    set<PqItem>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.count(root->children[i]) > 0)    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItem(root,i, dist));
    }
    int cur_probe = 0;
    while(!pq.empty() && cur_probe < nprobes){
        double top_dist;
        DumpyNode* node;
        top_dist = pq.begin()->dist;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->children[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);

        if(node->isInternalNode()){
            int len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    pq.insert(PqItem(node, i, dist));
                }
            }
        }else{
            node->search(k, queryTs, *heap, Const::fuzzyidxfn, hash_set);
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

bool comp_Dumpy(const DumpyNode* x, const DumpyNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::LowerBound_Paa_iSax(t_paa, x->sax, x->layer) < SaxUtil::LowerBound_Paa_iSax(t_paa, y->sax, y->layer);
}

vector<PqItemSeries *> * DumpySearcher::approxIncSearch(DumpyNode *root, float *query, int k, const string &index_dir,
                                                        int node_num) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNode(root, queryTs, sax, k, heap, index_dir, node_num);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void DumpySearcher::approxIncSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,int &node_num) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    DumpyNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNode(cur, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<DumpyNode*>candidates;
    unordered_set<DumpyNode*>cands;
    for(DumpyNode *node: parent->children)
        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();
    sort(candidates.begin(), candidates.end(), comp_Dumpy);


    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(!candidates[i]->isInternalNode()) {
            candidates[i]->search(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            approxIncSearchInterNode(candidates[i], queryTs, sax, k, heap, index_dir, node_num);
        }
    }

}


vector<PqItemSeries *> * DumpySearcher::approxIncSearchFuzzy(DumpyNode *root, float *query, int k, const string &index_dir,
                                                             int node_num) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();
    auto*hash_set = new unordered_set<float*, createhash, isEqual>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNodeFuzzy(root, queryTs, sax, k, heap, index_dir, node_num, hash_set);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}


void DumpySearcher::approxIncSearchInterNodeFuzzy(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                  vector<PqItemSeries *> *heap, const string &index_dir, int &node_num,
                                                  unordered_set<float*, createhash, isEqual>*hash_set) {
    if(root->isLeafNode() || node_num <= 0)  return;
    DumpyNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(cur->isLeafNode()){
            cur->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNode(cur, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<DumpyNode*>candidates;
    unordered_set<DumpyNode*>cands;
    for(DumpyNode *node: parent->children)
        if(node != nullptr && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();

    sort(candidates.begin(), candidates.end(), comp_Dumpy);

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(candidates[i]->isLeafNode()) {
            candidates[i]->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        else {
            approxIncSearchInterNodeFuzzy(candidates[i], queryTs, sax, k, heap, index_dir, node_num, hash_set);
        }
    }

}

// for ng-search fuzzy
void DumpySearcher::approxIncSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, unordered_set<DumpyNode*>&visit,
                                             unordered_set<float*, createhash, isEqual>*hash_set) {
    if(node_num <= 0)  return;
    DumpyNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->leaf_num > node_num) {
        parent = cur;
        cur = cur->route1step(sax);
    }

    if(cur!= nullptr){
        visit.insert(cur);
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }else{
            searchSubTree(cur, queryTs,k, heap, index_dir, node_num, hash_set);
        }
    }

    if(node_num <=0)    return;

    double bsf = (*heap)[0]->dist;
    vector<PqItem>candidates;
    int len = (1 << (parent->chosenSegments.size()));
    for(int i =0;i<len;++i){
        if(parent->children[i] == nullptr || parent->children[i] == cur)    continue;
        double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
        if(dist < bsf){
            candidates.emplace_back(parent, i, dist);
        }
    }
    sort(candidates.begin(),  candidates.end());

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        DumpyNode* node = parent->children[candidates[i].id];
        if(candidates[i].dist > bsf)    break;
        if(visit.count(node) > 0) continue;
        visit.insert(node);
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        else {
            searchSubTree(node, queryTs,k, heap, index_dir, node_num, hash_set);
        }
    }

}


// for ng-search
void DumpySearcher::approxIncSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, unordered_set<DumpyNode*>&visit) {
    if(node_num <= 0)  return;
    DumpyNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->leaf_num > node_num) {
        parent = cur;
        cur = cur->route1step(sax);
    }

    if(cur!= nullptr){
        visit.insert(cur);
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            searchSubTree(cur, queryTs,k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    double bsf = (*heap)[0]->dist;
    vector<PqItem>candidates;
    int len = (1 << (parent->chosenSegments.size()));
    for(int i =0;i<len;++i){
        if(parent->children[i] == nullptr || parent->children[i] == cur)    continue;
        double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
        if(dist < bsf){
            candidates.emplace_back(parent, i, dist);
        }
    }
    sort(candidates.begin(),  candidates.end());

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        DumpyNode* node = parent->children[candidates[i].id];
        if(candidates[i].dist > bsf)    break;
        if(visit.count(node) > 0) continue;
        visit.insert(node);
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            searchSubTree(node, queryTs,k, heap, index_dir, node_num);
        }
    }

}
