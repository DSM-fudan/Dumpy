//
// Created by Zeyu Wang on 2022/1/13.
//

#include <thread>
#include <cassert>
#include <cmath>
#include "../../include/DataStructures/DumpyNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/SaxUtil.h"

unsigned short *DumpyNode::saxes = nullptr;
float *DumpyNode::paas = nullptr;
int DumpyNode::a = MathUtil::nChooseK(Const::segmentNum, 1), DumpyNode::b = MathUtil::nChooseK(Const::segmentNum, 2), DumpyNode::c = MathUtil::nChooseK(Const::segmentNum, 3);
const int DumpyNode::power_2[]{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
int * DumpyNode::combine_num = nullptr;
int*** DumpyNode::combines = nullptr;

struct pack{
    vector<bool> cur_bits;
    vector<bool> cur_mask;
    int total_size;
    int pid{};
    int masked_bits_num;
    bool disabled;

    pack(){
        total_size = 0;
        masked_bits_num = 0;
        disabled = false;
    }
    pack(partUnit* node, int chosen_segment_num, int _pid){
        total_size = node->size;
        masked_bits_num = 0;
        cur_bits.resize(chosen_segment_num, false);
        cur_mask.resize(chosen_segment_num, false);
        int _id =  node->id;
        for(int i=0;i<chosen_segment_num;++i){
            cur_bits[chosen_segment_num - 1- i] = _id % 2;
            _id >>= 1;
        }
        pid = _pid;
        node->pid = pid;
        disabled = false;
    }

    int calc_cost(int _id, int chosen_seg_num){
        int cost = 0;
        for(int i=0;i<chosen_seg_num;++i){
            if(!cur_mask[chosen_seg_num-1-i] && cur_bits[chosen_seg_num - 1 -i]!= (_id %2)){
                ++cost;
            }
            _id >>=1;
        }
        return cost;
    }

    int calc_pack_merge_cost(const pack & p, int chosen_seg_num, int *cur_cost, int *tar_cost){
        *cur_cost = 0; *tar_cost = 0;
        int cost = 0;
        for(int i=0;i<chosen_seg_num;++i){
            if(cur_mask[i] && p.cur_mask[i])    continue;
            if(cur_mask[i] && !p.cur_mask[i]){
                (*tar_cost)++;
                ++cost;
            }
            else if (!cur_mask[i] && p.cur_mask[i]) { (*cur_cost)++; ++cost;}
            else if(cur_bits[i] != p.cur_bits[i])   {
                (*cur_cost)++;
                (*tar_cost)++;
                ++cost;
            }
        }
        return cost;
    }

    void merge_pack(pack* p, int chosen_seg_num){
        pack * dis_one, * res_one;
        if(pid < p->pid){
            dis_one = p;
            res_one = this;
        }else{
            dis_one = this;
            res_one = p;
        }
        dis_one->disabled = true;
        res_one->total_size = total_size + p->total_size;
        for(int i=0;i<chosen_seg_num;++i){
            if(dis_one->cur_mask[i] && res_one->cur_mask[i])    continue;
            if(dis_one->cur_mask[i] && !res_one->cur_mask[i]) { res_one->cur_mask[i] = true; res_one->masked_bits_num++;}
            else if(!dis_one->cur_mask[i] && res_one->cur_mask[i]) continue;
            else if(dis_one->cur_bits[i] != res_one->cur_bits[i]){
                res_one->masked_bits_num++;
                res_one->cur_mask[i] = true;
            }
        }

    }

    void insert(partUnit* node, int chosen_seg_num){
        node->pid = pid;
        int _id = node->id;
        total_size += node->size;
        for(int i=0;i<chosen_seg_num;++i){
            if(!cur_mask[chosen_seg_num-1-i] && cur_bits[chosen_seg_num - 1 -i]!= (_id %2)){
                cur_mask[chosen_seg_num-1-i] = true;
                masked_bits_num++;
            }
            _id>>=1;
        }
    }

    static bool comp_size(const pack&x, const pack&y){
        return x.total_size < y.total_size;
    }
};

void DumpyNode::loadCombines(){
    string base = "../combines/" + to_string(Const::segmentNum) + "-";
    auto ret = new int**[Const::segmentNum + 1];
    combine_num = new int[Const::segmentNum];
    ifstream ff("../combines/cnum-"+ to_string(Const::segmentNum) + ".txt", ios::in);
    for(int i=0;i<Const::segmentNum;++i){
        ff >> combine_num[i];
    }
    ff.close();

    for(int i=1;i<=Const::segmentNum - 1;++i){
        ret[i] = new int*[combine_num[i]];
        ifstream f(base + to_string(i) + ".txt", ios::in);
        for(int j=0;j<combine_num[i];++j){
            ret[i][j] = new int[i];
            for(int k=0;k<i;++k) {
                f >> ret[i][j][k];
            }
        }
        f.close();
    }
    ret[Const::segmentNum] = new int*[1];
    ret[Const::segmentNum][0] = new int[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i){
        ret[Const::segmentNum][0][i] = i;
    }
    combines = ret;

}

void materializeAllLeavesWithSax(string datafn, DumpyNode* root, int *navids, string index_dir, unsigned short*sax_tbl){
    auto start_sax = chrono::system_clock::now();
    Const::logPrint("Start move sax to disk file in 1st layer.");

    unordered_map<DumpyNode*, vector<unsigned short *>>sax_buffer;
    for(int i=0;i<root->size;++i){
        auto * sax = sax_tbl + i * Const::segmentNum;
        DumpyNode* node = root->route(sax);
        sax_buffer[node].push_back(sax);
    }
    for(auto &[node, buffer]:sax_buffer){
        string outfile = Const::idxfn + node->getFileName() + "_sax";
        if(node->partition_id == -1)    outfile += "_L";
        FILE *outf = fopen(outfile.c_str(), "a");
        for(auto sax:buffer)
            fwrite(sax, sizeof(unsigned short ), Const::segmentNum, outf);
        fclose(outf);
    }
    sax_buffer.clear();

    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<DumpyNode*, LBL_UNIT>lbl;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        lbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();

        //copy series to node to ensure write is only called once upon a node fbl
        for(long i = cur; i<cur+num;++i){
            DumpyNode* node = root->route(sax_tbl + i * Const::segmentNum);
            lbl[node].buffer.push_back(tss + (i-cur) * Const::tsLength);
        }

        // write series in order to node file from node fbl
        for(auto & [node,lbl_unit]:lbl){
            string outfile = Const::idxfn + node->getFileName();
            if(node->partition_id == -1)  outfile += "_L";
            FILE *outf = fopen(outfile.c_str(), "a");

            for(float *dat:lbl_unit.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now Materialize all leaves. Progress: " + to_string((double)cur / (double)total * 100) + "%");

    }

    fclose(f);
    delete[] navids;
}


DumpyNode* DumpyNode::route(const unsigned short *_sax){
    if(isLeafNode())
        return this;
    int nav_id = SaxUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    if(children[nav_id] == nullptr) return this;
    return children[nav_id]->route(_sax);
}

DumpyNode* DumpyNode::route1step(const unsigned short *_sax){
    assert(!isLeafNode());
    int nav_id;
    if(layer >= 1)
        nav_id = SaxUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    else
        nav_id = SaxUtil::invSaxHeadFromSax(_sax, Const::bitsCardinality, Const::segmentNum);
    return children[nav_id];
}

void DumpyNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

//    long fs = FileUtil::getFileSize(fn.c_str());
//    int series_num = fs / Const::tsLengthBytes;
//    assert(series_num == size);

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
//    for(int i=0;i<size;++i)
//        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);
    fread(ts, sizeof(float), size * Const::tsLength, f);

    for(int i=0;i<size;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void DumpyNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir, unordered_set<float*, createhash, isEqual>*hash_set) const{
    assert(isLeafNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    auto *ts = new float[size * Const::tsLength];
    fread(ts, sizeof(float), size * Const::tsLength, f);

    for(int i=0;i<size;++i){
        if(hash_set->find(ts + i * Const::tsLength) != hash_set->end()) continue;
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->erase(heap.back()->ts);
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void DumpyNode::searchDTW(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
    for(int i=0;i<size;++i)
        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);

    for(int i=0;i<size;++i){
        double dist = TimeSeriesUtil::dtw(queryTs->ts, ts + i * Const::tsLength, Const::tsLength,Const::dtw_window_size, bsf);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void DumpyNode::searchDTWSIMD(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
                              float* upperLemire, float* lowerLemire) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
    for(int i=0;i<size;++i)
        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);

    float cb[Const::tsLength];
    float cb1[Const::tsLength];

    int length = 2*Const::dtw_window_size+1;
    float tSum[length];
    // pre_cost
    float pCost[length];
    // raw distance
    float rDist[length];
    double dist;


    for(int i=0;i<size;++i){
        float dist2=TimeSeriesUtil::lb_keogh_data_bound(ts + i * Const::tsLength, upperLemire,lowerLemire,
                                                        cb1, Const::tsLength, bsf);
        if(dist2 < bsf) {
            cb[Const::tsLength - 1] = cb1[Const::tsLength - 1];
            for (int ii = Const::tsLength - 2; ii >= 0; ii--)
                cb[ii] = cb[ii + 1] + cb1[ii];

            dist = TimeSeriesUtil::dtwsimd(queryTs->ts, ts + i * Const::tsLength, cb, Const::tsLength,
                                           Const::dtw_window_size, bsf,
                                           tSum, pCost, rDist);

            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }

            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}


DumpyNode *DumpyNode::BuildIndex(string &datafn, string &saxfn) {
    Const::logPrint("Start building index.");
    FileUtil::checkDirClean(Const::idxfn.c_str());
    loadCombines();
    long series_num = generateSaxTbl();

    auto* root = new DumpyNode();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    partUnit nodeIn1stLayer[Const::vertexNum];
    int *navids = new int[series_num];
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].size=0, nodeIn1stLayer[i].pid = -1;

    // get 1st layer node size
    for(long i=0;i<series_num;++i){
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        navids[i] = nav_id;
        nodeIn1stLayer[nav_id].size++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition 1st layer
    int partNum = partition(nodeIn1stLayer, Const::segmentNum);
    Const::logPrint("Finish partition");
    DumpyNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new DumpyNode(1, i);
    root->children.resize(Const::vertexNum);
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size <= 0) continue;
        if(nodeIn1stLayer[i].size > Const::th) {
//            assert(nodeIn1stLayer[i].size > Const::th);
            root->children[i] = new DumpyNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }else if(nodeIn1stLayer[i].pid == -1){
            root->children[i] = new DumpyNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }
        else{
            int pid = nodeIn1stLayer[i].pid;
            root->children[i] = childrenList[pid];
            childrenList[pid]->size += nodeIn1stLayer[i].size;
            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");
    
    // put data offsets to internal nodes in 1st layer
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->offsets.reserve(nodeIn1stLayer[i].size);
    for(int i=0;i<series_num;++i){
        int nav_id = navids[i];
        root->children[nav_id]->offsets.push_back(i);
    }
    Const::logPrint("data offsets have been put into nodes in the 1st layer.");

    int j = 0;
    int milestone = 0.1 * Const::vertexNum;
    Const::logPrint("start grow the index structure");
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size > Const::th) {
            root->children[i]->growIndex();
        }
        if(++j%milestone == 0)
            Const::logPrint(to_string(j) + " nodes in the 1st layer has been processed.");
    }

    Const::logPrint("build index skeleton finished.");

    Const::logPrint("Start materialize leaves");
    materializeAllLeavesWithSax(datafn, root, navids, Const::idxfn, saxes);
    Const::logPrint("build index successfully!");
    delete[] saxes;

    return root;
}

void DumpyNode::growIndex() {
    if(size <= Const::th)   return;    
    determineSegments();
    int chosen_num = chosenSegments.size();
    
    // statistic children information in order to partition
    partUnit nodes[1<<chosen_num];
    for(int i=0;i<(1<<chosen_num);++i)
        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

    for(int i=0;i<size;++i){
        int new_id = SaxUtil::extendSax(DumpyNode::saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality, chosenSegments);
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
    }

    if(this->layer > 1) vector<int>().swap(offsets);

    int partNum = partition(nodes, chosen_num);     
    
    DumpyNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new DumpyNode(this, i);
    children.resize(1 << chosen_num);
    for(int i=0;i<(1 << chosen_num);++i){
        if(nodes[i].size <= 0)  continue;
        else if(nodes[i].size > Const::th) {
            children[i] = new DumpyNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            children[i]->offsets.resize(nodes[i].size);
            copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
            vector<int>().swap(node_offsets[i]);
        }else if(partition_id == -1){
            children[i] = new DumpyNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
        else{
            int _pid = nodes[i].pid;
            children[i] = childrenList[_pid];
            childrenList[_pid]->size += nodes[i].size;
            generateSaxAndCardinality4LeafNode(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
    }
    
    vector<vector<int>>().swap(node_offsets);
    
    for(auto &child: children){
        if(child!= nullptr && child->size > Const::th){
            child->growIndex();
        }
    }

}

void DumpyNode::determineFanout(int *lambda_min, int * lambda_max) const{
    if(size < 2 * Const::th)    {
        *lambda_min = 1;
        *lambda_max = 1;
        return;
    }
    *lambda_min = -1;
    *lambda_max = -1;
    double _min = size / (Const::th * Const::f_high);
    double _max = size / (Const::th * Const::f_low);
    for(int i = 1; i <= Const::segmentNum; ++i){
        if(*lambda_min == -1){
            if((1<< i) >= _min){
                *lambda_min = i;
            }
        }else{
            if((1<<i) == _max){
                *lambda_max = i;
                break;
            }else if((1<<i) > _max){
                *lambda_max = max(i-1,*lambda_min);
                break;
            }
        }
    }
}

void DumpyNode::determineSegments() {
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);

    vector<int>unit_size(Const::vertexNum, 0);

    vector<unordered_map<unsigned short, int>>data_seg_symbols(Const::segmentNum);
    for(int offset:offsets){
        unsigned short* cur_sax = saxes + offset * Const::segmentNum;
        for(int i=0;i<Const::segmentNum;++i){
            data_seg_symbols[i][cur_sax[i]]++;
        }
        int head = SaxUtil::extendSax(cur_sax, bits_cardinality);
        unit_size[head]++;
    }

    // compute stdev of each segment
    vector<double>data_seg_mean(Const::segmentNum, 0);
    vector<double>data_seg_stdev(Const::segmentNum, 0);
    for(int i=0;i<Const::segmentNum;++i){
        auto& map = data_seg_symbols[i];
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            data_seg_mean[i] += (SaxUtil::getMidLineFromSaxSymbolbc8(symbol) * iter.second);
        }
        data_seg_mean[i] /= size;
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            double mid_value = SaxUtil::getMidLineFromSaxSymbolbc8(symbol);
            data_seg_stdev[i] += (iter.second * ((mid_value - data_seg_mean[i]) * (mid_value - data_seg_mean[i])));
        }
//        data_seg_stdev[i] = sqrt(data_seg_stdev[i] / size);
        data_seg_stdev[i] /= size;
    }

    vector<double>().swap(data_seg_mean);
    vector<unordered_map<unsigned short, int>>().swap(data_seg_symbols);

    // start to compute the size of each node in each plan
    int plan_num;
    if(lambda_max < Const::segmentNum)
        plan_num = combine_num[lambda_max];
    else
        plan_num = 1;
    unordered_set<int>visited;
    double max_score = 0;
    vector<int> best_plan;
    for(int i=0;i<plan_num;++i){
        int *plan = combines[lambda_max][i];
        // first evaluate the whole plan
        vector<int>plan_node_sizes(1<<lambda_max, 0);
        int mask_code = MathUtil::generateMaskSettingKbits(plan, lambda_max, Const::segmentNum);
        map<int,int>max_node_size;
        for(int j=0;j<Const::vertexNum;++j){
            max_node_size[mask_code & j] += unit_size[j];
        }
//        assert(max_node_size.size() == (1 << lambda_max));
        int _ = 0;
        for(auto & iter: max_node_size){
            plan_node_sizes[_++] = iter.second;
        }
        map<int,int>().swap(max_node_size);

        double score = compute_score(plan_node_sizes, plan, lambda_max, data_seg_stdev);
        if(score > max_score){
            max_score = score;
            best_plan.clear();
            for(int j = 0; j<lambda_max;++j)
                best_plan.push_back(plan[j]);
        }

        if(lambda_min <= lambda_max - 1)
            visitPlanFromBaseTable(visited, lambda_max - 1, plan, plan_node_sizes,
                                   &max_score, best_plan, lambda_min, mask_code, data_seg_stdev, score);
        vector<int>().swap(plan_node_sizes);
    }

    unordered_set<int>().swap(visited);
    vector<int>().swap(unit_size);
    chosenSegments = best_plan;
}

double DumpyNode::compute_score(vector<int>&node_sizes, int *plan, int lambda, vector<double>&data_seg_stdev) const{
    if(size < 2*Const::th){
        if(node_sizes[0] > Const::th || node_sizes[1] > Const::th)
            return (double)min(node_sizes[0],node_sizes[1]) / Const::th;
        return data_seg_stdev[plan[0]] * 100;
    }
    int over_th_nodes_no = 0;
    for(int _:node_sizes){
        if(_ > Const::th)
            over_th_nodes_no++;
    }
    double w = ((double)over_th_nodes_no) / ((int)node_sizes.size());
    double sum_seg = 0;
    for(int i=0;i<lambda;++i){
        sum_seg += data_seg_stdev[plan[i]];
    }
    sum_seg  = sqrt(sum_seg / lambda);
    sum_seg = exp(1+sum_seg);
//    sum_seg  = sum_seg / lambda;
    auto *tmp = new double[node_sizes.size()];
    for(int i=0;i<node_sizes.size();++i){
        tmp[i] = ((double)node_sizes[i]) / Const::th;
    }
    double stdev_fill_factor = MathUtil::deviation(tmp, node_sizes.size());

//    double balance = ((1+w) * stdev_fill_factor);
//    double balance = 1.0 / stdev_fill_factor;
    double balance = exp(-(1+w) * stdev_fill_factor);
    double ret = sum_seg + Const::alpha* balance;
    delete[] tmp;
    return ret;
}

void
DumpyNode::visitPlanFromBaseTable(unordered_set<int> &visited, int cur_lambda, const int *plan, vector<int> &base_tbl,
                                  double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                  vector<double> &data_seg_stdev, double base_score) {
    // base mask is used to detect base tbl
    int base_mask = 1;
    for(int i=0;i<cur_lambda;++i)
        base_mask = (base_mask << 1) +1;

    for(int i=0;i<cur_lambda + 1 ;++i){
        int reset_pos = plan[i];
        // get the whole plan code
        int cur_whole_mask = mask_code - (1 << (Const::segmentNum - 1  - reset_pos));
        if(visited.contains(cur_whole_mask))
            continue;
        visited.insert(cur_whole_mask);
        // get the new plan
        int *new_plan = new int[cur_lambda];
        int _=0;
        for(int j=0;j<cur_lambda+1;++j){
            if(i != j)  new_plan[_++] = plan[j];
        }

        // get the current base mask
        int cur_base_mask = base_mask - (1 << (cur_lambda - i));
        map<int,int>node_size_map;
        for(int j=0;j<base_tbl.size();++j)
            node_size_map[cur_base_mask &j] += base_tbl[j];
//        assert(node_size_map.size() == (1<<cur_lambda));
        vector<int>new_tbl(1<<cur_lambda, 0);
        _ = 0;
        for(auto &iter:node_size_map)
            new_tbl[_++] = iter.second;
        map<int,int>().swap(node_size_map);
        double score  = compute_score(new_tbl, new_plan, cur_lambda, data_seg_stdev);
        if(score > *max_score){
            *max_score = score;
            best_plan.clear();
            for(_ = 0; _<cur_lambda;++_)
                best_plan.push_back(new_plan[_]);
        }
//        if(score > base_score && cur_lambda > 1)
        if(cur_lambda > lambda_min)
            visitPlanFromBaseTable(visited, cur_lambda - 1, new_plan, new_tbl,
                                   max_score, best_plan, lambda_min, cur_whole_mask, data_seg_stdev, score);
        vector<int>().swap(new_tbl);
        delete[] new_plan;
    }
}

// TODO: this function may be optimized with SIMD
PAA_INFO* DumpyNode::statPaa(){
    auto* r = new PAA_INFO();
    double split_line[Const::segmentNum],paa_max[Const::segmentNum],paa_min[Const::segmentNum],paa_mu[Const::segmentNum];
    // TODO: optimize
    double lb;  // ub is the new split line
    for(int i=0; i < Const::segmentNum; ++i)
        SaxUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &split_line[i]);
    for(auto &i:paa_max) i = - numeric_limits<double>::max();
    for(auto &i:paa_min) i = numeric_limits<double>::max();
    for(auto &i:r->paa_up_size) i = 0;
    for(auto &i:r->paa_below_size) i = 0;
    for(auto &i:r->paa_variance) i = 0;
    for(auto &i:paa_mu) i=0;
    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            paa_mu[i] += value;
            paa_min[i] = min(paa_min[i], value);
            paa_max[i] = max(paa_max[i], value);
            if(value > split_line[i]) {
                r->paa_up_size[i]++;
            }
            else {
                r->paa_below_size[i]++;
            }
        }
    }
    for(double & i : paa_mu) {
        i /= size;
    }

    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            r->paa_variance[i] += (value - paa_mu[i]) * (value - paa_mu[i]);
        }
    }
    return r;
}

struct tmp{
    int i{};
    double score{};
    tmp(int _i, double _score){i=_i;score = _score;}
    tmp(){;}

    static bool order(tmp a,tmp b){
        return a.score < b.score;
    }

    static bool orderdesc(tmp a,tmp b){
        return a.score > b.score;
    }
};

int DumpyNode::chooseOneSegment(PAA_INFO* node){
    int min = numeric_limits<int>::max(), min_index = -1;
    for(int i = 0; i < Const::segmentNum;++i){
        int big = max(node->paa_up_size[i], node->paa_below_size[i]);
        if(big < min){
            min = big;
            min_index = i;
        }
    }
    return min_index;
}

void DumpyNode::chooseSegment(PAA_INFO *paa, int chosen_num) {
    chosenSegments.resize(chosen_num);
    if(chosen_num == 1) { chosenSegments[0]=chooseOneSegment(paa) ; return;}

    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        if(bits_cardinality[i] >= Const::bitsCardinality)
            scores[i] = tmp(i, -1);
        else
            scores[i] = tmp(i, paa->paa_variance[i]);
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);
    for(int i=0;i<chosen_num;++i)
        chosenSegments[i] = scores[i].i;
    sort(chosenSegments.begin(), chosenSegments.end());
}

void DumpyNode::generateSaxAndCardIn1stLayer(int new_id){
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        unsigned short t = new_id % 2 ;
        new_id >>= 1;
        sax[i] =  t;
    }
}

void DumpyNode::generateSaxAndCardinality(DumpyNode* node, int new_id){
    copy(sax, sax + Const::segmentNum, node->sax);
    copy(bits_cardinality, bits_cardinality + Const::segmentNum, node->bits_cardinality);
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = (chosenSegments)[i];
        node->bits_cardinality[seg]++;
        int t = new_id % 2 ;
        new_id >>= 1;
        node->sax[seg] = (node->sax[seg] << 1) + t;
    }
}

void DumpyNode::generateSaxAndCardIn1stLayer4LeafNode(int new_id){
    if(bits_cardinality[0] == -1){
        for(int &i:bits_cardinality)    i=1;
        generateSaxAndCardIn1stLayer(new_id);
        return;
    }
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        int t = new_id % 2 ;
        new_id >>= 1;
        if(bits_cardinality[i] == 1 && sax[i] != t){
            bits_cardinality[i] = 0;
        }
    }
}

void DumpyNode::generateSaxAndCardinality4LeafNode(DumpyNode* node, int new_id){
    if(node->bits_cardinality[0] == -1){
        generateSaxAndCardinality(node, new_id);
        return;
    }
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = chosenSegments[i];
        int t = new_id % 2 ;
        new_id >>= 1;
        if(node->bits_cardinality[seg] == bits_cardinality[seg] + 1 && node->sax[seg] % 2 != t){
            node->bits_cardinality[seg]--;
            node->sax[seg] >>= 1;
        }
    }
}

struct failUnit{
    partUnit *node{};
    int neighbor_size{};

    failUnit(partUnit* a, int b){ node = a; neighbor_size = b;}
};

static bool comp_fail_node(failUnit *x, failUnit *y){
    return x->neighbor_size > y->neighbor_size;
}

int DumpyNode::partition(partUnit* nodes_map, int chosen_segment_number){
    vector<partUnit*>nodes;
    int input_node_number = 1 << chosen_segment_number;
    int total_size = 0, node_number;
    for(int i=0;i<input_node_number;++i){
        if(nodes_map[i].size < Const::th * Const::small_perc && nodes_map[i].size > 0) {
            nodes.push_back(&nodes_map[i]);
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number <= 3) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
//    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    int first_round_pack_num = (int)floor(total_size / (Const::th));
    if(node_number <= first_round_pack_num) return 0;
    int max_mask_num = floor(chosen_segment_number * Const::max_mask_bit_percentage);
    // 0.5 is a small number, generate more packs in the first round
    vector<pack>packs(first_round_pack_num);
    sort(nodes.begin(), nodes.end(), partUnit::comp_size);
    for(int i=0;i<first_round_pack_num;++i){
        packs[i] = pack(nodes[i], chosen_segment_number, i);
    }

    for(partUnit* cur_node:nodes) {
        if (cur_node->pid != -1) continue;   // the node has been packed
        int cur_id = cur_node->id;

        int min_cost = chosen_segment_number;
        int pack_id = -1;
        for (pack &p: packs) {
            if (p.total_size + cur_node->size > Const::th || p.masked_bits_num >= max_mask_num) continue;
            int cost = p.calc_cost(cur_id, chosen_segment_number);
            if(cost + p.masked_bits_num >= max_mask_num)    continue;
            if (cost < min_cost) {
                min_cost = cost;
                pack_id = p.pid;
            }
        }
        if (pack_id == -1) {
            packs.emplace_back(cur_node, chosen_segment_number, packs.size());
        }else {
            packs[pack_id].insert(cur_node, chosen_segment_number);
        }
    }

    // merge packs
    unordered_map<int,int>pid_map;
    for(int i = 0;i<packs.size();++i)
        pid_map[i] = i;
    sort(packs.begin(),packs.end(), pack::comp_size);
    for(int i=0;i<packs.size();++i){
        pack& cur_pack = packs[i];
        if(cur_pack.pid != pid_map[cur_pack.pid])    continue;

        int min_cost = chosen_segment_number;
        int min_size = numeric_limits<int>::max();
        int min_pack_id = -1;
        int cur_cost, tar_cost;
        for(int j=0;j<packs.size();++j){
            pack&target_pack = packs[j];
            if(target_pack.disabled || cur_pack.pid == target_pack.pid || cur_pack.total_size + target_pack.total_size > Const::th
               || cur_pack.masked_bits_num >= max_mask_num || target_pack.masked_bits_num >= max_mask_num) continue;
            cur_pack.calc_pack_merge_cost(target_pack, chosen_segment_number, &cur_cost, &tar_cost);
            if(cur_cost + cur_pack.masked_bits_num >= max_mask_num ||
               tar_cost + target_pack.masked_bits_num >= max_mask_num) continue;
            int cost = cur_cost + tar_cost;
            if(cost < min_cost || (cost == min_cost && cur_pack.total_size < min_size)){
                min_cost = cost;
                min_size = target_pack.total_size;
                min_pack_id = j;
            }
        }

        if(min_size < numeric_limits<int>::max()){
            cur_pack.merge_pack(&packs[min_pack_id], chosen_segment_number);
        }
    }

    // re-assign the pids to the nodes
    int max_pid = 0;
    for(partUnit* node:nodes) {
        node->pid = pid_map[node->pid];
        max_pid = max(max_pid, node->pid);
    }

    return max_pid + 1;
}

void DumpyNode::getIndexStats(){
    int total_leaf_node_num = getLeafNodeNum();
    int total_size = getTotalSize();
    cout << "Total size = " << total_size << endl;
    cout <<"Total nodes number = " << getNodeNum() << endl;
    cout << "Leaf node number = " << total_leaf_node_num << endl;
    cout << "1st layer node number = " << get1stLayerNodesNo() <<endl;
    cout << "1st layer internal node number = " << get1stLayerInterNodesNo() << endl;
    cout << "1st layer internal series number = " << get1stLayerInterNodeSeriesNo() << endl;
    cout << "Max. height = " << getMaxHeight() - 1 <<endl;
    cout << "Avg. Height = " << getSumHeight() / (double) total_leaf_node_num << endl;
    cout <<"Avg. Filling Factor = "<< total_size / (double)total_leaf_node_num / Const::th << endl;
    cout << "Bias leaf node ratio = " << (double)getBiasLeafNodeNum() / total_leaf_node_num << endl;
}

int DumpyNode::get1stLayerInterNodesNo(){
    unordered_set<DumpyNode*>node;
    for(DumpyNode* child:children){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int DumpyNode::get1stLayerNodesNo(){
    unordered_set<DumpyNode*>node;
    for(DumpyNode* child:children){
        if(child== nullptr || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int DumpyNode::get1stLayerInterNodeSeriesNo(){
    unordered_set<DumpyNode*>node;
    int ret = 0;
    for(DumpyNode* child:children){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
        ret += child->size;
    }
    return ret;
}

int DumpyNode::getMaxHeight(){
    if(isLeafNode())    return 1;
    int max_height = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        max_height = max(child->getMaxHeight(), max_height);
        hash_map.insert(child);
    }
    return max_height + 1;
}

int DumpyNode::getLeafNodeNum(){
    if(isLeafNode())    return 1;
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

// isax bias leaf nodes number
int DumpyNode::getBiasLeafNodeNum(){
    if(isLeafNode()){
        int max_b = 0, min_b=8;
        for(int &bc:bits_cardinality){
            max_b = max(max_b, bc);
            min_b = min(min_b, bc);
        }
        return (max_b - min_b >= 4);
    }
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getBiasLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

int DumpyNode::getTotalSize(){
    if(isLeafNode())    return size;
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getTotalSize();
        hash_map.insert(child);
    }
    return sum;
}

int DumpyNode::getNodeNum(){
    if(isLeafNode())    return 1;
    int sum = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getNodeNum();
        hash_map.insert(child);
    }
    return sum + 1;
}

int DumpyNode::getSumHeight(){
    if(isLeafNode())    return layer;
    int sum_height = 0;
    unordered_set<DumpyNode*>hash_map;
    for(DumpyNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum_height += child->getSumHeight();
        hash_map.insert(child);
    }
    return sum_height;
}

int DumpyNode::loadSax(const string & saxfn){
    long f_size = FileUtil::getFileSize(saxfn.c_str()), series_num = f_size / (sizeof(unsigned short) * Const::segmentNum);
    saxes = new unsigned short [f_size / sizeof(unsigned short )];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short ), f_size / sizeof(unsigned short), f);
    fclose(f);
    Const::logPrint("Finish loading sax");
    return series_num;
}

void DumpyNode::loadPaa(const string & paafn){
    long f_size = FileUtil::getFileSize(paafn.c_str());
    paas = new float [f_size / sizeof(float )];
    FILE *f = fopen(paafn.c_str(), "rb");
    fread(paas, sizeof(float ), f_size / sizeof(float ), f);
    fclose(f);
    Const::logPrint( "Finish loading paa");
}

long DumpyNode::generateSaxAndPaaTbl(){
    string fn = Const::datafn;
    long fs = FileUtil::getFileSize(fn.c_str());
    long series_num = fs / Const::tsLengthBytes;
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    paas = new float[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");


    while(rest > 0){
        unsigned num;
        if(rest > 4000000)    num = 4000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        fread(tss, sizeof(float),num * Const::tsLength,  f);

        for(int i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                SaxUtil::paaAndSaxFromTs(tss + i * Const::tsLength,
                                         paas + (cur+ i) * Const::segmentNum, saxes + (cur+ i) * Const::segmentNum,
                                         Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
    }

    fclose(f);
    return series_num;
}

long DumpyNode::generateSaxTbl(){
    string fn = Const::datafn;
    long series_num;
    if(Const::series_num == -1) {
        long fs = FileUtil::getFileSize(fn.c_str());
        series_num = fs / Const::tsLengthBytes;
    }else{
        series_num = Const::series_num;
    }
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");


    while(rest > 0){
        unsigned num;
        if(rest > 2000000)    num = 2000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        fread(tss, sizeof(float),num * Const::tsLength,  f);

        for(int i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                SaxUtil::saxFromTs(tss + i * Const::tsLength, saxes + (cur+ i) * Const::segmentNum,
                                   Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
    }

    fclose(f);
    return series_num;
}

DumpyNode *DumpyNode::loadFromDisk(const string &saxfn, const string &idxfn, bool need_sax) {
    if(need_sax)
        loadSax(saxfn);
    ifstream ifs(idxfn, ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new DumpyNode();
    ia >> (*g);
    ifs.close();
    return g;
}

int DumpyNode::assignLeafNum() {
    if(!isInternalNode()) {
        leaf_num = 1;
        return 1;
    }

    unordered_set<DumpyNode*>visited;
    for(DumpyNode* child: children){
        if(child == nullptr || visited.count(child) > 0)    continue;
        visited.insert(child);
        leaf_num += child->assignLeafNum();
    }

    return leaf_num;
}