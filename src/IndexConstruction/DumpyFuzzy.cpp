//
// Created by Zeyu Wang on 2022/1/30.
//

#include "../../include/DataStructures/DumpyNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/MathUtil.h"
#include <thread>

// We recommend building Dumpy-fuzzy with PAA table as follows.
// Only with SAX table is also viable.

static int fuzzy_num = 0;
int *DumpyNode::mask = nullptr;

// put actual series into disk file of nodes in 1st layer
void materialize1stLayerFuzzy(string datafn, DumpyNode* root, int *navids, string index_dir, unordered_map<DumpyNode*, NODE_RECORDER>* navigating_tbl){
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<DumpyNode*, LBL_UNIT>fbl;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];
        fread(tss, sizeof(float),num * Const::tsLength,  f);

        long bound = cur + num;
        // statistic the size of each node fbl size, and allocate memory for them
        for(long i=cur;i<bound;++i){
//            fbl[root->children[navids[i]]].offsets.push_back(i);
            fbl[root->children[navids[i]]].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        for(auto &iter : *navigating_tbl){
            if(iter.first->partition_id != -1){
                NODE_RECORDER& recorder = iter.second;
                int & cur_pos = recorder.actual_size;
                while(cur_pos < recorder.series_index_list.size() && recorder.series_index_list[cur_pos] < bound){
                    int index = recorder.series_index_list[cur_pos];
//                    fbl[iter.first].offsets.push_back(index);
                    fbl[iter.first].buffer.push_back(tss + (index - cur) * Const::tsLength);
                    ++cur_pos;
                }
            }
        }

        // write series in order to node file from node fbl
        for(auto & iter:fbl){
            string outfile = index_dir ;
            if(iter.first->partition_id == -1)
                outfile += "U_" + to_string(iter.first->id);
            else outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->partition_id);
            FILE *outf = fopen(outfile.c_str(), "a");

//            for(int offset: iter.second.offsets){
//                fwrite(tss + (offset - cur) * Const::tsLength, sizeof(float), Const::tsLength, outf);
//                fwrite(&offset, sizeof(int), 1, outf);
//            }
            for(int i=0;i<iter.second.buffer.size();++i)
                fwrite(iter.second.buffer[i], sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.(FIRST STAGE)");

    }

    fseek(f, 0, SEEK_SET);
    rest = root->size; cur = 0;
    while(rest > 0){
        fbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        long bound = cur + num;

        for(auto &iter : *navigating_tbl){
            // internal node
            if(iter.first->partition_id == -1){
                NODE_RECORDER& recorder = iter.second;
                int & cur_pos = recorder.actual_size;
                while(cur_pos < recorder.series_index_list.size() && recorder.series_index_list[cur_pos] < bound){
                    int index = recorder.series_index_list[cur_pos];
//                    fbl[iter.first].offsets.push_back(index);
                    fbl[iter.first].buffer.push_back(tss + (index - cur) * Const::tsLength);
                    ++cur_pos;
                }
            }
        }

        for(auto & iter:fbl){
            string outfile = index_dir ;
            assert(iter.first->partition_id == -1);
            outfile += "U_" + to_string(iter.first->id);
            FILE *outf = fopen(outfile.c_str(), "a");

//            for(int offset: iter.second.offsets){
//                fwrite(tss + (offset - cur) * Const::tsLength, sizeof(float), Const::tsLength, outf);
//                fwrite(&offset, sizeof(int), 1, outf);
//            }

            for(int i=0;i<iter.second.buffer.size();++i)
                fwrite(iter.second.buffer[i], sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }

        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.(SECOND STAGE)");
    }

    fclose(f);
    delete[] navids;
    navigating_tbl->clear();
}

void materializeInterNodeFuzzy(DumpyNode* node, unsigned short *saxes, int actual_size, unordered_map<DumpyNode*, NODE_RECORDER>& navigating_tbl){

    FILE *f = fopen((Const::fuzzyidxfn + "U_" + to_string(node->id)).c_str(), "r");
    long rest = node->size, cur = 0, num, bound;
    unordered_map<DumpyNode*, LBL_UNIT>lbl;

    while(rest > 0){
        lbl.clear();
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        bound = cur + num;

        // actual series routing and inserting to lbl
        for(long i = cur; i < bound && i < actual_size; ++i){
            DumpyNode* target = node->route(saxes + node->offsets[i] * Const::segmentNum);
            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        // TODO: this part can be re-written as multi-thread
        // fuzzy series fetching
        for(auto &iter : navigating_tbl){
            NODE_RECORDER& recorder = iter.second;
            int & cur_pos = recorder.actual_size;
            while(cur_pos < recorder.series_index_list.size() && recorder.series_index_list[cur_pos] < bound){
                int index = recorder.series_index_list[cur_pos];
                lbl[iter.first].buffer.push_back(tss + (index - cur) * Const::tsLength);
                ++cur_pos;
            }
        }

        // write one series one time
        for(auto &iter:lbl){
            string outfile = Const::fuzzyidxfn + iter.first->getFileName();
            FILE *outf = fopen(outfile.c_str(), "a");

//            LBL_UNIT &tmp = iter.second;
//            for(int i=0;i<tmp.offsets.size();++i){
//                fwrite(tmp.buffer[i], sizeof(float), Const::tsLength, outf);
//                fwrite(&tmp.offsets[i], sizeof(int), 1, outf);
//            }
            for(int i=0;i<iter.second.buffer.size();++i)
                fwrite(iter.second.buffer[i], sizeof(float), Const::tsLength, outf);
//            for(float *dat:iter.second.buffer)
//                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }

        delete[]tss;
        rest-=num;
        cur += num;
    }
    fclose(f);
    FileUtil::FileRemove((Const::fuzzyidxfn + "U_" + to_string(node->id)).c_str());

}

DumpyNode*  DumpyNode::BuildIndexFuzzy(const string & datafn, const string & saxfn, const string &paafn, vector<vector<int>>* g){
    FileUtil::checkDirClean(Const::fuzzyidxfn.c_str());
    Const::logPrint("start building index.");

//    long series_num = loadSax(saxfn);
//    loadPaa(paafn);
    long series_num = generateSaxAndPaaTbl();

    mask = MathUtil::generateMask(Const::segmentNum);
    auto* root = new DumpyNode();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    partUnit nodeIn1stLayer[Const::vertexNum];

    int *navids = new int[series_num];
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].size=0, nodeIn1stLayer[i].pid = -1;

    auto * paa_mu_part_units = new vector<vector<double> >(Const::vertexNum, vector<double>(Const::segmentNum, 0));
    // get 1st layer node size
    for(long i=0;i<series_num;++i){
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        float *paa = paas + i * Const::segmentNum;
        for(int j = 0;j<Const::segmentNum;++j)
            (*paa_mu_part_units)[nav_id][j] += paa[j];
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
    int internal_size = 0;
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size <= 0) continue;
        for(int j = 0;j<Const::segmentNum;++j)
            (*paa_mu_part_units)[i][j] /= nodeIn1stLayer[i].size;
        if(nodeIn1stLayer[i].size > Const::th) {
//            assert(nodeIn1stLayer[i].size > Const::th);
            root->children[i] = new DumpyNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
            internal_size += root->children[i]->size;
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
//        if(nodeIn1stLayer[i].pid == -1) {
//            root->children[i] = new DumpyNode(1, nodeIn1stLayer[i].size, i);
//            root->children[i]->generateSaxAndCardIn1stLayer(i);
//            internal_size += root->children[i]->size;
//        }
//        else{
//            int pid = nodeIn1stLayer[i].pid;
//            root->children[i] = childrenList[pid];
//            childrenList[pid]->size += nodeIn1stLayer[i].size;
//            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
//        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

    // put data offsets to internal nodes in 1st layer
    for(int i=0;i<Const::vertexNum;++i)
        if(root->children[i]!= nullptr)
            root->children[i]->offsets.reserve(nodeIn1stLayer[i].size);
    for(int i=0;i<series_num;++i){
        int nav_id = navids[i];
        root->children[nav_id]->offsets.push_back(i);
    }
    Const::logPrint("Start fuzzy 1st layer.");

    unordered_map<DumpyNode*, NODE_RECORDER>FLNT;
    root->fuzzyFirstLayer(nodeIn1stLayer, navids, FLNT, *paa_mu_part_units);
    delete paa_mu_part_units;
    Const::logPrint("1st layer fuzzy number is " + to_string(fuzzy_num));


    thread IO(materialize1stLayerFuzzy, datafn, root, navids, Const::fuzzyidxfn, &FLNT);

    int finished_series = 0, finished_percent = 0;
    Const::logPrint("start to grow the index structure");
    unordered_map<int, unordered_map<DumpyNode*, NODE_RECORDER>* >DPNT;
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size > Const::th) {
            DumpyNode * child = root->children[i];
            // ensure the offset and series index list in DPLT is consistent with FLNT(file)
            sort(child->offsets.begin() + nodeIn1stLayer[i].size, child->offsets.end());
            auto*navigating_tbl = new unordered_map<DumpyNode*, NODE_RECORDER>();
            navigating_tbl->emplace(root->children[i], NODE_RECORDER(nodeIn1stLayer[i].size, root->children[i]));
            DPNT.emplace(i, navigating_tbl);
            root->children[i]->growIndexFuzzy(*navigating_tbl, g);

            finished_series += nodeIn1stLayer[i].size;
            double percent = (double)finished_series / (double)internal_size;
            if(percent >= (finished_percent+1) * 0.1) {
                Const::logPrint(to_string(percent * 100) + "% internal series in the 1st layer has been processed.");
                finished_percent = percent * 10;
            }
        }
    }
    Const::logPrint("build index skeleton finished.");

    IO.join();
    // navids have been deleted in the process of materializing the 1st layer

    Const::logPrint("Start materialize internal nodes in the 1st layer");
    finished_series = 0; finished_percent = 0;
    for(int i=0;i<Const::vertexNum;++i) {
        if (nodeIn1stLayer[i].size > Const::th) {
            materializeInterNodeFuzzy(root->children[i], saxes, nodeIn1stLayer[i].size, *DPNT[i]);

            finished_series += nodeIn1stLayer[i].size;
            double percent = (double)finished_series / (double)internal_size;
            if(percent >= (finished_percent+1) * 0.1) {
                Const::logPrint(to_string(percent * 100) + "% internal series in the 1st layer has been written to disk.");
                finished_percent = percent * 10;
            }
        }
    }

    Const::logPrint("Fuzzy series number is " + to_string(fuzzy_num));
    Const::logPrint("build index successfully!");

    return root;
}

void DumpyNode::growIndexFuzzy(unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl, vector<vector<int>> *g) {
    if(size <= Const::th)   return;

    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    NODE_RECORDER& parent_recorder = navigating_tbl[this];
//    PAA_INFO* paa = statPaa();
//    chooseSegment(paa, chosen_num);
    determineSegments();

    // statistic children information in order to partition
    partUnit nodes[power_2[chosen_num]];
    for(int i=0;i<power_2[chosen_num];++i)
        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
    vector<vector<int>>node_offsets(power_2[chosen_num], vector<int>());
    vector<vector<int>>series_index_list(power_2[chosen_num], vector<int>());
    vector<int> node_actual_size(power_2[chosen_num],0);

    for(int i=0;i<parent_recorder.actual_size;++i){
        int new_id = SaxUtil::extendSax(DumpyNode::saxes + (long)offsets[i] * Const::segmentNum, bits_cardinality,chosenSegments);
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
        series_index_list[new_id].push_back(parent_recorder.series_index_list[i]);
    }

    bool isLeaf[power_2[chosen_num]];
    for(int i=0;i<power_2[chosen_num];++i){
        isLeaf[i] = false;
        if(nodes[i].size <= Const::th)  isLeaf[i] = true;
        node_actual_size[i] = nodes[i].size;
    }

    // split fuzzy series
    for(int i=parent_recorder.actual_size;i<size;++i){
//        if(offsets[i] == 10294926)
//            cout << i <<endl;
        int new_id = SaxUtil::extendSax(DumpyNode::saxes + (long)offsets[i] * Const::segmentNum, bits_cardinality,chosenSegments, sax);
        // if the target node needn't have split, then the fuzzy series should not be inserted into it which will result in unnecessary split
        if(isLeaf[new_id] && nodes[new_id].size >= Const::th)   continue;
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
        series_index_list[new_id].push_back(parent_recorder.series_index_list[i]);
    }
    if(layer > 1) vector<int>().swap(offsets);

//    int partNum = partition(nodes);
    int partNum = partition(nodes, chosen_num);
    // build rest data node if any
//    for(auto &node:nodes)
//        if(node.size <= Const::th && node.pid == -1)
//            node.pid = ++partNum;


    DumpyNode* leafChildrenList[partNum];
    for(int i=0;i<partNum;++i) {
        leafChildrenList[i] = new DumpyNode(this, i);
        navigating_tbl.emplace(leafChildrenList[i], NODE_RECORDER());
    }
    children.resize(power_2[chosen_num]);
    for(int i=0;i<power_2[chosen_num];++i){
        if(nodes[i].size <= 0)  continue;
        else if(nodes[i].size > Const::th || nodes[i].pid == -1) {
            children[i] = new DumpyNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            children[i]->offsets.resize(nodes[i].size);
            copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
            navigating_tbl.emplace(children[i], NODE_RECORDER(node_actual_size[i], series_index_list[i]));
        }
        else{
            int _pid = nodes[i].pid;
            children[i] = leafChildrenList[_pid];
            leafChildrenList[_pid]->size += nodes[i].size;
            generateSaxAndCardinality4LeafNode(children[i], i);
            navigating_tbl[children[i]].actual_size += node_actual_size[i];
            for(int j=node_actual_size[i];j<nodes[i].size;++j)
                navigating_tbl[children[i]].series_index_list.push_back(series_index_list[i][j]);
        }
    }


    // note that for internal nodes, tbl stores all series index list inside while for leaf nodes, only fuzzy series index list are in the tbl
    fuzzy(nodes, node_actual_size, node_offsets, series_index_list, chosen_num, navigating_tbl);

    vector<vector<int>>().swap(node_offsets);
    vector<vector<int>>().swap(series_index_list);
    vector<int>().swap(node_actual_size);
    navigating_tbl.erase(this);

    for(auto &child: children){
        if(child!= nullptr){
            if(child->size > Const::th)
                child->growIndexFuzzy(navigating_tbl, g);
            else{   // this may be executed many times for the same node, but no problem
                if(navigating_tbl[child].series_index_list.empty())
                    navigating_tbl.erase(child);
                else{
                    navigating_tbl[child].actual_size = 0;
                    sort(navigating_tbl[child].series_index_list.begin(),  navigating_tbl[child].series_index_list.end());
                }
            }
        }
    }
}

void DumpyNode::fuzzy(partUnit* part_units, vector<int>& actual_sizes,
                      vector<vector<int>>&node_offsets, vector<vector<int>>&series_index_list,
                      int chosen_num, unordered_map<DumpyNode*, NODE_RECORDER>& navigating_tbl) const{
    for(int i=0;i<power_2[chosen_num];++i){
        fuzzySeriesInPartUnit(part_units, actual_sizes[i],chosen_num, node_offsets[i], series_index_list[i], navigating_tbl,i);
    }
}

struct CAND{
    int id;
    double score;

    CAND(int _i, double _score){ id = _i; score = _score;}

    static bool order(CAND &a, CAND &b){return a.score < b.score;}
};


void DumpyNode::fuzzySeriesInPartUnit(partUnit *part_units, int actual_size, int chosen_num, vector<int> &node_offsets,
                                      vector<int> &series_index_list,
                                      unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl, int _id) const{
    float *paa;
    double range, lb, ub;
    DumpyNode*temp_node;
    vector<CAND>candidates; // max is chosen num
    int new_id, seg, sax_symbol_last_bit, sax_symbol;

    for(int j=0;j<actual_size;++j){
        paa = DumpyNode::paas + (long)node_offsets[j] * Const::segmentNum;
        for(int i=0;i<chosenSegments.size();++i){
            seg = chosenSegments[i];
            new_id = _id ^ mask[Const::segmentNum - chosen_num + i];
            sax_symbol_last_bit = (_id & mask[Const::segmentNum - chosen_num + i]) > 0 ? 1 : 0;
            sax_symbol = (sax[seg] << 1) + sax_symbol_last_bit;
            if(part_units[new_id].size > 0 && children[_id] != children[new_id]){
                temp_node = children[new_id];
                // if target leaf node full, continue
                if(temp_node->partition_id != -1 && temp_node->size >= Const::th)   continue;

                // calculate the range of new partition unit in this segment
                // if any bound is infinity, then select the neighbor range as range
                if( sax_symbol == (power_2[bits_cardinality[seg] + 1] - 1) )
                    SaxUtil::getValueRange(sax[seg] << 1, bits_cardinality[seg] + 1, &lb, &ub);
                else if(sax[seg] == 0 && sax_symbol_last_bit == 0)
                    SaxUtil::getValueRange(1, bits_cardinality[seg] + 1, &lb, &ub);
                else
                    SaxUtil::getValueRange(sax_symbol , bits_cardinality[seg] + 1, &lb, &ub);
                range = (ub -lb) * Const::fuzzy_f;
                // check whether this series can be put into the candidates list
                // Single precision paa may lead to conflicts with sax, so we need 1e-4 and abs. series are routed by sax to target node.
                if(sax_symbol == (power_2[bits_cardinality[seg] + 1] - 1)){
//                    assert(paa[seg] >= ub - 1e-4);
                    if(abs(paa[seg] - ub) <= range) candidates.emplace_back(new_id, abs(paa[seg] - ub));
                }else if(sax_symbol == 0){
//                    assert(lb >= paa[seg] - 1e-4);
                    if(abs(lb - paa[seg]) <= range)  candidates.emplace_back(new_id, abs(lb -paa[seg]));
                }else if(sax_symbol_last_bit == 0){
//                     TODO: check
//                    assert(ub >= paa[seg] - 1e-4);
                    if(abs(ub - paa[seg]) <= range) candidates.emplace_back(new_id, abs(ub - paa[seg]));
                }else{
//                    assert(paa[seg] >= lb - 1e-4);
                    if(abs(paa[seg] - lb) <= range)    candidates.emplace_back(new_id, abs(paa[seg] - lb));
                }
            }
        }
        sort(candidates.begin(), candidates.end(), CAND::order);
        int n= 0;
        // no need to update node_offsets and series_index_lists
        for(int i = 0; i< candidates.size() && n < Const::delta; ++i){
            temp_node = children[candidates[i].id];
            if(temp_node->partition_id != -1){   // temp node is a leaf node
                if(temp_node->size < Const::th){
                    temp_node->size++;
                    navigating_tbl[temp_node].series_index_list.push_back(series_index_list[j]);
                    ++n;
                }
            }else{  // temp node is an internal node
                assert(navigating_tbl[temp_node].actual_size > Const::th);
                temp_node->size++;
                temp_node->offsets.push_back(node_offsets[j]);
                navigating_tbl[temp_node].series_index_list.push_back(series_index_list[j]);
                ++n;
            }
        }
        candidates.clear();
        fuzzy_num +=n;
    }
}

void DumpyNode::fuzzySeriesInPartUnitInFirstLayer(partUnit *part_units, vector<int> &node_offsets, int _id,
                                                  unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl,
                                                  vector<vector<double >> &paa_mu_part_units) const{
    float *paa, range;
    DumpyNode* temp_node;
    int new_id;
    vector<CAND>candidates;
    for(int j=0; j < part_units[_id].size; ++j){
//        if(j == 31068)
//            cout << j << endl;
        paa = DumpyNode::paas + (long)node_offsets[j] * Const::segmentNum;
        for(int i=0;i<Const::segmentNum;++i){
            new_id = _id ^ mask[i];
            if(part_units[new_id].size > 0 && children[_id] != children[new_id]){
                if(part_units[new_id].size < Const::th && children[new_id]->size >= Const::th)   continue;
                if(paa[i] > 0){
                    range = paa_mu_part_units[_id][i] * Const::fuzzy_f_1st;
                    if(paa[i] <= range) candidates.emplace_back(new_id, paa[i]);
                }else{
                    range = (-paa_mu_part_units[_id][i]) * Const::fuzzy_f_1st;
                    if(-paa[i] <= range)    candidates.emplace_back(new_id, -paa[i]);
                }
            }
        }
        sort(candidates.begin(), candidates.end(), CAND::order);
        int n= 0;
        for(int i = 0; i< candidates.size() && n < Const::delta; ++i){
            temp_node = children[candidates[i].id];
            if(temp_node->partition_id != -1){   // temp node is a leaf node
                if(temp_node->size < Const::th){
                    temp_node->size++;
                    navigating_tbl[temp_node].series_index_list.push_back(node_offsets[j]);
                    ++n;
                }
            }else{  // temp node is an internal node
                temp_node->size++;
                temp_node->offsets.push_back(node_offsets[j]);
                navigating_tbl[temp_node].series_index_list.push_back(node_offsets[j]);
                ++n;
            }
        }
        candidates.clear();
        fuzzy_num += n;
    }//end loop series in part unit

}


void DumpyNode::fuzzyFirstLayer(partUnit *part_units, const int *nav_ids,
                                unordered_map<DumpyNode *, NODE_RECORDER> &navigating_tbl,
                                vector<vector<double>> &paa_mu_part_units) const {
    // data offsets are stored in internal node rather than leaf part units
    // node_offsets store data offsets of leaf part units
    unordered_map<int, vector<int> > node_offsets;
    for(int i=0;i<Const::vertexNum;++i){
        if(part_units[i].size <= Const::th)
            node_offsets[i].reserve(part_units[i].size);
    }
    for(int i=0;i<this->size; ++i){
//        if(i == 4492048){
//            cout << "here"<< endl;
//        }
        if(part_units[nav_ids[i]].size <= Const::th)
            node_offsets[nav_ids[i]].push_back(i);
    }
    for(int i=0;i<Const::vertexNum;++i){
        if(children[i] == nullptr)  continue;
//        if(i == 1020)
//            cout <<i << endl;
        if(this->children[i]->partition_id == -1)
            fuzzySeriesInPartUnitInFirstLayer(part_units, this->children[i]->offsets, i, navigating_tbl,
                                              paa_mu_part_units);
        else
            fuzzySeriesInPartUnitInFirstLayer(part_units, node_offsets[i], i, navigating_tbl, paa_mu_part_units);
    }
    node_offsets.clear();
    for(auto &iter: navigating_tbl){
        iter.second.actual_size = 0;
        sort(iter.second.series_index_list.begin(),  iter.second.series_index_list.end());
    }
}