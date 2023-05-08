//
// Created by Zeyu Wang on 2021/8/7.
//

#include "../../include/Utils/FileUtil.h"
#include "../../include/DataStructures/TimeSeries.h"
#include <fstream>
#include <unistd.h>
#include <unordered_set>
#include <cstdio>
#include <filesystem>
#include <dirent.h>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

bool FileUtil::checkFileExists(const char *name) {
    ifstream f(name);
    return f.good();
}

int FileUtil::FileRemove(const char* fname)
{
    return remove(fname);
}

void FileUtil::deleteFiles(const string fs[], int num){
    for(int i=0;i<num;++i)
        FileRemove(fs[i].c_str());
}

long FileUtil::getFileSize(const char* fname)
{
    struct stat statbuf{};
    if(stat(fname,&statbuf)==0)
        return statbuf.st_size;
    return -1;
}

void FileUtil::Getfilepath(const char *path, const char *filename,  char *filepath)
{
    strcpy(filepath, path);
    if(filepath[strlen(path) - 1] != '/')
        strcat(filepath, "/");
    strcat(filepath, filename);
}

bool FileUtil::createDir(const char * path){
    return !mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);///home/newdir
}

bool FileUtil::checkDirClean(const char* path)
{
    if(!checkFileExists(path))
        createDir(path);
    DIR *dir;
    struct dirent *dirinfo;
    struct stat statbuf{};
    char filepath[256] = {0};
    lstat(path, &statbuf);

    if (S_ISREG(statbuf.st_mode))//判断是否是常规文件
        {
        remove(path);
        }
    else if (S_ISDIR(statbuf.st_mode))//判断是否是目录
        {
        if ((dir = opendir(path)) == nullptr)
            return true;
        while ((dirinfo = readdir(dir)) != nullptr)
        {
            Getfilepath(path, dirinfo->d_name, filepath);
            if (strcmp(dirinfo->d_name, ".") == 0 || strcmp(dirinfo->d_name, "..") == 0)//判断是否是特殊目录
                continue;
            checkDirClean(filepath);
            rmdir(filepath);
        }
        closedir(dir);
        }
    return false;
}

void FileUtil::mergeFiles(const string sources[], const string& dest, int num ){
    FILE * out = fopen(dest.c_str(), "wb");
    for(int i=0;i<num;++i){
        long size = getFileSize(sources[i].c_str());
        char * tmp = new char[size];
        FILE *in = fopen(sources[i].c_str(), "rb");
        fread(tmp, 1, size, in);
        fclose(in);
        fwrite(tmp, 1, size, out);
        free(tmp);
    }
    fclose(out);
}

float* FileUtil::readQueries(){
    auto queries = new float[Const::tsLength * Const::query_num];
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    for(int i=0;i<Const::query_num;++i){
        fread(queries + i * Const::tsLength, sizeof(float ), Const::tsLength, f);
    }
    fclose(f);
    return queries;
}

void FileUtil::getFiles(const string& path, vector<string>& files )
{
    files.clear();
    for(auto &entry:filesystem::directory_iterator(path))
        if(entry.is_directory())    continue;
        else files.push_back(entry.path());
}

float * FileUtil::readSeries(FILE* f){
    auto *ts = new float[Const::tsLength];
    fread(ts, sizeof(float), Const::tsLength, f);
    return ts;
}
