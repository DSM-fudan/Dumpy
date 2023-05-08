//
// Created by Zeyu Wang on 2021/8/7.
//

#ifndef DUMPY_FILEUTIL_H
#define DUMPY_FILEUTIL_H
#include <string>
#include <vector>
#include "../DataStructures/PqItemSeries.h"

using namespace std;

class FileUtil {

public:
    static bool checkFileExists(const char *name);

    static void Getfilepath(const char *path, const char *filename, char *filepath);

    static bool checkDirClean(const char *path);

    static long getFileSize(const char *fname);

    static void getFiles(const string& path, vector<string>& files );

    static float *readSeries(FILE *f);

    static int FileRemove(const char *fname);

    static void mergeFiles(const string sources[], const string &dest, int num);

    static void deleteFiles(const string *fs, int num);

    static bool createDir(const char *path);

    static float* readQueries();
};


#endif //DUMPY_FILEUTIL_H
