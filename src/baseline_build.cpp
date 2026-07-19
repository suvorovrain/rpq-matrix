#include <sys/stat.h>
#include <dirent.h>
#include <cstring>
#include <string>
#include <vector>
#include <omp.h>

extern "C" {
#include "baseline/matrix.h"
}

int main(int argc, char **argv) {
    uint logside, l, n;
    uint *coords;
    FILE *f;
    struct stat st;
    matrix M;
    char sys[1024], fname[1024];
    DIR *dir;
    struct dirent *entry;

    if (argc < 4) {
        fprintf(stderr, "Usage: %s <pairs dir> <matrices dir> <side>\n", argv[0]);
        exit(1);
    }

    sprintf(sys, "mkdir -p %s", argv[2]); 
    system(sys);

    dir = opendir(argv[1]);
    if (dir == NULL) {
        fprintf(stderr, "Error: cannot open directory %s\n", argv[1]);
        exit(1);
    }
    uint side = atoi(argv[3]);

    std::vector<std::string> files;
    while ((entry = readdir(dir))) {
        l = strlen(entry->d_name);
        if (strcmp(".pairs", entry->d_name + l - strlen(".pairs"))) continue;
        files.push_back(entry->d_name);
    }
    closedir(dir);

    #pragma omp parallel for private(fname, coords, n, f, M, st) schedule(dynamic)
    for (size_t i = 0; i < files.size(); i++) {
        sprintf(fname, "%s/%s", argv[1], files[i].c_str());
        if (stat(fname, &st) != 0) {
            fprintf(stderr, "Error: cannot stat file %s\n", fname);
            continue;
        }

        n = st.st_size / sizeof(uint);
        coords = (uint*)myalloc(n * sizeof(uint));

        #pragma omp critical
        {
            f = fopen(fname, "r");
            if (f == NULL) {
                fprintf(stderr, "Error: cannot open file %s\n", fname);
                myfree(coords);
                continue;
            }
            fread(coords, sizeof(uint), n, f);
            fclose(f);
        }

        printf("Generating matrix for %s...", fname); 
        fflush(stdout);

        M = matCreate(side + 1, side + 1, n / 2, coords);

        #pragma omp critical
        {
            sprintf(fname, "%s/%s", argv[2], files[i].c_str());
            strcpy(fname + strlen(fname) - strlen(".pairs"), ".mat");
            f = fopen(fname, "w");
            if (f == NULL) {
                fprintf(stderr, "Error: cannot create file %s\n", fname);
                matDestroy(M);
                myfree(coords);
                continue;
            }
            matSave(M, f);
            fclose(f);
            matDestroy(M);
            myfree(coords);
            printf("written\n");
        }
    }

    return 0;
}