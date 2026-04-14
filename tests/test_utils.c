#include <assert.h>
#include "../utils.h"

int main(){
    // test of read_from_file and write_to_file_iter
    char* name = "test_read_write_from_file";
    int n = 1000;
    
    double* res = malloc(n*sizeof(double));

    for (int i = 0; i < n; i++){
        res[i] = (double) rand() / (double) RAND_MAX;
    }

    write_to_file_iter(res, n, name);

    double* read = read_from_file(name, n);

    if (n != 0){
        for(int i = 0; i < n; i++){
            assert(res[i] - read[i] <= 1);
        }
    }
    return 0;
}