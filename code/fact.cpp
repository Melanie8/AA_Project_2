/*
 * Written by :
 * Melanie Sedda <sedda@kth.se>
 * Catarina Vaz  <acvaz@kth.se>
 * October 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdbool.h>
#include<iostream>
#include <gmp.h>
using namespace std;

int max_len = 29;
bool print = true;

long timediff(clock_t t1, clock_t t2) {
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000000;
    return elapsed;
}

mpz_t pollard_rho(mpz_t N){
    vector<mpz_t> x;
    x.push_back[2];
    int i = 0;
    while (true){
        x.push_back(x[i]*x[i])
        i++;
    }
}

int main(int argc, char *argv[]) {
    clock_t start = clock();
    srand(time(NULL));
    int i, j, k,r;
    char *s = (char *)malloc(max_len*sizeof(char));

    mpz_t N;
    mpz_init(N);
    /* Get arguments */
    scanf("%s",s);
    mpz_set_str(N, s, 10);
    if (print) {
        cout<<"Value: ";
        mpz_out_str(stdout, 10, N); //Stream, numerical base, var
        cout<<endl;
    }
    while(strncmp(s, "0", max_len)) {
        printf("Do stuff\n");

        // Get next number
        scanf("%s",s);
        mpz_set_str(N, s, 10);
        if (print) {
            cout<<"Value :";
            mpz_out_str(stdout, 10, N); //Stream, numerical base, var
            cout<<endl;
        }
    }

    return EXIT_SUCCESS;
}
