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
#include <string.h>
#include<iostream>
#include <gmp.h>
using namespace std;

bool print = false;
bool print_time = true;

long long int MOD = 100000000;
int L = 5;


/*
 * Time profiling
 */
long timediff(clock_t t1, clock_t t2){
    long elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000000;
    return elapsed;
}


/*
 * Big integers
 */
typedef struct bigint{
	long long int d[5];

	bigint& operator =(const bigint& a){
		for(int i = 0; i < L; i++) d[i] = a.d[i];
		return *this;
	}
} bigint;

bigint zero, one, two;
bigint power10[63];

void read(bigint* x){
	string s;
	*x = zero;
	cin >> s;
	int where = 0;
	int pow10 = 1;
	for(int i = 0; i < s.size(); i++){
		int digit = s[s.size()-1-i] - '0';
		x->d[where] += pow10 * digit;
		if(i % 8 == 7){
			where++;
			pow10 = 1;
		}
		else pow10 *= 10;
	}
}

void write(const bigint x){
	bool leadingzeros = false;
	for(int i = L-1; i >= 0; i--){
		if(leadingzeros) printf("%08d", (int)x.d[i]);
		else if(x.d[i] != 0 || i == 0){
			printf("%d", (int)x.d[i]);
			leadingzeros = true;
		}
	}
}

void writeln(const bigint x){
	write(x);
	printf("\n");
}

int size(const bigint x){
	for(int i = L-1; i >= 0; i--){
		if(x.d[i] != 0){
			int res = 0;
			int pow10 = 1;
			while(x.d[i] >= pow10){
				pow10 *= 10;
				res++;
			}
			return res + 8*i;
		}
	}
	return 0;
}

bool operator==(const bigint x, const bigint y){
	for(int i = 0; i < L; i++){
		if(x.d[i] != y.d[i]) return false;
	}
	return true;
}

bool operator<(const bigint x, const bigint y){
	for(int i = L-1; i >= 0; i--){
		if(x.d[i] < y.d[i]) return true;
		else if(x.d[i] > y.d[i]) return false;
	}
	return false;
}

bool operator>(const bigint x, const bigint y){
	return (y < x);
}

bigint min_bigint(bigint x, bigint y){
	if (x<y) return x;
	return y;
}

bool operator<=(const bigint x, const bigint y){
	return !(y < x);
}

bool operator>=(const bigint x, const bigint y){
	return !(x < y);
}

bigint operator+(const bigint x, const bigint y){
	bigint z;
	long long int rest = 0;
	for(int i = 0; i < L; i++){
		z.d[i] = rest + x.d[i] + y.d[i];
		rest = 0;
		while(z.d[i] > MOD){
			z.d[i] -= MOD;
			rest++;
		}
	}
	return z;
}

bigint operator-(const bigint x, const bigint y)
{
	bigint z;
	long long int rest = 0;
	for(int i = 0; i < L; i++)
	{
		z.d[i] = rest + x.d[i] - y.d[i];
		rest = 0;
		while(z.d[i] < 0)
		{
			z.d[i] += MOD;
			rest--;
		}
	}
	if(rest < 0) fprintf(stderr, "NEGATIVE NUMBER ERROR\n");
	return z;
}

bigint operator*(const bigint x, const bigint y)
{
	bigint z;
	long long int rest = 0;
	for(int i = 0; i < L; i++)
	{
		z.d[i] = rest;
		for(int j = 0; j <= i; j++)
		{
			z.d[i] += x.d[j] * y.d[i-j];
		}
		rest = z.d[i] / MOD;
		z.d[i] = z.d[i] % MOD;
	}
	return z;
}

void div10(bigint* x)
{
	for(int i = 0; i < L; i++)
	{
		int r = x->d[i] % 10;
		if(i > 0) x->d[i-1] += r*10000000;
		x->d[i] /= 10;
	}
}

bigint operator/(const bigint x, const bigint y)
{
	bigint r = x, q = zero;
	if(y == zero)
	{
		fprintf(stderr, "DIVISION BY ZERO ERROR\n");
		return x;
	}
	int t = size(x) - size(y);
	bigint start = power10[t] * y;
	while(t >= 0)
	{
		while(start <= r)
		{
			r = r-start;
			q = q+power10[t];
		}
		t--;
		div10(&start);
	}
	return q;
}

bigint operator%(const bigint x, const bigint y)
{
	bigint r = x, q = zero;
	if(y == zero)
	{
		fprintf(stderr, "MODULO ZERO ERROR\n");
		return x;
	}
	int t = size(x) - size(y);
	bigint start = power10[t] * y;
	while(t >= 0)
	{
		while(start <= r)
		{
			r = r-start;
			q = q+power10[t];
		}
		t--;
		div10(&start);
	}
	return r;
}

bool even(const bigint x){
	if(x.d[0] & 1) return false;
	else return true;
}

bigint power(bigint x, bigint n){
	if(n == zero) return one;
	if(n == one) return x;
	bigint y = power(x, n/two);
	if(even(n)) return y*y;
	else return (y*y)*x;
}

bigint powermod(bigint x, bigint n, bigint mod){
	if(n == zero) return one;
	if(n == one) return x;
	bigint y = powermod(x, n/two, mod);
	if(even(n)) return (y*y)%mod;
	else return (((y*y)%mod)*x)%mod;
}

void initialize(){
	for(int i = 0; i < L; i++) zero.d[i] = 0;
	one = zero;
	one.d[0] = 1;
	two = zero;
	two.d[0] = 2;
	for(int i = 0; i < L; i++){
		int p = 1;
		for(int j = 0; j < 8; j++){			power10[8*i+j] = zero;
			power10[8*i+j].d[i] = p;
			p *= 10;
		}
	}
}


// Random bigint between 2 and 1000000
void random_bigint(bigint* x, bigint N){
	*x = zero;
	bigint ma = zero;
	ma.d[0] = 99999999;
	bigint mi = min_bigint(N, ma);
	x->d[0] = (rand() % (mi.d[0]-1)) + 1;
}

bigint gcd(bigint a, bigint b){
	while (!(b==zero)){
		bigint tmp = b;
		b = a % b;
		a = tmp;
	}
	return a;
}

/*
 * Primality test
 */
bool miller_rabin(bigint N){
	if (N==one) return false;
	if (even(N)){
		if (N == two) return true;
		return false;
	}
	if (print) printf("N is odd\n");
	//N − 1 = t*2^s with t odd
	int s = 0;
	bigint t = N-one;
	while (even(t)){
		t = t/two;
		s++;
	}
	if (print){
		printf("s = %d\n", s);
		writeln(t);
	}
	// Loop
	bigint u;
	for (int i = 0; i < 20; i++){
		if (print) printf("Test number %d\n", i);
		bool prob_prime = false;
		bigint a;
		random_bigint(&a, N);
		if (print) {
			printf("Random number a = ");
			writeln(a);
		}
		u = powermod(a, t, N);
		if (print) {
			printf("a^t mod N = ");
			writeln(u);
		}
		if (u==one) prob_prime = true;
		if (u==N-one) prob_prime = true;
		for (int j=0; j < s-1; j++){
		 	u = powermod(u, two, N);
		 	if (print) {
		 		printf("u^2 mod N = ");
		 		writeln(u);
		 	}
			if (u==N-one) prob_prime = true;
		}
		if (print) printf("prob_prime : %s\n", prob_prime ? "true" : "false");
		if (prob_prime == false) return false;
	}
	return true;
}

/*
 * Factorisation
 */
pair<bigint, bigint> pollard_rho(bigint N, bigint param){
    bigint x = two;
    bigint temp = zero;
    bigint y = two;
    bigint d = one;
    while (d == one){
        x = (power(x, two) + param) % N;
        temp = (power(y, two) + param) % N;
        y = (power(temp, two) + param) % N;
        if (x > y) d = gcd(x-y, N);
        else d = gcd(y-x, N);
    }
    if (d == N) return pollard_rho(N, param + one);
    return make_pair(d, N/d);
}

vector<bigint> factorize(bigint N){
	vector<bigint> factors;
	vector<bigint> tofactor;
	while (even(N)){
		N = N/two;
		factors.push_back(two);
	}
	if (print){
		printf("N without all the factors 2 = ");
		writeln(N);
	}
	if (miller_rabin(N)) factors.push_back(N);
	else {
		tofactor.push_back(N);
		while(!tofactor.empty()){
			vector<bigint> newtofactor;
			for(int i = 0; i < tofactor.size(); i++){
				if (print){
					printf("To factor : %d/%lu", i, tofactor.size());
					writeln(tofactor[i]);
				}
				pair<bigint, bigint> p = pollard_rho(tofactor[i], one);
				if (miller_rabin(p.first)) factors.push_back(p.first);
				else newtofactor.push_back(p.first);
				if (miller_rabin(p.second)) factors.push_back(p.second);
				else newtofactor.push_back(p.second);
			}
			tofactor = newtofactor;
		}
	}
	sort(factors.begin(), factors.end());
	/*printf("All the factors \n");
	for (int i = 0; i < factors.size(); i++){
		writeln(factors[i]);
	}*/
	return factors;
}

int main(int argc, char *argv[]){
    initialize(); // DO NOT REMOVE!

	bigint x;
	read(&x);
	while (!(x==zero)){
		//clock_t start = clock();
		vector<bigint> v = factorize(x);
		/*if (print_time) {
			clock_t end = clock();
			long elapsed = timediff(start, end);
			printf("Took %ld microseconds\n", elapsed);
		}*/

		bigint current = v[0];
		int exp = 1;
		bool newfact = false;
		for (int i = 1; i < v.size(); i++){
			if (v[i]==current) exp++;
			else{
				write(current);
				printf("^%d ", exp);
				current = v[i];
				exp = 1;
			}
			if (i == v.size()-1){
				write(current);
				printf("^%d\n", exp);
			}
		}
		read(&x);
	}

    return EXIT_SUCCESS;
}
