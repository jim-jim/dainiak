#include <string>
#include <windows.h>

#define MAX_STACK	44444444
#define MAX_VARC	4444
#define FILTER		4
#define OUTLEN		31
typedef unsigned char sample;
typedef long long ngtype;
struct Ngram{
	ngtype a;
	ngtype b;
	byte len;
	byte count;
	char mnm;
	char mxm;
	int eq1, eq2;
	Ngram **next;
	Ngram *ng2;
	void print();
	char* StrA(char *buf);
};
struct Sol{
	int varc;
	int k;
	int nglen;
	int maxlen;
	sample *sol;
	Ngram *ng;
	void print();
};
int main(int argc, char **argv, char **envp);
int first_step(int k, int delta);
int next_step(Sol *sol);
bool solve_eq(int *eqX, int **eqs, int eqc, int &solc);
void normalize_eqX(int *eqX, int eqvarc, int **eqs, int eqc);
int get_codes(Sol *sol);
int _get_codes(Sol *sol, Ngram *ngram, int pos, int mnm, int mxm);

#define _CRT_SECURE_NO_WARNINGS
#define LOG_PATH	"log.txt"
//#define OUTPUT
void msg(int x);
void msg(char *x);
void printline();
void print_eq(int *eq, int mode);
void print_eqs(int **eqs, int eqc);
void write(char *s, int len);
void print_ng(Ngram *ng, int count);
void print_ints(int *data, int count);