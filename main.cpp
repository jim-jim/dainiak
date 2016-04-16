#include "main.h"

int **VarEqs, Varc, VarX;
long long *Stack, *StackEnd;
sample *Solution, *SolutionEnd, *CurrSol;
bool MaxStackExc = false;
int BestLen = 0;
int *AcMin, *AcMax;
char *CurrCode;
int OutLen;
template<typename T>
inline T* alloc(int count){
	int size = (count * sizeof(T) + sizeof(long long) - 1) / sizeof(long long);
	T* res = (T*)Stack;
	Stack += size;
	return res;
}
int main(int argc, char **argv, char **envp){
	int k = (argc > 1 ? atoi(argv[1]) : 4);
	int delta = (argc > 2 ? atoi(argv[2]) : -2);
	if (k < 1 || k > 44 || 4 * k + delta < 1 || delta > 4 * k) return 0;
	DeleteFileA(LOG_PATH);
	CurrSol = new sample[MAX_VARC];
	VarEqs = new int*[MAX_VARC];
	Stack = new long long [MAX_STACK];
	StackEnd = Stack + MAX_STACK;
	AcMin = new int[MAX_VARC];
	AcMax = new int[MAX_VARC];
	CurrCode = new char[MAX_VARC];
	first_step(k, delta);
	if (MaxStackExc) msg("MAX");
	msg(BestLen);
	return 0;
}
int first_step(int k, int delta){
	//#1=k, #1+#2+#3+#4=2k, #1+2#2+3#3+4#4=4k+delta
	//#2+#3+#4=k, #3+2#4=k+delta
	//#4=#2+delta, #3=k-#2-#4=k-2#2-delta
	long long *stack = Stack;
	int n = 4 * k + delta;
	int ac = 0;
	AcMin[0] = AcMax[0] = n;
	for (int i = 1, im = n / 2; i <= im; i++){
		if (i == 4) ac += 2;
		//if (i == 7) ac += 2;
		AcMin[i] = AcMin[n - i] = -ac;
		AcMax[i] = AcMax[n - i] = ac;
	}
	for (int i = 0; i < 4; i++){
		AcMin[n + i] = AcMin[i];
		AcMax[n + i] = AcMax[i];
	}
	Sol sol;
	sol.k = k;
	sol.nglen = 0;
	sol.varc = 4;
	sol.sol = alloc<sample>(4);
	sol.sol[0] = k;
	sol.ng = alloc<Ngram>(4);
	sol.maxlen = min(n - 3, sizeof(ngtype) * 4 - 1);
	OutLen = min(OUTLEN, sol.maxlen);
	Ngram **solnext = alloc<Ngram*>(5);
	for (int j = 0; j < 4; j++) solnext[j] = sol.ng + j;
	solnext[4] = NULL;
	for (int i = 0; i < 4; i++){
		Ngram *ngram = sol.ng + i;
		ngram->count = 2;
		ngram->eq1 = ngram->eq2 = -1;
		ngram->len = 0;
		ngram->mnm = ngram->mxm = i + 1;
		ngram->a = i;
		ngram->b = (1 << (i+1)) - 1;
		ngram->next = solnext;
		ngram->ng2 = NULL;
	}
	int count = 0;
	for (int x2 = max(0,-delta), x2m = (k-delta)/2; x2 <= x2m; x2++){
		sol.sol[1] = x2;
		sol.sol[2] = k - 2 * x2 - delta;
		sol.sol[3] = x2 + delta;
		count += next_step(&sol);
	}
	Stack = stack;
	return count;
}
int next_step(Sol *sol){
	static int table[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };
	long long *stack = Stack;
	int count = 0;
	Sol sol2;
	sol2.k = sol->k;
	sol2.nglen = sol->nglen + 1;
	sol2.maxlen = sol->maxlen;
	sol2.ng = alloc<Ngram>(sol->varc * 4);
	int **eqs = alloc<int*>(sol->varc * 4);
	int *eqa = alloc<int>(sol->varc * 4 * 6);
	if (Stack >= StackEnd) goto end;
	int eqc = 0;
	int eq0v = 2 * sol->k;
	int eqvarc = 0;
	int ac = sol2.nglen + 1;
	Varc = 0;
	for (int i = 0; i < sol->varc; i++){
		if (sol->sol[i] == 0) continue;
		Ngram *ngram = sol->ng + i;
		int a_first = (ngram->a & 3) + 1;
		if (ngram->len == sol->nglen){
			bool posit = (ngram->count % 4 == 0);
			int jm = (posit ? FILTER - ngram->mxm : FILTER + ngram->mnm);
			if (jm == 0) goto end;
			if (jm > 4) jm = 4;
			int eq1 = -1;
			int eqi1 = 0;
			ngtype mask = ((ngtype)1 << ngram->count) - 1;
			for (int j = 1; j <= jm; j++){
				ngtype _a = (ngram->a | ((ngtype)(j - 1) << (ngram->count))) >> 2;
				int eq2 = -1;
				bool was = false;
				Ngram **next = ngram->next;
				Ngram *node;
				while (*next){
					if (((*next)->a & mask) == _a){
						if ((*next)->eq2 >= 0){ node = *next; eq2 = node->eq2; was = true; break; }
						else if (eq2 < 0){
							eq2 = eqc;
							eqs[eq2] = eqa + 6 * eqc;
							eqc++;
							eqs[eq2][1] = eqs[eq2][0] = 0;
							node = *next;
							node->eq2 = eq2;
						}
						eqs[eq2][1] += sol->sol[(int)((*next)-sol->ng)];
					}
					next++;
				}
				if (eq2 < 0 || eq2 == INT_MAX) continue;
				if (eqs[eq2][1] == 0){ eqc--; node->eq2 = INT_MAX; continue; }
				if (j == 1 && !was){
					if (posit) eq0v -= eqs[eq2][1];
					else eq0v += eqs[eq2][1];
				}
				sol2.ng[eqvarc] = *ngram;
				Ngram *_ngram = &sol2.ng[eqvarc];
				if (posit){
					_ngram->mnm = min(j,_ngram->mnm+j);
					_ngram->mxm = max(j,_ngram->mxm+j);
					_ngram->a |= (ngtype)(j - 1) << (_ngram->count);
					_ngram->b |= (ngtype)((1 << j) - 1) << (_ngram->len + a_first);
				}
				else{
					_ngram->mnm = min(-j,_ngram->mnm-j);
					_ngram->mxm = max(-j,_ngram->mxm-j);
					_ngram->a |= (ngtype)(j - 1) << (_ngram->count);
					_ngram->b |= (ngtype)1 << (_ngram->len + a_first + j);
				}
				_ngram->len += j;
				_ngram->count += 2;
				if (eq1 < 0){
					eq1 = eqc;
					eqs[eq1] = eqa + 6 * eqc;
					eqc++;
					eq0v -= (table[((1 << a_first) - 1) & (_ngram->b ^ (_ngram->b >> ac))] + (j == 1 ? (posit ? -1 : 1) : 0)) * sol->sol[i];
					eqs[eq1][eqi1++] = 0;
					eqs[eq1][eqi1++] = sol->sol[i];
					ngram->eq1 = eq1;
				}
				eqs[eq1][0]++;
				eqs[eq1][eqi1++] = eqvarc;
				eqs[eq2][++(eqs[eq2][0]) + 1] = eqvarc;
				VarEqs[2 * eqvarc] = eqs[eq1];
				VarEqs[2 * eqvarc + 1] = eqs[eq2];
				eqvarc++;
			}
			if (eq1 < 0) goto end;
		}
		else{
			eq0v -= (table[((1 << a_first) - 1) & (ngram->b ^ (ngram->b >> ac))]) * sol->sol[i];
			Varc++;
		}
	}
	if (eq0v < AcMin[ac] || eq0v > AcMax[ac]) goto end;
	Varc += eqvarc;
	sol2.varc = eqvarc;
	for (int i = 0; i < sol->varc; i++){
		Ngram *ngram = sol->ng + i;
		if (sol->sol[i] == 0) continue;
		if (ngram->len != sol->nglen){
			sol2.ng[sol2.varc] = *ngram;
			CurrSol[sol2.varc] = sol->sol[i];
			ngram->ng2 = sol2.ng + sol2.varc;
			sol2.varc++;
		}
	}
	Ngram **sol2next = alloc<Ngram*>(Varc * 44);
	if (Stack >= StackEnd) goto end;
	Ngram **sol2ni = sol2next;
	for (int i = 0; i < sol2.varc; i++){
		Ngram *ngram = sol2.ng + i;
		ngram->eq1 = ngram->eq2 = -1;
		ngram->ng2 = NULL;
		Ngram **next = ngram->next;
		ngtype mask = ((ngtype)1 << (ngram->count - 2)) - 1;
		ngtype _a = (ngram->a >> 2);
		int j = 0;
		while (*next){
			if ((*next)->ng2){
				if (((*next)->a & mask) == _a){
					sol2ni[j] = (*next)->ng2;
					j++;
				}
			}
			else if ((*next)->eq1 >= 0){
				int *eq = eqs[(*next)->eq1];
				for (int jj = 0; jj < eq[0]; jj++){
					Ngram *_ng = sol2.ng + eq[jj + 2];
					if ((_ng->a & mask) == _a){
						sol2ni[j] = _ng;
						j++;
					}
				}
			}
			next++;
		}
		sol2ni[j] = NULL;
		ngram->next = sol2ni;
		sol2ni += j + 1;
	}
	Stack = 1 + (long long *)sol2ni;
	int *eqX = alloc<int>(eqvarc + 2);
	if (Stack >= StackEnd) goto end;
	memset(eqX, 0, (eqvarc + 2) * sizeof(int));
	eqX[1] = 2 * sol->k - AcMin[ac + 1];
	VarX = AcMax[ac + 1] - AcMin[ac + 1];
	for (int i = 0; i < sol2.varc; i++){
		Ngram *ngram = sol2.ng + i;
		int a_first = (ngram->a & 3) + 1;
		ngtype _b = ngram->b;
		if (ngram->len == sol2.nglen){
			bool posit = (ngram->count % 4 == 0);
			if (posit) _b |= ((ngtype)3 << (ngram->len + a_first));
			ngtype _a = (ngram->a >> 2);
			ngtype mask = ((ngtype)1 << ngram->count) - 1;
			Ngram **next = ngram->next;
			while (*next){
				if (((*next)->a & mask) == _a){
					if ((*next)->eq2 >= 0) break;
					(*next)->eq2 = 0;
					int vi = (int)((*next)-sol2.ng);
					if (posit){
						if (vi < eqvarc) eqX[2 + vi]++;
						else eqX[1] -= CurrSol[vi];
					}
					else{
						if (vi < eqvarc) eqX[2 + vi]--;
						else eqX[1] += CurrSol[vi];
					}
				}
				next++;
			}
		}
		int v = (table[((1 << a_first) - 1) & (_b ^ (_b >> (sol2.nglen+2)))]);
		if (i < eqvarc) eqX[2 + i] += v;
		else eqX[1] -= v * CurrSol[i];
	}
	for (int i = 0; i < sol2.varc; i++){
		Ngram *ngram = sol2.ng + i;
		ngram->eq2 = -1;
	}
	for (int i = 0; i < eqvarc; i++){
		if (eqX[i + 2] != 0) eqX[0]++;
	}
#ifdef OUTPUT
	sol->print();
	print_eqs(eqs, eqc);
	print_eq(eqX, 0);
	print_ng(sol2.ng, sol2.varc);
#endif
	normalize_eqX(eqX, eqvarc, eqs, eqc);
	if (eqX[1] < 0) goto end;
	int solc = 0;
	Solution = (sample*)Stack;
	SolutionEnd = (sample*)StackEnd;
	sample *sol2sol = Solution;
	solve_eq(eqX, eqs, eqc, solc);
	Stack = Stack + (solc * sol2.varc * sizeof(sample) + sizeof(long long) - 1) / sizeof(long long);
	if (solc > 0) BestLen = max(BestLen, sol2.nglen);
#ifdef OUTPUT
	int outp[2];
	outp[0] = solc;
	outp[1] = sol2.nglen;
	print_ints(outp, 2);
#endif
	if (sol2.nglen < sol2.maxlen){
		for (int i = 0; i < solc; i++){
			sol2.sol = sol2sol + i * sol2.varc;
			count += next_step(&sol2);
		}
	}
	if (sol2.nglen == OutLen){
		for (int i = 0; i < solc; i++){
			sol2.sol = sol2sol + i * sol2.varc;
			get_codes(&sol2);
		}
	}
end:
	for (int i = 0; i < sol->varc; i++){
		Ngram *ngram = sol->ng + i;
		ngram->eq1 = ngram->eq2 = -1;
		ngram->ng2 = NULL;
	}
	if (Stack >= StackEnd) MaxStackExc = true;
	Stack = stack;
	return count;
}
void normalize_eqX(int *eqX, int eqvarc, int **eqs, int eqc){
	for (int i = 0; i < eqvarc; i++){
		if (eqX[2 + i] >= 0) continue;
		int **zz = VarEqs + (i << 1);
		if ((*zz)[1] > (*(zz + 1))[1]) zz++;
		int *eq = *zz;
		int add = -eqX[2 + i];
		for (int j = 0; j < eq[0]; j++){
			int ind = 2 + eq[2 + j];
			if (eqX[ind] == 0) eqX[0]++;
			eqX[ind] += add;
			if (eqX[ind] == 0) eqX[0]--;
		}
		eqX[1] += eq[1] * add;
	}
	if (eqX[0] == 0 || eqX[1] <= 0) return;
start:
	int opti = 0, optv = 0;
	int **_eqs = eqs;
	for (int i = 0; i < eqc; i++, _eqs++){
		int *eq = *_eqs;
		int _optv = eq[1];
		for (int j = 0; j < eq[0]; j++){
			int vi = eq[2 + j];
			if (eqX[2 + vi] == 0){
				int **zz = VarEqs + (vi << 1);
				if ((*zz) == eq) zz++;
				_optv -= (*zz)[1];
			}
		}
		if (_optv > optv){ optv = _optv; opti = i; }
	}
	if (optv == 0) return;
	int *eq = eqs[opti];
	eqX[1] -= eq[1];
	for (int j = 0; j < eq[0]; j++){
		int vi = eq[2 + j];
		if (eqX[2 + vi] == 0){
			int **zz = VarEqs + (vi << 1);
			if ((*zz) == eq) zz++;
			int *_eq = *zz;
			eqX[1] += _eq[1];
			for (int jj = 0; jj < _eq[0]; jj++){
				if (_eq[2 + jj] == vi) continue;
				int ind = 2 + _eq[2 + jj];
				if (eqX[ind] == 0) eqX[0]++;
				eqX[ind]++;
			}
		}
		else{
			if (--(eqX[2 + vi]) == 0) eqX[0]--;
		}
	}
	if (eqX[0] == 0 || eqX[1] <= 0) return;
	goto start;
}
bool solve_eq(int *eqX, int **eqs, int eqc, int &solc){
	static int table[4] = {1, 1, 2, 6};
	if (eqc <= 0){
		if (eqX[1] > VarX) return true;
		sample *s = Solution + Varc;
		if (s > SolutionEnd){ MaxStackExc = true; return false; }
		memcpy(Solution, CurrSol, Varc * sizeof(sample));
		Solution = s;
		solc++;
		return true;
	}
	int opti = 0, optv = INT_MAX;
	int **_eqs = eqs;
	for (int i = 0; i < eqc; i++, _eqs++){
		int *eq = *_eqs;
		int c = eq[0] - 1, r = eq[1] + c;
		int _optv = r;
		for (int j = 0; j < c; j++){ r--; _optv *= r; }
		_optv /= table[c];
		if (_optv < optv){ optv = _optv; opti = i; if (optv == 0) break; }
	}
	eqc--;
	int *eq = eqs[opti];
	eqs[opti] = eqs[eqc];
	eqs[eqc] = eq;
	int c = *eq++;
	int r = *eq++;
	int *e[4];
	for (int i = 0; i < c; i++){
		int **zz = VarEqs + (eq[i] << 1);
		if (*zz == eq - 2) zz++;
		e[i] = *zz;
		for (int j = 0; j < e[i][0]; j++){
			if (e[i][j+2] == eq[i]){
				int tmp = e[i][j+2];
				e[i][j+2] = e[i][e[i][0]+1];
				e[i][e[i][0]+1] = tmp;
				break;
			}
		}
	}
	if (r == 0){
		for (int i = 0; i < c; i++){ e[i][0]--; CurrSol[eq[i]] = 0; }
		if (!solve_eq(eqX, eqs, eqc, solc)) return false;
		for (int i = 0; i < c; i++) e[i][0]++;
	}
	else if (c == 0){
		if (eqs[opti][1] != 0) return true;
		if (!solve_eq(eqX, eqs, eqc, solc)) return false;
	}
	else if (c == 1){
		if (r > e[0][1]) return true;
		e[0][0]--; e[0][1] -= r;
		CurrSol[eq[0]] = r;
		int xadd = r * eqX[2 + eq[0]];
		eqX[1] -= xadd;
		if (eqX[1] >= 0){
			if (!solve_eq(eqX, eqs, eqc, solc)) return false;
		}
		eqX[1] += xadd;
		e[0][0]++; e[0][1] += r;
	}
	else if (c == 2){
		for (int i0 = max(0, r - e[1][1]), i0m = min(r, e[0][1]); i0 <= i0m; i0++){
			e[0][0]--; e[0][1] -= i0;
			CurrSol[eq[0]] = i0;
			int i1 = r - i0;
			e[1][0]--; e[1][1] -= i1;
			CurrSol[eq[1]] = i1;
			int xadd = i0 * eqX[2 + eq[0]] + i1 * eqX[2 + eq[1]];
			eqX[1] -= xadd;
			if (eqX[1] >= 0){
				if (!solve_eq(eqX, eqs, eqc, solc)) return false;
			}
			eqX[1] += xadd;
			e[1][0]++; e[1][1] += i1;
			e[0][0]++; e[0][1] += i0;
		}
	}
	else if (c == 3){
		for (int i0 = 0, i0m = min(r, e[0][1]); i0 <= i0m; i0++){
			e[0][0]--; e[0][1] -= i0;
			CurrSol[eq[0]] = i0;
			int j1 = r - i0;
			for (int i1 = max(0, j1 - e[2][1]), i1m = min(j1, e[1][1]); i1 <= i1m; i1++){
				e[1][0]--; e[1][1] -= i1;
				CurrSol[eq[1]] = i1;
				int i2 = j1 - i1;
				e[2][0]--; e[2][1] -= i2;
				CurrSol[eq[2]] = i2;
				int xadd = i0 * eqX[2 + eq[0]] + i1 * eqX[2 + eq[1]] + i2 * eqX[2 + eq[2]];
				eqX[1] -= xadd;
				if (eqX[1] >= 0){
					if (!solve_eq(eqX, eqs, eqc, solc)) return false;
				}
				eqX[1] += xadd;
				e[2][0]++; e[2][1] += i2;
				e[1][0]++; e[1][1] += i1;
			}
			e[0][0]++; e[0][1] += i0;
		}
	}
	else{
		for (int i0 = 0, i0m = min(r, e[0][1]); i0 <= i0m; i0++){
			e[0][0]--; e[0][1] -= i0;
			CurrSol[eq[0]] = i0;
			for (int i1 = 0, i1m = min(r - i0, e[1][1]); i1 <= i1m; i1++){
				e[1][0]--; e[1][1] -= i1;
				CurrSol[eq[1]] = i1;
				int j2 = r - i0 - i1;
				for (int i2 = max(0, j2 - e[3][1]), i2m = min(j2, e[2][1]); i2 <= i2m; i2++){
					e[2][0]--; e[2][1] -= i2;
					CurrSol[eq[2]] = i2;
					int i3 = j2 - i2;
					e[3][0]--; e[3][1] -= i3;
					CurrSol[eq[3]] = i3;
					int xadd = i0 * eqX[2 + eq[0]] + i1 * eqX[2 + eq[1]] + i2 * eqX[2 + eq[2]] + i3 * eqX[2 + eq[3]];
					eqX[1] -= xadd;
					if (eqX[1] >= 0){
						if (!solve_eq(eqX, eqs, eqc, solc)) return false;
					}
					eqX[1] += xadd;
					e[3][0]++; e[3][1] += i3;
					e[2][0]++; e[2][1] += i2;
				}
				e[1][0]++; e[1][1] += i1;
			}
			e[0][0]++; e[0][1] += i0;
		}
	}
	return true;
}
int get_codes(Sol *sol){
	int i;
	int count = 0;
	for (i = 0; i < sol->varc; i++){
		if (sol->sol[i] != 0){
			Ngram *ngram = sol->ng + i;
			sol->sol[i]--;
			count = _get_codes(sol, ngram, 0, 0, 0);
			sol->sol[i]++;
			break;
		}
	}
	write(0, 0);
	return count;
}
int _get_codes(Sol *sol, Ngram *ngram, int pos, int mnm, int mxm){
	char a = (ngram->a & 3) + 1;
	CurrCode[pos++] = '0' + a;
	if (pos & 1) a = -a;
	mnm = min(a, mnm+a);
	mxm = max(a, mxm+a);
	if (mnm < -FILTER || mxm > FILTER) return 0;
	if (pos >= 2 * sol->k){
		int sum = 0;
		for (int i = 0; i < pos; i++){
			int _a = CurrCode[i] - '0';
			if (i & 1) _a = -_a;
			sum += _a;
			mnm = min(_a, mnm+_a);
			mxm = max(_a, mxm+_a);
			if (mnm < -FILTER || mxm > FILTER) return 0;
		}
		if (sum != 0) return 0;
		write(CurrCode, pos);
		return 1;
	}
	int count = 0;
	Ngram **next = ngram->next;
	while (*next){
		int i = (int)((*next) - sol->ng);
		if (sol->sol[i] > 0){
			sol->sol[i]--;
			count += _get_codes(sol, *next, pos, mnm, mxm);
			sol->sol[i]++;
		}
		next++;
	}
	return count;
}
void Ngram::print(){
	char buf[44];
	char *dst = StrA(buf);
	write(buf, (int)(dst - buf));
}
char* Ngram::StrA(char *buf){
	int c = count >> 1;
	int _a = a;
	for (int i = 0; i < c; i++, buf++){
		*buf = '1' + (_a & 3);
		_a >>= 2;
	}
	return buf;
}
void Sol::print(){
	static char buf[MAX_VARC * 12];
	char *dst = buf;
	for (int i = 0; i < varc; i++){
		if (sol[i] == 0) continue;
		if (dst != buf) *dst++ = ',';
		*dst++ = 'X';
		dst = ng[i].StrA(dst);
		*dst++ = '=';
		_itoa(sol[i], dst, 10);
		while (*dst) dst++;
	}
	*dst++ = '\n';
	write(buf, (int)(dst - buf));
}
void msg(int x){
	static char buf[20];
	_itoa(x, buf, 10);
	msg(buf);
}
void msg(char *x){
	MessageBoxA(NULL, x, "Dioph", MB_OK);
}
void print_ints(int *data, int count){
	int chs = 12;
	char empty = ' ';
	int buflen = count * chs;
	char *buf = new char[buflen];
	memset(buf, empty, buflen);
	for (int i = 0; i < count; i++, data++){
		char *dst = _itoa(*data, buf + i * chs, 10);
		dst[strlen(dst)] = empty;
	}
	write(buf, buflen);
}
void print_ng(Ngram *ng, int count){
	for (int i = 0; i < count; i++){
		ng[i].print();
	}
}
void print_eqs(int **eqs, int eqc){
	for (int i = 0; i < eqc; i++){
		print_eq(eqs[i], 1);
	}
}
void printline(){
	write("\r\n", 2);
}
void print_eq(int *eq, int mode){
	int varc = *eq++;
	int right = *eq++;
	int chs = 7;
	char empty = ' ';
	int buflen = (varc + 1) * chs;
	char *buf = new char[buflen];
	memset(buf, empty, buflen);
	char *dst = buf;
	if (mode == 0){
		int *_eq = eq;
		for (int i = 0; i < varc; eq++){
			if (!*eq) continue;
			_itoa(*eq, dst, 10);
			int l = strlen(dst);
			*(dst + l) = 'X';
			_itoa((int)(eq - _eq), dst + l + 1, 10);
			dst[strlen(dst)] = empty;
			i++;
			dst += chs;
		}
	}
	else{
		for (int i = 0; i < varc; i++, eq++, dst += chs){
			*dst = '+';
			*(dst + 1) = 'X';
			_itoa(*eq, dst + 2, 10);
			dst[strlen(dst)] = empty;
		}
	}
	*dst = '=';
	_itoa(right, dst + 1, 10);
	write(buf, buflen);
}
void write(char *s, int len){
	static FILE *f = NULL;
	if (s){
		if (!f) f = fopen(LOG_PATH, "ab");
		fwrite(s, len, 1, f);
		fwrite("\r\n", 2, 1, f);
	}
	else if (f) {
		if (len >= 0) fflush(f);
		else{ fclose(f); f = NULL; }
	}
}