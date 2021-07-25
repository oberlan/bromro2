#ifndef ___UTIL_H__
#define ___UTIL_H__

#include <vector>
#include <iostream>
#include <cstring>
#include <string>
#include <iomanip>
#include <time.h>
#include <utility>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits.h>
#include <sys/time.h>

using namespace std;

/************************************************************************/
/* DEBUG                                                                */
/************************************************************************/

//#define IMPRIMIR_INFO
//#define IMPRIMIR_INSTANCIA
#define IMPRIMIR_CPLEX_LOG
#define IMPRIMIR_SOLUCAO
#define PASSEIAQUI(i) cerr << "ok" << i << endl;

/************************************************************************/
/** Variáveis da PDA                                                    */
/************************************************************************/
#define NUMERO_MAX_ITERACAO 500
//#define NUMERO_MAX_ITERACAO_SEM_MELHORA 50
static int NUM_MAX_ITERACAO_SEM_MELHORA = 50;
static double valorInicialDelta = 0.90;
#define PODE_ZERAR_DELTAS
#define AREANULA_DELTAZERO
static double valorSigma = 0.90;
static int iteracao = 1;

#define SIGMA (pow(valorSigma, iteracao))

typedef double tipoCusto;


/************************************************************************/
/** FUNCOES AUXILIARES                                                  */
/************************************************************************/
#ifdef _DEBUG
#define vetor(TIPO) scpp::vector<TIPO>
#define GET(VETOR, POS) VETOR.at(POS)
#else
#define vetor(TIPO) std::vector<TIPO>
#define GET(VETOR, POS) VETOR[POS]
#endif

#define INF INT_MAX

#define ZERO 0.0000001
#define DOUBLE_EQUALS(X, Y) (fabs((X)-(Y)) < ZERO)
#define COMPARA_DELTAS(DELTA_NOVO, DELTA_VELHO) ( ( fabs( (DELTA_NOVO)-(DELTA_VELHO) ) )  <= 0.01 * (DELTA_VELHO) )
//#define COMPARA_DELTAS(X, Y) ((fabs((X)-(Y))) < 0.01)

//#define DOUBLE_TO_INT(X) ( (int) ((X) + ZERO) )
#define DOUBLE_TO_INT(X) ( (int) ((X) + 0.0001) )
#define MIN(X,Y) ((X)<(Y) ? (X) : (Y))
#define MAX(X,Y) ((X)>(Y) ? (X) : (Y))

#define IGUAL(X, Y) fabs((X)-(Y))<ZERO

//#define MOD(X, Y) ((X)-((X)/(Y))*(Y))


//#ifdef __APPLE__
#ifdef _WIN32
#define PAUSE system("pause");
#else
#define PAUSE printf("Pressione uma tecla para continuar...\n"); getchar();
#endif

#define ASSERT(x) \
                 if (! (x)) \
                { \
                  cerr << "ERROR!! Assert " << #x << " failed\n"; \
                  cerr << " on line " << __LINE__  << "\n"; \
                  cerr << " in file " << __FILE__ << "\n";  \
				  exit(0);	\
                }

 /************************************************************************/
/* CLASSE PARA IMPRIMIR VIRGULA EM NUMEROS DECIMAIS                      */
/************************************************************************/
template<typename CharT>
class DecimalSeparator : public std::numpunct<CharT> {
public:
	DecimalSeparator(CharT Separator)
		: m_Separator(Separator)
	{}

protected:
	CharT do_decimal_point()const {
		return m_Separator;
	}

private:
	CharT m_Separator;
};

template<class InIt>
void print_range(InIt first, InIt last, char const* delim = "\n"){
    --last;
    for(; first != last; ++first){
        std::cout << *first << delim;
    }
    std::cout << *first;
};

class MyError : public logic_error {
public:
  string errorMsg;
  explicit MyError(const string& msg = "") : logic_error(msg), errorMsg(msg) {}
  virtual const char * what() const throw() {
    return errorMsg.c_str();
  }
  ~MyError() throw() {}
};


#if defined (MSVC)
#define INLINE __inline
#else
#define INLINE inline
#endif

INLINE long double MyClock(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + tv.tv_usec;
};

// static double __tempoInicio = MyClock();
//static bool usarValoresDuais = false;

// #define ACABOU_TEMPO ((MyClock() - __tempoInicio)/1000000.0) >= tempoMaximo

struct TuplaCor{
    string nome, tipo, valor;
    TuplaCor(string n, string t, string v):nome(n), tipo(t), valor(v){}
};

ostream& operator<<(ostream &out, const TuplaCor &t){
    out << t.nome;
    return out;
}

struct Par{
    int periodo;
    int objeto;
    Par(int p=0, int i=0):periodo(p),objeto(i){}
};

ostream& operator<<(ostream &out, const Par &p){
    out << "(" << p.periodo << ", " << p.objeto << ")";
    return out;
}


#endif
