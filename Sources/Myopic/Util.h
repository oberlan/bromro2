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
#include "scpp_assert.h"
#include "scpp_vector.h"

using namespace std;

typedef unsigned long ulong;
typedef unsigned long long ullong;
typedef long long llong;
typedef double tipoCusto;

/************************************************************************/
/** DEBUG                                                               */
/************************************************************************/

//#define IMPRIMIR_INFO
//#define IMPRIMIR_SAIDA_SIMULACAO
//#define IMPRIMIR_INSTANCIA
#define IMPRIMIR_CPLEX_LOG
#define IMPRIMIR_SOLUCAO

/************************************************************************/
/** Variáveis da PDA                                                    */
/************************************************************************/
#define NUMERO_MAX_ITERACAO 500
//#define NUMERO_MAX_ITERACAO_SEM_MELHORA 50
static int NUM_MAX_ITERACAO_SEM_MELHORA = 50;
static double valorInicialDelta = 0.8;
#define PODE_ZERAR_DELTAS
#define AREANULA_DELTAZERO
static double valorSigma = 0.50;//0.90
static ulong iteracao = 1;

#define LIMITAR_CUSTO_OBJETOS
#define CTE_MULT_CUSTOOBJETO 2


#define ALPHA0 (0.5)
#define ALPHA1 (pow(valorSigma, iteracao))
#define ALPHA2 (pow(1.0/iteracao, 0.8))
#define ALPHA3 ((double)iteracao/(1.0+iteracao-0.0))


#define ALPHA ALPHA1


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
#define DOUBLE_TO_INT(X) ( (ulong) ((X) + 0.0001) )
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
/* CLASSE PARA IMPRIMIR VIRGULA EM NUMEROS DECIMAIS                     */
/************************************************************************/
template<typename CharT>
class DecimalSeparator : public std::numpunct<CharT> {
public:
	DecimalSeparator(CharT Separator)
		: m_Separator(Separator) {}
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

INLINE double MyClock() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec * 1000000) + tv.tv_usec;
};


struct TuplaCor{
    string nome, tipo, valor;
    TuplaCor(string n, string t, string v):nome(n), tipo(t), valor(v){}
};

struct Par{
    ulong periodo;
    ulong objeto;
    Par(ulong p=0, ulong i=0):periodo(p),objeto(i){}
};

static ostream& operator<<(ostream &out, const Par &p){
    out << "(" << p.periodo << ", " << p.objeto << ")";
    return out;
}

struct Distribuicao{
    double media;
    double desvioPadrao;
    Distribuicao(double m=0.0, double sd=0.0):media(m),desvioPadrao(sd){}
};



#endif
