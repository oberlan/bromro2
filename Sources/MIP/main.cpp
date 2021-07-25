
#include <iostream>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <list>
#include <set>
#include <limits>
#include <stdexcept>
#include <string>
#include <cmath>
#include <sstream>
#include <sys/time.h>

#include "Util.h"

#define CPLEX_GET(X) DOUBLE_TO_INT(cplex.getValue(X))
#define CPLEX_GETD(X) cplex.getValue(X)
#define THREADS 2
#define EXPORTAR_MODELO
#define EXECUTABLE_PATH "./"
#define DIR "../../Instances/bromro/"
#define INSTANCIA "2.txt"

using namespace std;

typedef std::numeric_limits< double > dbl;

typedef IloArray<IloNumVarArray> NumVarMatrix2D;
typedef IloArray<NumVarMatrix2D> NumVarMatrix3D;
typedef IloArray<IloNumArray>    FloatMatrix;
typedef IloArray<IloExprArray>   ExprMatrix;
typedef unsigned long long ullong;

struct Retangulo {
    int w, h, c, p;
    double x, y;
    Retangulo(int _w=0, int _h=0, int _c=0, int _x=0, int _y=0, int _p=0):w(_w),h(_h),c(_c),x(_x),y(_y),p(_p)
    {}
    bool operator<(const Retangulo& r) const {
        if(w < r.w) return true;
        else if(w == r.w) return (h < r.h);
        else return false;
    }
    int getCusto() const {
        return w*h*c;
    }
    void print() const {
        cout << w << " x " << h << "\t(" << x << ", " << y << ")\n";
    }
    int getArea() const{
        return w*h;
    }
};

struct Ponto{
    double x, y;
    Ponto(double _x=0, double _y=0):x(_x),y(_y){}
    bool operator==(const Ponto& pt) const{
        return (DOUBLE_EQUALS(x, pt.x) && DOUBLE_EQUALS(y, pt.y));
    }
};

void leEntrada(const string &arquivoEntrada);
void resolveModelo();
void setParametros(IloCplex &cplex);


bool is_file_exist(const char *fileName);
void adicionaPontoDeCorte(vetor(Ponto)&, Ponto&);
void leArgumentos(int argc, char** argv);
bool cabeItemCatalogo(vetor(int) &w_bar, vetor(int) &h_bar, Retangulo& sobra);

//Variáveis (dados) de entrada
int P, d, xi;
vetor(int) n;
vetor(int) m;
vetor(int) m_bar;
//vetor(int) c_bar;
vetor(vetor(int)) c_bar;
vetor(int) p1;
vetor(int) w_bar;
vetor(int) h_bar;
vetor( vetor(int) ) n_q;
vetor( vetor(int) ) o_q;

vetor( vetor(Retangulo) ) item;
vetor( vetor(Retangulo) ) objeto;

vetor( vetor(Retangulo) ) objetosSobras;
vetor( vetor( vetor(Retangulo) ) ) listaItensDoObjeto;
vetor( vetor( vetor(Ponto) ) ) posicaoCorteDoObjeto;
vetor( Retangulo ) sobrasFinais;
vetor( Par ) refSobrasFinais;
vector< vector<bool> > objetoComItem;


bool adicionarRestricaoArea = true;

//Indica de qual objeto e período a sobra foi originada
//Objetos que podem ser comprados (1...m[s]) referenciam a sí próprio
vetor( vetor(Par) ) refObjPer;

int WMax = 0;
int HMax = 0;

string nomeInstancia;
string nomeArquivoInstancia;
int tempoLimiteCplex;
bool imprimirLogCplex;
bool setLimiteThreadsCplex;
int numThreadsCplex;
bool setLimiteMemoriaCplex;
int maxMemoriaCplex;
int mipStrategyFile;

vector< TuplaCor > cores;
vector< pair<string, string> > padrao;

int main(int argc, char** argv) {
	std::cout.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));
    auto tempoInicio = MyClock();

    leArgumentos(argc, argv);

	leEntrada(nomeInstancia);

    cout << "\nItens do catalogo: " << d << endl;
    for(int i = 1; i <= d; i++){
        cout << w_bar[i] << " x " << h_bar[i] << endl;
    }
    cout << endl;

	resolveModelo();

    cout << "Tempo Total Gasto: " << setprecision(5) << fixed << (MyClock()-tempoInicio)/CLOCKS_PER_SEC << endl << endl;

	return 0;
}

void leArgumentos(int argc, char** argv) {
    tempoLimiteCplex = 7200;
    imprimirLogCplex = false;
    setLimiteMemoriaCplex = true;
    maxMemoriaCplex = 16000;
    setLimiteThreadsCplex = false;
    adicionarRestricaoArea = true;
    xi = 4;
    mipStrategyFile = 1;

    if(argc == 1) {
        nomeInstancia = string(DIR) + string(INSTANCIA);
        tempoLimiteCplex = 7200;
        imprimirLogCplex = true;
        setLimiteMemoriaCplex = true;
        maxMemoriaCplex = 16000;
        setLimiteThreadsCplex = true;
        numThreadsCplex = 1;
        adicionarRestricaoArea = true;
        xi = 4;
    }

    for(int i = 1 ; i < argc ; i++) {
        if (string(argv[i]) == "-i") {
            nomeInstancia = string(argv[i + 1]);
            continue;
        }
        if (string(argv[i]) == "-tl") {
            tempoLimiteCplex = atoi(argv[i+1]);
            continue;
        }
        if (string(argv[i]) == "-xi") {
            xi = atoi(argv[i+1]);
            continue;
        }
        if (string(argv[i]) == "-logCplex") {
            imprimirLogCplex = true;
            continue;
        }
        if(string(argv[i]) == "-mem"){
            setLimiteMemoriaCplex = true;
            maxMemoriaCplex = atoi(argv[i+1]);
            continue;
        }
        if(string(argv[i]) == "-threads"){
            setLimiteThreadsCplex = true;
            numThreadsCplex = atoi(argv[i+1]);
            continue;
        }
        if(string(argv[i]) == "-nodefileind") {
            mipStrategyFile = atoi(argv[i + 1]);
        }
    }

    //separando o nome da instancia
    string::size_type loc  = nomeInstancia.find_last_of( "/", nomeInstancia.size() );
    string::size_type loc2 = nomeInstancia.find_last_of( ".", nomeInstancia.size() );
    if( !((loc != string::npos ) && (loc2 != string::npos)) ){
        cout << "ERRO!" << endl;
        exit(0);
    }
    nomeArquivoInstancia = string("");
    nomeArquivoInstancia.append(nomeInstancia, loc+1, loc2-loc-1);

    cout << "\n\nInstancia: " << nomeInstancia << endl;
    cout << "xi: " << xi << endl;
    cout << "Tempo Limite Cplex: " << tempoLimiteCplex << endl;
    cout << "Imprimir Log Cplex? " << std::boolalpha << imprimirLogCplex << endl;
}

void resolveModelo() {
    auto inicio = MyClock();
	IloEnv env;
	cout << env.getVersion() << endl;
	int i, j, k, l, p, s1, j1, j2;
	char str[30];

	try {
		IloModel modelo(env, "Cutting");
		IloCplex cplex(env);

		//Variáveis auxiliares
        auto L = (unsigned int) (floor(log2(WMax)) + 1);
		vector<int> pot(L+1);
        for(k=1; k<=L; ++k)
            pot[k] = (1<<(k-1));

		//Variáveis de decisão
		NumVarMatrix2D u(env, P+1);
		NumVarMatrix3D v(env, P+1);
		NumVarMatrix2D x(env, P+1);
		NumVarMatrix2D y(env, P+1);
		NumVarMatrix3D pi(env, P+1);
		NumVarMatrix3D tau(env, P+1);
		NumVarMatrix2D r(env, P+1);
		NumVarMatrix2D t(env, P+1);
		NumVarMatrix2D W_bar(env, P+1);
		NumVarMatrix2D H_bar(env, P+1);
		IloNumVarArray gamma(env, m_bar[P]+1, 0.0, IloInfinity, ILOFLOAT);
		NumVarMatrix2D eta(env, P+1);
        NumVarMatrix2D theta(env, m_bar[P]+1);
        NumVarMatrix2D mu(env, m_bar[P]+1);
        NumVarMatrix2D q(env, m_bar[P]+1);
        //IloNumVarArray omega(env, m_bar[P]+1, 0.0, IloInfinity, ILOFLOAT);
        NumVarMatrix2D omega(env, m_bar[P]+1);
        NumVarMatrix2D zeta(env, m_bar[P]+1);

        for(p = 0; p <= P; ++p) {
            u[p] = IloNumVarArray(env, m_bar[p]+1, 0, 1, ILOBOOL);
            v[p] = NumVarMatrix2D(env, n[p]+1);
            pi[p] = NumVarMatrix2D(env, n[p]+1);
            tau[p] = NumVarMatrix2D(env, n[p]+1);
            r[p] = IloNumVarArray(env, m_bar[p]+1, 0.0, IloInfinity, ILOFLOAT);
            t[p] = IloNumVarArray(env, m_bar[p]+1, 0.0, IloInfinity, ILOFLOAT);
            x[p] = IloNumVarArray(env, n[p]+1, -IloInfinity, IloInfinity, ILOFLOAT);
            y[p] = IloNumVarArray(env, n[p]+1, -IloInfinity, IloInfinity, ILOFLOAT);
            W_bar[p] = IloNumVarArray(env, m_bar[p]+1, 0.0, IloInfinity, ILOFLOAT);
            H_bar[p] = IloNumVarArray(env, m_bar[p]+1, 0.0, IloInfinity, ILOFLOAT);
            eta[p] = IloNumVarArray(env, m_bar[p]+1, 0, 1, ILOBOOL);
            for(i = 1; i <= n[p]; ++i){
                  v[p][i] = IloNumVarArray(env, m_bar[p]+1, 0, 1, ILOBOOL);
                 pi[p][i] = IloNumVarArray(env, n[p]+1, 0, 1, ILOBOOL);
                tau[p][i] = IloNumVarArray(env, n[p]+1, 0, 1, ILOBOOL);
            }
        }

        for(j = 1; j <= m_bar[P]; ++j){
            theta[j] = IloNumVarArray(env, L+1, 0, 1, ILOBOOL);
            omega[j] = IloNumVarArray(env, L+1, 0.0, IloInfinity, ILOFLOAT);
             zeta[j] = IloNumVarArray(env, d+1, 0, 1, ILOBOOL);
               mu[j] = IloNumVarArray(env, d+1, -IloInfinity, IloInfinity, ILOFLOAT);
                q[j] = IloNumVarArray(env, d+1, 0, 1, ILOBOOL);
        }

        for(p = 0; p <= P; ++p) {
            char nome[30];
            for(j = 1; j <=  m_bar[p]; ++j){
                sprintf(nome, "u[%d][%d]", p, j);
                u[p][j].setName(nome);

                sprintf(nome, "H_bar[%d][%d]", p, j);
                H_bar[p][j].setName(nome);
                sprintf(nome, "W_bar[%d][%d]", p, j);
                W_bar[p][j].setName(nome);

                //if(p != P){
                    sprintf(nome, "r[%d][%d]", p, j);
                    r[p][j].setName(nome);
                    sprintf(nome, "t[%d][%d]", p, j);
                    t[p][j].setName(nome);
                //}
            }
            for(i = 1; i <= n[p]; ++i){
                sprintf(nome, "x[%d][%d]", p, i);
                x[p][i].setName(nome);
                sprintf(nome, "y[%d][%d]", p, i);
                y[p][i].setName(nome);
                for(int il = i + 1; il <= n[p]; ++il) {
                    sprintf(nome, "pi[%d][%d][%d]", p, i, il);
                    pi[p][i][il].setName(nome);
                    sprintf(nome, "tau[%d][%d][%d]", p, i, il);
                    tau[p][i][il].setName(nome);
                }
                for(j=1; j<= m_bar[p]; ++j){
                    sprintf(nome, "v[%d][%d][%d]", p, i, j);
                    v[p][i][j].setName(nome);
                }
            }
        }
        for(j = 1; j <= m_bar[P]; ++j){
            char nome[30];
            sprintf(nome, "gamma[%d]", j);
            gamma[j].setName(nome);
            for(l = 1; l <= L; ++l) {
                sprintf(nome, "theta[%d][%d]", j, l);
                theta[j][l].setName(nome);

                sprintf(nome, "omega[%d][%d]", j,l);
                omega[j][l].setName(nome);
            }
            for(l = 1; l <= d; ++l){
                sprintf(nome, "zeta[%d][%d]", j, l);
                zeta[j][l].setName(nome);

                sprintf(nome, "mu[%d][%d]", j, l);
                mu[j][l].setName(nome);

                sprintf(nome, "q[%d][%d]", j, l);
                q[j][l].setName(nome);
            }
        }


		IloExpr objetosUsados(env);
		IloExpr sobrasUsadas(env);
		IloExpr objectCust(env);
		IloExpr leftoversValue(env);
		IloExpr funcaoObjetivo(env);
		ullong PESO = 0;
		for(p = 0; p < P; ++p){
            for(j = 1; j <= m[p]; ++j){
                objectCust += (objeto[p][j].getCusto()*u[p][j]);
                objetosUsados += (u[p][j]);
                PESO += objeto[p][j].getCusto();
            }
            for(j = m[p] + 1; j <= m_bar[p]; ++j)
                sobrasUsadas += (u[p][j]);
        }
        for(j = 1; j <= m_bar[P]; ++j) {
            p = refObjPer[P][j].periodo;
            int obj = refObjPer[P][j].objeto;
            leftoversValue += (objeto[p][obj].c * gamma[j]);
        }


		funcaoObjetivo = (PESO*objectCust - leftoversValue);

		//=====================================================Funcao Objetivo
		modelo.add(IloMinimize(env, funcaoObjetivo));

		//=====================================================Restricoes
		IloRangeArray restricoes(env);


        //=====================================================
        for(p = 0; p < P; ++p)
            for(i = 1; i <= n[p]; ++i){
                IloExpr somatorio(env);
                for(j = 1; j <= m_bar[p]; ++j)
                    somatorio += v[p][i][j];
                restricoes.add(somatorio == 1);
                somatorio.end();
            }

		//=====================================================
        for(p = 0; p < P; ++p)
            for(i = 1; i <= n[p]; ++i)
                for(j = 1; j <= m_bar[p]; ++j)
                    restricoes.add(u[p][j] - v[p][i][j] >= 0);

        //****************************************************************************************//
        //Objetos só podem ser selecionados se tiverem algum item atribuido  //
        for(p = 0; p < P; ++p)
            for(j = 1; j <= m_bar[p]; ++j){
                IloExpr somatorio(env);
                for(i = 1; i <= n[p]; ++i)
                    somatorio += v[p][i][j];
                restricoes.add(u[p][j] - somatorio <= 0);
                somatorio.end();
            }
        //****************************************************************************************//


        //Atribui o tamanho de cada objeto
        for(p = 0; p < P; ++p)
            for(j = 1; j <= m[p]; ++j){
                restricoes.add(W_bar[p][j] - objeto[p][j].w == 0);
                restricoes.add(H_bar[p][j] - objeto[p][j].h == 0);
            }


        //=====================================================
        for(p = 0; p < P; ++p) {
            for (j = 1; j <= m_bar[p]; ++j) {
                restricoes.add(t[p][j] - H_bar[p][j] <= 0);
                restricoes.add(r[p][j] - W_bar[p][j] <= 0);
            }
        }


        //========================================== Simetria do eta
        for(p = 0; p < P; ++p)
            for(j = 1; j <= m_bar[p]; ++j){
                restricoes.add(eta[p][j] - u[p][j] <= 0);
            }

        //=====================================================
        for(p = 0; p < P; ++p)
            for(j = 1; j <= m[p]; ++j){
                j1 = m[p + 1] + 2 * j - 1;
                j2 = m[p + 1] + 2 * j;
                if(j1 <= m_bar[p+1]) {
                    restricoes.add(H_bar[p + 1][j1] - (u[p][j] * HMax) <= 0);
                    restricoes.add(t[p][j] - (1 - u[p][j]) * HMax - H_bar[p + 1][j1] <= 0);
                    restricoes.add(H_bar[p + 1][j1] - (t[p][j] + (1 - u[p][j]) * HMax) <= 0);

                    restricoes.add(W_bar[p + 1][j1] - (u[p][j] * WMax) <= 0);
                    restricoes.add(W_bar[p][j] - r[p][j] - (1 - eta[p][j]) * WMax - (1 - u[p][j]) * WMax - (W_bar[p + 1][j1]) <= 0);
                    restricoes.add(W_bar[p + 1][j1] - (W_bar[p][j] - r[p][j] + (1 - eta[p][j]) * WMax + (1 - u[p][j]) * WMax) <= 0);
                    restricoes.add(W_bar[p][j] - eta[p][j] * WMax - (1 - u[p][j]) * WMax - (W_bar[p + 1][j1]) <= 0);
                    restricoes.add(W_bar[p + 1][j1] - (W_bar[p][j] + eta[p][j] * WMax + (1 - u[p][j]) * WMax) <= 0);
                }
                if(j2 <= m_bar[p+1]) {
                    restricoes.add(W_bar[p + 1][j2] - (u[p][j] * WMax) <= 0);
                    restricoes.add(r[p][j] - (1 - u[p][j]) * WMax - (W_bar[p + 1][j2]) <= 0);
                    restricoes.add(W_bar[p + 1][j2] - (r[p][j] + (1 - u[p][j]) * WMax) <= 0);

                    restricoes.add(H_bar[p + 1][j2] - (u[p][j] * HMax) <= 0);
                    restricoes.add(H_bar[p][j] - (1 - eta[p][j]) * HMax - (1 - u[p][j]) * HMax - (H_bar[p + 1][j2]) <= 0);
                    restricoes.add(H_bar[p + 1][j2] - (H_bar[p][j] + (1 - eta[p][j]) * HMax + (1 - u[p][j]) * HMax) <= 0);
                    restricoes.add(H_bar[p][j] - t[p][j] - eta[p][j] * HMax - (1 - u[p][j]) * HMax - (H_bar[p + 1][j2]) <= 0);
                    restricoes.add(H_bar[p + 1][j2] - (H_bar[p][j] - t[p][j] + eta[p][j] * HMax + (1 - u[p][j]) * HMax) <= 0);
                }
            }

        //=====================================================
        for(p = 0; p < P; ++p)
            for(j = m[p] + 1; j <= m_bar[p]; ++j) {
                j1 = m[p + 1] + 2 * j - 1;
                j2 = m[p + 1] + 2 * j;
                if(j1 <= m_bar[p+1]) {
                    restricoes.add(H_bar[p][j] - u[p][j] * HMax - (H_bar[p + 1][j1]) <= 0);
                    restricoes.add(H_bar[p + 1][j1] - (H_bar[p][j] + u[p][j] * HMax) <= 0);
                    restricoes.add(t[p][j] - (1 - u[p][j]) * HMax - (H_bar[p + 1][j1]) <= 0);
                    restricoes.add(H_bar[p + 1][j1] - (t[p][j] + (1 - u[p][j]) * HMax) <= 0);

                    restricoes.add(W_bar[p][j] - u[p][j] * WMax - (W_bar[p + 1][j1]) <= 0);
                    restricoes.add(W_bar[p + 1][j1] - (W_bar[p][j] + u[p][j] * WMax) <= 0);
                    restricoes.add(W_bar[p][j] - r[p][j] - (1 - eta[p][j]) * WMax - (1 - u[p][j]) * WMax - (W_bar[p + 1][j1]) <= 0);
                    restricoes.add(W_bar[p + 1][j1] - (W_bar[p][j] - r[p][j] + (1 - eta[p][j]) * WMax + (1 - u[p][j]) * WMax) <= 0);
                    restricoes.add(W_bar[p][j] - eta[p][j] * WMax - (1 - u[p][j]) * WMax - (W_bar[p + 1][j1]) <= 0);
                    restricoes.add(W_bar[p + 1][j1] - (W_bar[p][j] + eta[p][j] * WMax + (1 - u[p][j]) * WMax) <= 0);
                }
                if(j2 <= m_bar[p+1]) {
                    restricoes.add(W_bar[p + 1][j2] - (u[p][j] * WMax) <= 0);
                    restricoes.add(r[p][j] - (1 - u[p][j]) * WMax - (W_bar[p + 1][j2]) <= 0);
                    restricoes.add(W_bar[p + 1][j2] - (r[p][j] + (1 - u[p][j]) * WMax) <= 0);

                    restricoes.add(H_bar[p + 1][j2] - (u[p][j] * HMax) <= 0);
                    restricoes.add(H_bar[p][j] - (1 - eta[p][j]) * HMax - (1 - u[p][j]) * HMax - (H_bar[p + 1][j2]) <= 0);
                    restricoes.add(H_bar[p + 1][j2] - (H_bar[p][j] + (1 - eta[p][j]) * HMax + (1 - u[p][j]) * HMax) <= 0);
                    restricoes.add(H_bar[p][j] - t[p][j] - eta[p][j] * HMax - (1 - u[p][j]) * HMax - (H_bar[p + 1][j2]) <= 0);
                    restricoes.add(H_bar[p + 1][j2] - (H_bar[p][j] - t[p][j] + eta[p][j] * HMax + (1 - u[p][j]) * HMax) <= 0);
                }
            }


		//=====================================================
        for(p = 0; p < P; ++p)
            for(i=1; i <= n[p]; ++i){
                restricoes.add(x[p][i] - item[p][i].w/2.0 >= 0);
                restricoes.add(y[p][i] - item[p][i].h/2.0 >= 0);
            }
		for(p = 0; p < P; ++p)
            for(i = 1; i <= n[p]; ++i)
                for(j = 1; j <= m_bar[p]; ++j){
                    restricoes.add(x[p][i] + item[p][i].w/2.0 - (W_bar[p][j] - r[p][j] + (WMax)*(1-v[p][i][j])) <= 0);
                    restricoes.add(y[p][i] + item[p][i].h/2.0 - (H_bar[p][j] - t[p][j] + (HMax)*(1-v[p][i][j])) <= 0);
                }

        //=====================================================

        for(p = 0; p < P; ++p)
            for(j = 1; j <= m_bar[p]; ++j)
                for(i = 1; i <= n[p]; ++i)
                    for(int il = i+1; il <= n[p]; ++il){
                        restricoes.add( x[p][i] - x[p][il]  - (0.5*(item[p][i].w + item[p][il].w) - WMax*((1-v[p][i][j]) + (1-v[p][il][j]) + pi[p][i][il] + tau[p][i][il])) >= 0);
                        restricoes.add(-x[p][i] + x[p][il]  - (0.5*(item[p][i].w + item[p][il].w) - WMax*((1-v[p][i][j]) + (1-v[p][il][j]) + pi[p][i][il] + (1 - tau[p][i][il]))) >= 0);

                        restricoes.add( y[p][i] - y[p][il]  - (0.5*(item[p][i].h + item[p][il].h) - HMax*((1-v[p][i][j]) + (1-v[p][il][j]) + (1 - pi[p][i][il]) + tau[p][i][il])) >= 0);
                        restricoes.add(-y[p][i] + y[p][il]  - (0.5*(item[p][i].h + item[p][il].h) - HMax*((1-v[p][i][j]) + (1-v[p][il][j]) + (1 - pi[p][i][il]) + (1 - tau[p][i][il]))) >= 0);
                    }

        for(j = 1; j <= m_bar[P]; ++j){
            IloExpr somatorio(env);
            for(l = 1; l <= L; ++l){
                somatorio += pow(2, l-1)*theta[j][l];
            }
            restricoes.add(W_bar[P][j] - somatorio == 0);
        }

		//=====================================================11
        for(j = 1; j <= m_bar[P]; ++j)
            for(l = 1; l <= L; ++l){
                restricoes.add(omega[j][l] - H_bar[P][j] <= 0);
                restricoes.add(H_bar[P][j] - (1 - theta[j][l])*HMax - omega[j][l] <= 0 );
                restricoes.add(omega[j][l] - (theta[j][l]*HMax) <= 0 );
//                restricoes.add(omega[j][l] - W_bar[P][j] <= 0);
//                restricoes.add(W_bar[P][j] - (1 - theta[j][l])*WMax - omega[j][l] <= 0 );
//                restricoes.add(omega[j][l] - (theta[j][l]*WMax) <= 0 );
            }

        //=====================================================
        for(j = 1; j <= m_bar[P]; ++j)
            for(i = 1; i <= d; ++i){
                restricoes.add(w_bar[i] - (W_bar[P][j] + WMax*(1 - zeta[j][i])) <= 0);
                restricoes.add(h_bar[i] - (H_bar[P][j] + HMax*(1 - zeta[j][i])) <= 0);
            }

        //=====================================================
        for(j = 1; j <= m_bar[P]; ++j){
            IloExpr somatorio(env);
            for(l = 1; l <= L; ++l){
                somatorio += pow(2, l-1)*omega[j][l];
            }
            restricoes.add(gamma[j] - somatorio <= 0);
        }
        for(j = 1; j <= m_bar[P]; ++j){
            IloExpr somatorio(env);
            for(i = 1; i <= d; ++i){
                somatorio += zeta[j][i];
            }
            restricoes.add(gamma[j] - somatorio*WMax*HMax <= 0);
        }

		modelo.add(restricoes);
		cplex.extract(modelo);

        //cplex.exportModel("modelo.mps");

		if(!imprimirLogCplex)
            cplex.setOut(env.getNullStream());


		setParametros(cplex);
        string MIPSolution = string(string(DIR) + "/MIPSolutions/" + nomeArquivoInstancia + ".mst");
        #ifdef CARREGAR_SOLUCAO_CPLEX
        if(is_file_exist(MIPSolution.c_str()))
            cplex.readMIPStarts(MIPSolution.c_str());
        #endif
        cout << " Numero de variaveis binarias: " << cplex.getNbinVars() << endl;
        cout << " Numero de variaveis inteiras: " << cplex.getNintVars() << endl;
        cout << "Numero de variaveis continuas: " << cplex.getNcols() - cplex.getNbinVars() - cplex.getNintVars() << endl;
        cout << "         Numero de restricoes: " << cplex.getNrows() << endl;


		cplex.solve();

        cout << "[###] Solver status: " << cplex.getStatus() << endl;
   		IloInt accumulatedIterations = cplex.getNiterations();
   		IloInt accumulatedNodes = cplex.getNnodes();
   		IloNum fObjBestBound = cplex.getBestObjValue();
   		IloNum fObjValue = cplex.getObjValue();
   		double cplexGap = (fObjValue - fObjBestBound)/(0.0000000001 + fObjValue);
        cout << fixed << setprecision(10);

   		//cout.precision(dbl::digits10);
  		cout << "[###] Best objective bound: " << fObjBestBound << endl;
  		cout << "[###] Cplex GAP: "            << cplexGap*100.0 << "%" << endl;
   		cout << "[###] Objective obtained: "   << fObjValue << endl;
        cout << "[###] MIP Iterations: "       << accumulatedIterations << endl;
        cout << "[###] B&B Nodes: "            << accumulatedNodes << endl << endl;

        cout << "----------------SOLUTION FOUND----------------\n";
        cout << "Objects cost: "    << CPLEX_GETD(objectCust) << endl;
        cout << "Leftovers value: " << CPLEX_GETD(leftoversValue) << endl;
        cout << "Objects used: "    << CPLEX_GETD(objetosUsados) << endl;
        cout << "Leftovers used: "  << CPLEX_GETD(sobrasUsadas) << endl;
        cout << "Cplex time: "      << (MyClock()-inicio)/CLOCKS_PER_SEC << endl << endl;

        cout << xi << "\t" << fObjBestBound << "\t" << fObjValue << "\t" << CPLEX_GET(objectCust) << "\t" << CPLEX_GET(leftoversValue)
             << "\t" << cplexGap*100 << "%\t" << accumulatedIterations << "\t" << accumulatedNodes << "\t"
             << (MyClock()-inicio)/CLOCKS_PER_SEC << endl << endl;

        cerr << fixed << setprecision(10) <<  nomeArquivoInstancia << "\t" << xi << "\t" << fObjBestBound << "\t" << fObjValue << "\t"
        << CPLEX_GET(objectCust) << "\t" << CPLEX_GET(leftoversValue) << "\t" << cplexGap*100
        << "%\t" << accumulatedIterations << "\t" << accumulatedNodes << "\t" << (MyClock()-inicio)/CLOCKS_PER_SEC << endl;

	}
	catch (IloException& ex){
		cerr << "Error: " << ex << endl;
		env.end();
		exit(0);
	}
	catch (...){
		cerr << "Error" << endl;
		env.end();
		exit(0);
	}

	env.end();

}

void print(int *v, int ini, int fim){
    for(int i = ini; i<=fim; ++i)
        cout << v[i] << "\t";
    cout << endl;
}

void print(vetor(int) v){
    for(auto x: v)
        cout << x << "\t";
    cout << endl;
}

/*
* P
* n1 n2 ... nP
* w[1][1] w[1][2] ... w[1][n1]
* h[1][1] h[1][2] ... h[1][n2]
* ...
* w[nP][1] w[nP][2] ... w[nP][n1]
* h[nP][1] h[nP][2] ... h[nP][n2]
*
* d
* w_bar[1] w_bar[2] ... w_bar[d]
* h_bar[1] w_bar[2] ... w_bar[d]
*
* m1 m2 ... mP
* W[1][1] W[1][2] ... W[1][m]
* H[1][1] H[1][2] ... H[1][m]
* c[1][1] c[1][2] ... c[1][m]
* ...
* W[mP][1] W[mP][2] ... W[mP][m]
* H[mP][1] H[mP][2] ... H[mP][m]
* c[mP][1] c[mP][2] ... c[mP][m]
*/

void leEntrada(const string &arquivoEntrada){
	int i, j, p;
	try{
		ifstream file(arquivoEntrada.c_str());
		if(!file){
			cerr << "Erro ao abrir o arquivo!!!" << endl;
			exit(0);
		}
		//Número de periodos
		file >> P;
		//Número de itens por periodo
//		n = new int[P+1];
        n.resize((unsigned int)(P+1));
		for(p = 0; p < P; ++p)
            file >> n[p];
//        print(n, 1, P);
//        item = new Retangulo*[P+1];
//        for(p = 1; p<=P; ++p)
//            item[p] = new Retangulo[n[p]+1];
        item.resize((unsigned int)(P));
        for(p = 0; p < P; ++p)
            item[p].resize((unsigned int)(n[p]+1));

        for(p = 0; p < P; ++p){
            for(i = 1; i <= n[p];  ++i)
                file >> item[p][i].w;
            for(i = 1; i <= n[p];  ++i)
                file >> item[p][i].h;
            for(i = 1; i <= n[p];  ++i)
                item[p][i].p = p;
        }

        //Ordena os itens em cada periodo
        for(p = 0; p < P; ++p)
		sort(item[p].begin()+1, item[p].end(),
            [](Retangulo i1, Retangulo i2){return (i1.w < i2.w?true:(i1.w==i2.w?(i1.h<i2.h):false));});

        //Calcula os valores de n_q e o_q
        p1.resize((unsigned int)(P));
        n_q.resize((unsigned int)(P));
        o_q.resize((unsigned int)(P));

        for(p = 1; p < P; ++p){
            p1[p] = 1;
            for(i = 1; i <= n[p];  ++i){
                if(item[p][i].w == item[p][i+1].w && item[p][i].h == item[p][i+1].h)
                    continue;
                p1[p]++;
            }
            n_q[p].resize((unsigned int)(p1[p]+1));
            o_q[p].resize((unsigned int)(p1[p]+1));
            int o = 1;
            int k = 1;
            for(i = 1; i <= n[p];  ++i){
                if(item[p][i].w == item[p][i+1].w && item[p][i].h == item[p][i+1].h)
                    o++;
                else{
                    n_q[p][k++] = o;
                    o = 1;
                }
            }
            n_q[p][k] = o;
            o_q[p][1] = 0;
            for(i = 2; i <= p1[p]; ++i){
                o_q[p][i] = o_q[p][i-1] + n_q[p][i-1];
            }
        }

        file >> d;
        w_bar.resize((unsigned int)(d+1));
        h_bar.resize((unsigned int)(d+1));
		for(i = 1; i<=d; ++i)
			file >> w_bar[i];
		for(i = 1; i<=d; ++i)
			file >> h_bar[i];

        m.resize((unsigned int)(P+1));
		for(p = 0; p < P; ++p)
            file >> m[p];
        m[P] = 0;

        m_bar.resize((unsigned int)(P+1));
        int soma;
        for(p = 0; p <= P; ++p){
            soma = 0;
            for(int l = 1; l <= MIN(p,xi); l++){
                soma += (int)pow(2, l)*m[p-l];
            }
            m_bar[p] = m[p] + soma;
        }

        objeto.resize((unsigned int)(P+1));
        for(p = 0; p <= P; ++p)
            objeto[p].resize((unsigned int)(m_bar[p]+1));
        for(p = 0; p < P; ++p){
            for(j = 1; j <= m[p];  ++j){
                file >> objeto[p][j].w;
                WMax = (objeto[p][j].w>WMax?objeto[p][j].w:WMax);
            }
            for(j = 1; j <= m[p];  ++j){
                file >> objeto[p][j].h;
                HMax = (objeto[p][j].h>HMax?objeto[p][j].h:HMax);
            }
            for(j = 1; j <= m[p];  ++j)
                file >> objeto[p][j].c;
        }

        refObjPer.resize((unsigned int)(P+1));
        for( p=0; p <= P; p++){
            refObjPer[p].resize((unsigned int)(m_bar[p]+1));
            for( j = 1; j <= m[p]; ++j){
                refObjPer[p][j] = Par{p, j};
            }
        }
        for( p = 0; p < P; p++) {
            for (j = 1; j <= m_bar[p]; j++) {
                int j1 = m[p + 1] + 2 * j - 1;
                int j2 = j1 + 1;
                if(j1 <= m_bar[p+1]){
                    refObjPer[p + 1][j1] = refObjPer[p][j];
                     objeto[p + 1][j1].c = objeto[p][j].c;
                }
                if(j2 <= m_bar[p+1]){
                    refObjPer[p + 1][j2] = refObjPer[p][j];
                     objeto[p + 1][j2].c = objeto[p][j].c;
                }
            }
        }
        c_bar.resize((unsigned int)(P+1));
        for( p=0; p <= P; p++) {
            c_bar[p].resize((unsigned int) (m_bar[p] + 1));
        }
        for( p=0; p <= P; p++) {
            for (j = 1; j <= m_bar[p]; j++) {
                int obj = refObjPer[p][j].objeto;
                int per = refObjPer[p][j].periodo;
                c_bar[p][j] = objeto[per][obj].c;
            }
        }
	}
	catch(exception &e){
		cout << "Error: " << e.what() << endl;
	}
}

void setParametros(IloCplex &cplex){
    if(setLimiteMemoriaCplex) {
        cplex.setParam(IloCplex::WorkDir, ".");
        cplex.setParam(IloCplex::WorkMem, maxMemoriaCplex);
        cplex.setParam(IloCplex::NodeFileInd, mipStrategyFile);
    }
    if(setLimiteThreadsCplex)
        cplex.setParam(IloCplex::Threads, numThreadsCplex);

    cplex.setParam(IloCplex::TiLim, tempoLimiteCplex);
    cplex.setParam(IloCplex::ClockType, 2); //1 CPU time. 2 Wall clock time (total physical time elapsed).
    //cplex.setParam(IloCplex::EpOpt, 1e-9);
    cplex.setParam(IloCplex::EpAGap, 0.999999); //default: 1e-06
    cplex.setParam(IloCplex::EpGap, 0.0); //default: 1e-04
    //cplex.setParam(IloCplex::NumericalEmphasis, true);
    //cplex.setParam(IloCplex::PreInd, false);
    //cplex.setParam(IloCplex::AdvInd, 2);
    //cplex.setParam(IloCplex::Symmetry, 0);
    //cplex.setParam(IloCplex::Threads, THREADS);

    //cplex.setParam(IloCplex::MIPEmphasis, 2);
}



inline bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

void adicionaPontoDeCorte(vetor(Ponto)& v, Ponto& pt){
    int tam = (int) v.size();
    for(int i=0; i<tam; ++i)
        if(v[i] == pt)
            return;
    v.push_back(pt);
}

bool cabeItemCatalogo(vetor(int) &w_bar, vetor(int) &h_bar, Retangulo& sobra){
    if(sobra.getArea() <= 0) return false;
    unsigned long tam = w_bar.size();
    for(unsigned long i=1; i<tam; ++i){
        if(w_bar[i] <= sobra.w && h_bar[i] <= sobra.h)
            return true;
    }
    return false;
}