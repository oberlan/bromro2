#include <iostream>
#include <fstream>
#include <ilcplex/ilocplex.h>
#include <sys/time.h>
#include <chrono>
#include <vector>
#include <algorithm>

#include "InstanciaPorPeriodo.h"
#include "Util.h"
#include "Objeto.h"
#include "Item.h"

#define DIRETORIO_EXEC "../../Instances/bromro2/"
#define INSTANCIA "2.txt"
#define CPLEX_GET(X) DOUBLE_TO_INT(cplex.getValue(X))
#define CPLEX_GETD(X) cplex.getValue(X)

using namespace std;

typedef std::numeric_limits< double > dbl;
typedef vetor(InstanciaPorPeriodo) Instancia;
typedef IloArray<IloNumVarArray> NumVarMatrix2D;
typedef IloArray<NumVarMatrix2D> NumVarMatrix3D;
typedef IloArray<IloNumArray>    FloatMatrix;
typedef IloArray<IloExprArray>   ExprMatrix;
enum Metodo {PDA=1, GULOSO, HR, RF};


ulong P = 4;
ulong xi;
Instancia instancia;
string nomeInstancia;
string nomeArquivoInstancia;
int HMax = 0, WMax = 0;
double tempoLimiteCplex;
bool imprimirLogCplex;
bool desenharSolucaoPDF;
bool preSolveCplex;
bool setLimiteThreadsCplex;
int numThreadsCplex;
bool setLimiteMemoriaCplex;
int maxMemoriaCplex;
vector< pair<string, string> > padrao;

vector< Objeto > ListaObjetos;
vector< Item > ListaItens;
vector< Item > ListaItensCatalogo;


vetor( vetor( vetor(Item) ) ) listaItensDoObjeto;
vetor( vetor( vetor(Item) ) ) listaItensDoObjetoCpy;
vetor(double) tempoIteracaoPDA;
vetor( vetor(ulong) ) numObjetosDisponiveis;
vetor(double) valorDeltaInicial;
bool deltaIniVariavelPorPeriodo;

vetor( vetor(int) ) numObjetosPeriodo;
vetor( vetor(int) ) numSobrasPeriodo;


inline void leEntrada();
inline void metodoGuloso();
inline bool resolveModelo(ulong periodo);
inline bool cabeAlgumItemNoObjeto(Objeto& obj, vetor(Item)& item);
inline void setParametros(IloCplex &cplex, bool solucaoInteira = true);
inline double custoObjetos(ulong periodo);
inline void leArgumentos(int argc, char** argv);

int main(int argc, char** argv){

    leArgumentos(argc, argv);
    leEntrada();

    numObjetosPeriodo.resize(P+1);
    numSobrasPeriodo.resize(P+1);

    metodoGuloso();

	return 0;
}

void leArgumentos(int argc, char** argv) {
    tempoLimiteCplex = 600; 
    imprimirLogCplex = false;
    desenharSolucaoPDF = false;
    preSolveCplex = true;
    setLimiteMemoriaCplex = true;
    maxMemoriaCplex = 12000;
    setLimiteThreadsCplex = false;
    nomeInstancia = string(DIRETORIO_EXEC) + string(INSTANCIA);
    deltaIniVariavelPorPeriodo = false;
    xi = 1;

    if(argc == 1) {
        nomeInstancia = string(DIRETORIO_EXEC) + string(INSTANCIA);
        tempoLimiteCplex = 60;
        imprimirLogCplex = true;
        desenharSolucaoPDF = true;
        setLimiteMemoriaCplex = true;
        maxMemoriaCplex = 6000;
        setLimiteThreadsCplex = true;
        numThreadsCplex = 1;
        preSolveCplex = true;
        deltaIniVariavelPorPeriodo = false;
        xi = 1;
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
        if (string(argv[i]) == "-logCplex") {
            imprimirLogCplex = true;
            continue;
        }
        if (string(argv[i]) == "-solPDF") {
            desenharSolucaoPDF = true;
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
        if(string(argv[i]) == "-deltaIni"){
            valorInicialDelta = atof(argv[i+1]);
        }
        if(string(argv[i]) == "-deltaIniVariavel"){
            deltaIniVariavelPorPeriodo = true;
        }
        if(string(argv[i]) == "-sigma"){
            valorSigma = atof(argv[i+1]);
        }
        if(string(argv[i]) == "-xi"){
            xi = atoi(argv[i+1]);
        }
    }
    xi = xi<0?0:xi;

    //separando o nome da instancia
    string::size_type loc  = nomeInstancia.find_last_of( "/", nomeInstancia.size() );
    string::size_type loc2 = nomeInstancia.find_last_of( ".", nomeInstancia.size() );
    if( !((loc != string::npos ) && (loc2 != string::npos)) ){
        cout << "ERRO!" << endl;
        exit(0);
    }
    nomeArquivoInstancia = string("");
    nomeArquivoInstancia.append(nomeInstancia, loc+1, loc2-loc-1);

}

void leEntrada() {
    ulong i, s, j, j1, k;
    int minw=INF, minh=INF;
    try {
        ifstream file(nomeInstancia.c_str());
        if(!file){
            cerr << "Arquivo " << nomeInstancia << " nao encontrado" << endl;
            exit(0);
        }
        //Número de periodos
        file >> P;
        instancia.resize((unsigned int)P+2);
        //Número de itens por periodo
        for(s = 1; s<=P; ++s)
            file >> instancia[s].n;

        //item.resize((unsigned long)(P+1));
        for(s = 1; s<=P; ++s)
            instancia[s].item.resize((unsigned int)(instancia[s].n+1));

        for(s = 1; s<=P; ++s){
            instancia[s].periodo = s;
            for(j = 1; j <= instancia[s].n;  ++j) {
                file >> instancia[s].item[j].w;
                minw = (instancia[s].item[j].w<minw?instancia[s].item[j].w:minw);
            }
            for(j = 1; j <=instancia[s].n; ++j) {
                file >> instancia[s].item[j].h;
                minh = (instancia[s].item[j].h<minh?instancia[s].item[j].h:minh);
            }
            for(j = 1; j <=instancia[s].n; ++j)
                instancia[s].item[j].id = s;
        }

        //Ordena os itens pela área
        for(s = 1; s<=P; ++s) {
            sort(instancia[s].item.begin()+1, instancia[s].item.end(),
                 [](Item i1, Item i2){return (i1.getArea() > i2.getArea()?true:(i1.getArea()==i2.getArea()?(i1.w>i2.w):false));});
        }

        //Calcula os valores de n_q e o_q
        for(s = 1; s<=P; ++s){
            instancia[s].p = 1;
            for(j = 1; j <instancia[s].n; ++j){
                if(instancia[s].item[j] == instancia[s].item[j+1])
                    continue;
                instancia[s].p++;
            }
            instancia[s].n_q.resize((unsigned long)(instancia[s].p+1));
            instancia[s].o_q.resize((unsigned long)(instancia[s].p+1));
            ulong o = 1;
            k = 1;
            for(j = 1; j <instancia[s].n; ++j){
                if(instancia[s].item[j] == instancia[s].item[j+1])
                    o++;
                else{
                    instancia[s].n_q[k++] = o;
                    o = 1;
                }
            }

            instancia[s].n_q[k] = o;
            instancia[s].o_q[1] = 0;
            for(j = 2; j <=instancia[s].p; ++j){
                instancia[s].o_q[j] = instancia[s].o_q[j -1] + instancia[s].n_q[j -1];
            }
        }

        //Itens do catálogo usados apenas no último período
        file >> instancia[P].d;
        instancia[P].w_bar.resize((unsigned long)(instancia[P].d+1));
        instancia[P].h_bar.resize((unsigned long)(instancia[P].d+1));
        for(j = 1; j <=instancia[P].d; ++j)
            file >> instancia[P].w_bar[j];
        for(j = 1; j <=instancia[P].d; ++j)
            file >> instancia[P].h_bar[j];

        for(s = 1; s<P; ++s){
            instancia[s].d = instancia[P].d;
            instancia[s].w_bar.resize((unsigned long)(instancia[s].d+1));
            instancia[s].h_bar.resize((unsigned long)(instancia[s].d+1));
            for(i = 1; i<=instancia[P].d; ++i){
                instancia[s].w_bar[i] = (instancia[P].w_bar[i]);
                instancia[s].h_bar[i] = (instancia[P].h_bar[i]);
            }
        }

        for(s = 1; s<=P; ++s) {
            file >> instancia[s].m;
        }
        instancia[1].objeto.resize((unsigned long) (instancia[1].m + 1));
        instancia[1].m_bar = instancia[1].m;
        instancia[P + 1].m = 0;
        ulong soma;
        for(s = 2; s <= P + 1; ++s) {
            soma = 0;
            for(ulong l = 1; l <= MIN(s,xi); l++){
                soma += (ulong)pow(2, l)*instancia[s - l].m;
            }
            instancia[s].m_bar = instancia[s].m + soma;
            instancia[s].objeto.resize((unsigned long) (instancia[s].m_bar + 1));
        }

        for(s = 1; s <= P + 1; ++s){
//            cerr << "Periodo: " << s << endl;
            if (s <= P) {
                for (j = 1; j <= instancia[s].m; ++j) {
                    file >> instancia[s].objeto[j].w;
                    WMax = (instancia[s].objeto[j].w > WMax ? instancia[s].objeto[j].w : WMax);
                }
                for (j = 1; j <= instancia[s].m; ++j) {
                    file >> instancia[s].objeto[j].h;
                    HMax = (instancia[s].objeto[j].h > HMax ? instancia[s].objeto[j].h : HMax);
                }
                for (j = 1; j <= instancia[s].m; ++j) {
                    file >> instancia[s].objeto[j].c;
                    instancia[s].objeto[j].c_bar = instancia[s].objeto[j].c;
                }
            }
            for(j = 1; j <= instancia[s].m; ++j) {
                instancia[s].objeto[j].id = Par{s, j};
                instancia[s].objeto[j].codigo = instancia[s].objeto[j].codigoObjetoOrigem = Par{s, j};
                instancia[s].objeto[j].areaSobraDireita = instancia[s].objeto[j].w*instancia[s].objeto[j].h;
                instancia[s].objeto[j].areaSobraSuperior = instancia[s].objeto[j].w*instancia[s].objeto[j].h;
                instancia[s].objeto[j].delta = Delta(valorInicialDelta, valorInicialDelta, valorInicialDelta, valorInicialDelta);
            }
            if(s > 1) {
                for (j = instancia[s].m + 1, j1 = 1; j <= instancia[s].m_bar; j += 2, j1++) {
                    instancia[s].objeto[j].id = instancia[s - 1].objeto[j1].id;
                    instancia[s].objeto[j].c = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j].c_bar = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j].codigo = Par{s, j};
                    instancia[s].objeto[j].codigoObjetoOrigem = Par{s - 1, j1};
                    instancia[s].objeto[j].ehSobra = true;
                    instancia[s].objeto[j].delta = Delta(valorInicialDelta, valorInicialDelta, valorInicialDelta,
                                                         valorInicialDelta);

                    instancia[s].objeto[j + 1].id = instancia[s - 1].objeto[j1].id;
                    instancia[s].objeto[j + 1].c = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j + 1].c_bar = instancia[s - 1].objeto[j1].c;
                    instancia[s].objeto[j + 1].codigo = Par{s, j + 1};
                    instancia[s].objeto[j + 1].codigoObjetoOrigem = Par{s - 1, j1};
                    instancia[s].objeto[j + 1].ehSobra = true;
                    instancia[s].objeto[j + 1].delta = Delta(valorInicialDelta, valorInicialDelta, valorInicialDelta,
                                                             valorInicialDelta);

                }
            }
        }

        instancia[P+1].m = 0;
        instancia[P+1].m_bar = instancia[P+1].m + 2 * instancia[P].m_bar; 
        instancia[P+1].objeto.resize((unsigned long) (instancia[P+1].m_bar + 1));

    }catch(exception &e){
        cout << "Error: " << e.what() << endl;
    }
}

void metodoGuloso(){
        chrono::time_point<std::chrono::system_clock> inicioSim, fimSim;
        ulong s, j;
        bool encontrouSolucaoViavel;

        //Calcula a soma do custo dos objetos de cada tipo em cada período
        long double PESO = 0;
        for (s = 1; s <= P; s++) {
            for (j = 1; j <= instancia[s].m; ++j) {
                PESO += instancia[s].objeto[j].getCusto();
            }
        }


        inicioSim = chrono::system_clock::now();

#ifdef IMPRIMIR_SAIDA_SIMULACAO
        cerr << "------------------------- EXECUTANDO -------------------------" << endl;
#endif

        bool modeloViavel;
        encontrouSolucaoViavel = true;
        listaItensDoObjeto.clear();
        listaItensDoObjeto.resize(P + 1);
        double limiteCustoObj = INF;
        for (s = 1; s <= P; s++) {

#ifdef IMPRIMIR_SAIDA_SIMULACAO
            cerr << setw(16) << s;
#endif
            listaItensDoObjeto[s].resize(instancia[s].m + 1);
            modeloViavel = resolveModelo(s);
            if (!modeloViavel) {
                encontrouSolucaoViavel = false;
                break;
            }
        }
        if (!encontrouSolucaoViavel) {
            cerr << endl << "*** Solucao não encontrada ***" << endl;
            return;
        }

        long double custoTotal = 0;
        long double valorSobrasRemanescentes = 0;
        for (s = 1; s <= P; s++) {
            for (j = 1; j <= instancia[s].m; ++j) {
                if (instancia[s].objeto[j].utilizado) {
                    custoTotal += instancia[s].objeto[j].getCusto();
                    valorSobrasRemanescentes += instancia[s].objeto[j].getValorSobraRemanescente();
                }
            }
        }
        long double FO = PESO * custoTotal - valorSobrasRemanescentes;
        fimSim = chrono::system_clock::now();
        chrono::duration<double> tempoGasto = fimSim - inicioSim;

        cerr << fixed << setprecision(2);
        cerr << nomeArquivoInstancia << "\t" << xi << "\t" << FO
             << "\t"<< custoTotal << "\t" << valorSobrasRemanescentes << "\t" << tempoGasto.count() << endl;

}

bool resolveModelo(ulong periodo) {
    double inicio = MyClock();
    IloEnv env;
    ulong i, j, k, l, s, j1, j2, il;
    //int m = instancia[periodo].m;
    ulong numObjetos = 0;
    ulong numSobrasValidas = 0;
    ulong m = 0;
    ulong m_bar = instancia[periodo].m_bar;
    ulong n = instancia[periodo].n;
    ulong p = instancia[periodo].p;
    ulong d = instancia[periodo].d;
    vetor(ulong) &n_q = instancia[periodo].n_q;
    vetor(ulong) &o_q = instancia[periodo].o_q;
    vetor(ulong) &w_bar = instancia[periodo].w_bar;
    vetor(ulong) &h_bar = instancia[periodo].h_bar;
    vetor(Item) &item = instancia[periodo].item;
#ifdef IMPRIMIR_INFO
    cout << "--- Periodo: " << periodo << endl;
    cout << "Tempo Limite: " << tempoLimiteCplex << endl;
#endif
    vetor(Objeto) objeto;
    objeto.push_back(Objeto());
    HMax = 0;
    WMax = 0;
#ifdef IMPRIMIR_INFO
    cout << "\tObjetos:\n";
#endif
    //Copia apenas os objetos que possuiem dimensões válidas e que podem conter ao menos um item do período
    for(j=1; j<=m_bar; ++j){
        if(instancia[periodo].objeto[j].ehObjetoValido() && periodo - instancia[periodo].objeto[j].id.periodo <= xi){
            Objeto obj = instancia[periodo].objeto[j];
#ifdef IMPRIMIR_INFO
            cout << "\t\t";
            cout << obj.id << " ";
            cout << obj.w << " x " << obj.h
                 << " " << (obj.ehSobra?"Sobra":"Objeto");
#endif
            if(cabeAlgumItemNoObjeto(obj, item)) {
#ifdef IMPRIMIR_INFO
                cout << " (Incluido)\n";
#endif
                objeto.push_back(obj);
                m++;

                WMax = (obj.w>WMax?obj.w:WMax);
                HMax = (obj.h>HMax?obj.h:HMax);
                if(obj.ehSobra)
                    numSobrasValidas++;
                else
                    numObjetos++;
            }
            else{
#ifdef IMPRIMIR_INFO
                cout << " (Nao incluido)\n";
#endif
                if(periodo < P && periodo - obj.id.periodo < xi){
                    if(obj.id.periodo != periodo) {
                        instancia[periodo + 1].adicionaObjetoSobraDireita(obj);
                        listaItensDoObjeto[obj.id.periodo][obj.id.objeto].push_back(Item(obj.w, obj.h, obj.x, obj.y, periodo, true));
                    }
                }
                else{
                    if((obj.id.periodo != P)) { //Objeto originado de sobra mas não utilizado
                        listaItensDoObjeto[obj.id.periodo][obj.id.objeto].push_back(Item(obj.w, obj.h, obj.x, obj.y, periodo, true));
                    }
                }
            }

        }
    }

    numObjetosPeriodo[periodo].push_back(numObjetos);
    numSobrasPeriodo[periodo].push_back(numSobrasValidas);
#ifdef IMPRIMIR_INFO
    cout << "\tItens: " << endl;
    for(i = 1; i <= n; ++i)
        cout << "\t\t" << item[i].w << " x " << item[i].h << endl;
    cout << "\tItens catálogo: " << endl;
    for(i = 1; i <= d; ++i)
        cout << "\t\t" << w_bar[i] << " x " << h_bar[i] << endl;
    cout << endl;
    cout << "WMax = " << WMax << "\n" << "HMax = " << HMax << endl;
#endif

    try {
        IloModel modelo(env, "Cutting");
        IloCplex cplex(env);

        m_bar = 2*m + 1; /******************************************************************/

        auto L = (unsigned int) (floor(log2(WMax)) + 1);
        vector<int> pot(L+1);
        for(k=1; k<=L; ++k)
            pot[k] = (1<<(k-1));

        //Variáveis de decisão
        IloNumVarArray u(env, m + 1, 0, 1, ILOBOOL);
        NumVarMatrix2D v(env, n + 1);
        IloNumVarArray x(env, n + 1, 0.0, IloInfinity, ILOFLOAT);
        IloNumVarArray y(env, n + 1, 0.0, IloInfinity, ILOFLOAT);
        NumVarMatrix2D pi(env, n + 1);
        NumVarMatrix2D tau(env, n + 1);
        IloNumVarArray r(env, m + 1, 0, IloInfinity, ILOFLOAT);
        IloNumVarArray t(env, m + 1, 0, IloInfinity, ILOFLOAT);
        NumVarMatrix2D W_bar(env, 2);
        NumVarMatrix2D H_bar(env, 2);
        IloNumVarArray eta(env, m + 1, 0, 1, ILOBOOL);

        NumVarMatrix2D theta(env, m_bar+1);
        IloNumVarArray gamma(env, m_bar+1, 0.0, IloInfinity, ILOFLOAT);
        NumVarMatrix2D mu(env, m_bar+1);
        NumVarMatrix2D q(env, m_bar+1);
        NumVarMatrix2D omega(env, m_bar+1);
        NumVarMatrix2D zeta(env, m_bar+1);


//        IloNumVarArray alphaR(env, m + 1, 0, IloInfinity, ILOFLOAT);
//        IloNumVarArray alphaT(env, m + 1, 0, IloInfinity, ILOFLOAT);
//        IloNumVarArray betaR(env, m + 1, 0, IloInfinity, ILOFLOAT);
//        IloNumVarArray betaT(env, m + 1, 0, IloInfinity, ILOFLOAT);
//        IloNumVarArray gammaT(env, m + 1, 0, IloInfinity, ILOFLOAT);
//        IloNumVarArray gammaR(env, m + 1, 0, IloInfinity, ILOFLOAT);
//        NumVarMatrix2D theta(env, Nbits + 1);
//        NumVarMatrix2D omega(env, Nbits + 1);
//        NumVarMatrix2D alphaR_bar(env, d + 1);
//        NumVarMatrix2D alphaT_bar(env, d + 1);
//        NumVarMatrix2D betaR_bar(env, d + 1);
//        NumVarMatrix2D betaT_bar(env, d + 1);
        for (i = 1; i <= n; ++i) {
            v[i] = IloNumVarArray(env, m + 1, 0, 1, ILOBOOL);
            pi[i] = IloNumVarArray(env, n + 1, 0, 1, ILOBOOL);
            tau[i] = IloNumVarArray(env, n + 1, 0, 1, ILOBOOL);
        }
//        for (i = 1; i <= d; ++i) {
//            alphaR_bar[i] = IloNumVarArray(env, m + 1, 0, 1, ILOBOOL);
//            alphaT_bar[i] = IloNumVarArray(env, m + 1, 0, 1, ILOBOOL);
//            betaR_bar[i] = IloNumVarArray(env, m + 1, 0, 1, ILOBOOL);
//            betaT_bar[i] = IloNumVarArray(env, m + 1, 0, 1, ILOBOOL);
//        }
        for(j = 1; j <= m_bar; ++j){
            theta[j] = IloNumVarArray(env, L+1, 0, 1, ILOBOOL);
            omega[j] = IloNumVarArray(env, L+1, 0.0, IloInfinity, ILOFLOAT);
            zeta[j] = IloNumVarArray(env, d+1, 0, 1, ILOBOOL);
            mu[j] = IloNumVarArray(env, d+1, -IloInfinity, IloInfinity, ILOFLOAT);
            q[j] = IloNumVarArray(env, d+1, 0, 1, ILOBOOL);
        }
        W_bar[0] = IloNumVarArray(env, m+1, 0.0, IloInfinity, ILOFLOAT);
        H_bar[0] = IloNumVarArray(env, m+1, 0.0, IloInfinity, ILOFLOAT);
        W_bar[1] = IloNumVarArray(env, m_bar+1, 0.0, IloInfinity, ILOFLOAT);
        H_bar[1] = IloNumVarArray(env, m_bar+1, 0.0, IloInfinity, ILOFLOAT);

        //Número de objeto utilizados
        IloExpr objetosUsados(env);
        //Custo dos objetos utilizados
        IloExpr custoObjetos(env);
        //Valor da proporção de utilização das sobras dos objetos
        IloExpr deltasObjetos(env);
        //Custo dos objetos de entrada utilizados
        IloExpr custoObjetosEntrada(env);
        //Valor das sobras dos objetos utilizados
        IloExpr valorSobras(env);
        //Valor das sobras dos objetos usados mais as sobras provenientes de outros períodos não utilizadas
        IloExpr totalValorSobras(env);
        //IloExpr leftoversValue2(env);
        IloExpr numeroSobras(env);
        IloExpr objetivoOriginal(env);
        long double PESO = 0.0;


        PESO = 0;
        for (j = 1; j <= m; ++j) {
            j1 = 2 * j - 1;
            j2 = 2 * j;
            if (!objeto[j].ehSobra) {
                custoObjetos += (objeto[j].getCusto() * u[j]);
            }
            objetosUsados += u[j];
            if (periodo - objeto[j].id.periodo < xi) {
                valorSobras += (objeto[j].getCustoSobraPorUnidadeArea() * (gamma[j1] + gamma[j2]));
            }
            PESO += objeto[j].getCustoSobra();
        }
        objetivoOriginal = PESO * custoObjetos - valorSobras;


        //=====================================================Funcao Objetivo
        modelo.add(IloMinimize(env, objetivoOriginal));

        //=====================================================Restricoes
        IloRangeArray restricoes(env);


        //=====================================================
        /*
        forall(i in ITENS)
            sum(j in OBJETOS) v[i][j] == 1;
        */
        for (i = 1; i <= n; ++i) {
            IloExpr somatorio(env);
            for (j = 1; j <= m; ++j)
                somatorio += v[i][j];
            restricoes.add(somatorio == 1);
            somatorio.end();
        }

        //=====================================================
        /*
        forall(i in ITENS, j in OBJETOS)
            u[j] >= v[i][j];
        */
        for (i = 1; i <= n; ++i)
            for (j = 1; j <= m; ++j)
                restricoes.add(u[j] - v[i][j] >= 0);

        //=====================================================
        for (j = 1; j <= m; ++j) {
            IloExpr somatorio(env);
            for (i = 1; i <= n; ++i)
                somatorio += v[i][j];
            restricoes.add(u[j] - somatorio <= 0);
            somatorio.end();
        }

        //Atribui o tamanho de cada objeto
        for (j = 1; j <= m; ++j) {
            restricoes.add(W_bar[0][j] - objeto[j].w == 0);
            restricoes.add(H_bar[0][j] - objeto[j].h == 0);
        }

        //=====================================================
        for (j = 1; j <= m; ++j) {
            restricoes.add(t[j] - H_bar[0][j] <= 0);
            restricoes.add(r[j] - W_bar[0][j] <= 0);
        }

        //=====================================================
        for(j = 1; j <= m; ++j){
            j1 = 2 * j - 1;
            j2 = 2 * j;
            restricoes.add(H_bar[1][j1] - (u[j] * HMax) <= 0);
            restricoes.add(t[j] - (1 - u[j]) * HMax - H_bar[1][j1] <= 0);
            restricoes.add(H_bar[1][j1] - (t[j] + (1 - u[j]) * HMax) <= 0);

            restricoes.add(W_bar[1][j1] - (u[j] * WMax) <= 0);
            restricoes.add(W_bar[0][j] - r[j] - (1 - eta[j]) * WMax - (1 - u[j]) * WMax - (W_bar[1][j1]) <= 0);
            restricoes.add(W_bar[1][j1] - (W_bar[0][j] - r[j] + (1 - eta[j]) * WMax + (1 - u[j]) * WMax) <= 0);
            restricoes.add(W_bar[0][j] - eta[j] * WMax - (1 - u[j]) * WMax - (W_bar[1][j1]) <= 0);
            restricoes.add(W_bar[1][j1] - (W_bar[0][j] + eta[j] * WMax + (1 - u[j]) * WMax) <= 0);


            restricoes.add(W_bar[1][j2] - (u[j] * WMax) <= 0);
            restricoes.add(r[j] - (1 - u[j]) * WMax - (W_bar[1][j2]) <= 0);
            restricoes.add(W_bar[1][j2] - (r[j] + (1 - u[j]) * WMax) <= 0);

            restricoes.add(H_bar[1][j2] - (u[j] * HMax) <= 0);
            restricoes.add(H_bar[0][j] - (1 - eta[j]) * HMax - (1 - u[j]) * HMax - (H_bar[1][j2]) <= 0);
            restricoes.add(H_bar[1][j2] - (H_bar[0][j] + (1 - eta[j]) * HMax + (1 - u[j]) * HMax) <= 0);
            restricoes.add(H_bar[0][j] - t[j] - eta[j] * HMax - (1 - u[j]) * HMax - (H_bar[1][j2]) <= 0);
            restricoes.add(H_bar[1][j2] - (H_bar[0][j] - t[j] + eta[j] * HMax + (1 - u[j]) * HMax) <= 0);

        }

        //=====================================================
        for (i = 1; i <= n; ++i) {
            restricoes.add(x[i] - item[i].w / 2.0 >= 0);
            restricoes.add(y[i] - item[i].h / 2.0 >= 0);
        }
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= m; ++j) {
                restricoes.add(x[i] + item[i].w / 2.0 - (W_bar[0][j] - r[j] + (WMax) * (1 - v[i][j])) <= 0);
                restricoes.add(y[i] + item[i].h / 2.0 - (H_bar[0][j] - t[j] + (HMax) * (1 - v[i][j])) <= 0);
            }
        }

        //=====================================================
        for(j = 1; j <= m; ++j)
            for(i = 1; i <= n; ++i)
                for(il = i+1; il <= n; ++il){
                    restricoes.add( x[i] - x[il]  - (0.5*(item[i].w + item[il].w) - WMax*((1-v[i][j]) + (1-v[il][j]) + pi[i][il] + tau[i][il])) >= 0);
                    restricoes.add(-x[i] + x[il]  - (0.5*(item[i].w + item[il].w) - WMax*((1-v[i][j]) + (1-v[il][j]) + pi[i][il] + (1 - tau[i][il]))) >= 0);

                    restricoes.add( y[i] - y[il]  - (0.5*(item[i].h + item[il].h) - HMax*((1-v[i][j]) + (1-v[il][j]) + (1 - pi[i][il]) + tau[i][il])) >= 0);
                    restricoes.add(-y[i] + y[il]  - (0.5*(item[i].h + item[il].h) - HMax*((1-v[i][j]) + (1-v[il][j]) + (1 - pi[i][il]) + (1 - tau[i][il]))) >= 0);
                }

        //=====================================================
        // Restrição de simetria de itens iguais
        for (ulong q = 1; q <= p; ++q)
            for (i = o_q[q] + 1; i <= o_q[q] + n_q[q]; ++i)
                for (ulong i1 = i + 1; i1 <= o_q[q] + n_q[q]; ++i1)
                    for (j = 1; j <= m; ++j) {
                        restricoes.add(x[i1] - item[i1].w / 2.0 - (x[i] + item[i].w / 2.0 - WMax * (1 - v[i][j]) - WMax * (1 - v[i1][j]) - WMax * pi[i][i1]) >= 0);
                        restricoes.add(y[i1] - item[i1].h / 2.0 - (y[i] + item[i].h / 2.0 - HMax * (1 - v[i][j]) - HMax * (1 - v[i1][j]) - HMax * (1 - pi[i][i1])) >= 0);
                    }
        for (ulong q = 1; q <= p; ++q)
            for (i = o_q[q] + 1; i <= o_q[q] + n_q[q]; ++i)
                for (ulong i1 = o_q[q] + n_q[q] + 1; i1 <= n; ++i1)
                    for (j = 1; j <= m; ++j) {
                        restricoes.add(x[i1] + item[i1].w / 2.0 - (x[i] - item[i].w / 2.0 + WMax * (1 - v[i][j]) + WMax * (1 - v[i1][j]) + WMax * pi[i][i1] + WMax * tau[i][i1]) <= 0);
                        restricoes.add(x[i1] - item[i1].w / 2.0 - (x[i] + item[i].w / 2.0 - WMax * (1 - v[i][j]) - WMax * (1 - v[i1][j]) - WMax * pi[i][i1] - WMax * (1 - tau[i][i1])) >= 0);
                        restricoes.add(y[i1] + item[i1].h / 2.0 - (y[i] - item[i].h / 2.0 + HMax * (1 - v[i][j]) + HMax * (1 - v[i1][j]) + HMax * (1 - pi[i][i1]) + HMax * tau[i][i1]) <= 0);
                        restricoes.add(y[i1] - item[i1].h / 2.0 - (y[i] + item[i].h / 2.0 - HMax * (1 - v[i][j]) - HMax * (1 - v[i1][j]) - HMax * (1 - pi[i][i1]) - HMax * (1 - tau[i][i1])) >= 0);
                    }


        //=====================================================
        for(j = 1; j <= m_bar; ++j){
            IloExpr somatorio(env);
            for(l = 1; l <= L; ++l){
                somatorio += pow(2, l-1)*theta[j][l];
            }
            restricoes.add(W_bar[1][j] - somatorio == 0);
        }

        //=====================================================
        for(j = 1; j <= m_bar; ++j) {
            for (l = 1; l <= L; ++l) {
                restricoes.add(omega[j][l] - H_bar[1][j] <= 0);
                restricoes.add(H_bar[1][j] - (1 - theta[j][l]) * HMax - omega[j][l] <= 0);
                restricoes.add(omega[j][l] - (theta[j][l] * HMax) <= 0);
            }
        }

        //=====================================================
        for(j = 1; j <= m_bar; ++j) {
            for (i = 1; i <= d; ++i) {
                restricoes.add(w_bar[i] - (W_bar[1][j] + WMax * (1 - zeta[j][i])) <= 0);
                restricoes.add(h_bar[i] - (H_bar[1][j] + HMax * (1 - zeta[j][i])) <= 0);
            }
        }

        //=====================================================
        for(j = 1; j <= m_bar; ++j) {
            IloExpr somatorio(env);
            for (l = 1; l <= L; ++l) {
                somatorio += pow(2, l - 1) * omega[j][l];
            }
            restricoes.add(gamma[j] - somatorio <= 0);
        }
        for(j = 1; j <= m_bar; ++j) {
            IloExpr somatorio(env);
            for (i = 1; i <= d; ++i) {
                somatorio += zeta[j][i];
            }
            restricoes.add(gamma[j] - somatorio * WMax * HMax <= 0);
        }

        modelo.add(restricoes);
        cplex.extract(modelo);

        if (!imprimirLogCplex)
            cplex.setOut(env.getNullStream());

        //setParametros(cplex, periodo == P);
        setParametros(cplex);

        cplex.solve();

        bool finished;
#ifdef IMPRIMIR_INFO
        cout << "[###] Solver status: " << cplex.getStatus() << endl;
#endif
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            finished = true;
        } else if (cplex.getStatus() == IloAlgorithm::Feasible) {
            finished = true;
        } else {
#ifdef IMPRIMIR_INFO
            cerr << "[###] Solver status: " << cplex.getStatus() << endl;
#endif
            finished = false;
        }
        if(!finished) {
            env.end();
            return finished;
        }

        IloInt accumulatedIterations = cplex.getNiterations();
        IloInt accumulatedNodes = cplex.getNnodes();
        IloNum fObjBestBound = cplex.getBestObjValue();
        IloNum fObjValue = cplex.getObjValue();
        double cplexGap = cplex.getMIPRelativeGap();//(fObjValue - fObjBestBound) / (0.00000001 + fObjValue);

//        cout << "custoObjetos: " << cplex.getValue(custoObjetos) << endl;
//        cout << "ValorDeltas : " << cplex.getValue(deltasObjetos) << endl;
//        cout << "Sobras      : " << cplex.getValue(valorSobras) << endl;
//        cout << "PESO        : " << PESO << endl;
//        cout << "FO          : " << cplex.getValue(objetivoOriginal) << endl;
//        PAUSE;

#ifdef IMPRIMIR_INFO
        cout.precision(dbl::digits10);
        cout << setprecision(2) << fixed;
        cout << "[###] Best objective bound: " << fObjBestBound << endl;
        cout << "[###] Cplex GAP: " << fabs(cplexGap) * 100.0 << "%" << endl;
        cout << "[###] Objective obtained: " << fObjValue << endl;
//        cout << "[###] MIP Iterations: " << accumulatedIterations << endl;
//        cout << "[###] B&B Nodes: " << accumulatedNodes << endl << endl;

//        cout << "----------------SOLUTION FOUND----------------\n";
       // cout << "Objects cost: " << CPLEX_GETD(custoObjetos) << endl;
        cout << "Custo dos objetos: " << (periodo!=P?CPLEX_GETD(custoObjetosEntrada):CPLEX_GETD(custoObjetos)) << endl;
        cout << "Valor das sobras: " << (periodo!=P?CPLEX_GETD(totalValorSobras):CPLEX_GETD(valorSobras)) << endl;
        //cout << "Leftovers value: " << CPLEX_GETD(valorSobras) << endl;
        cout << "Objects used: " << CPLEX_GETD(objetosUsados) << endl;
        cout << "Leftovers: " << CPLEX_GETD(numeroSobras) << endl;
        cout << "Cplex time: " << (MyClock() - inicio) / CLOCKS_PER_SEC << endl;
#endif
        double cte = 0.0;
        double FO;

#ifdef IMPRIMIR_SAIDA_SIMULACAO
//            for (j = 1; j <= m; ++j) cte += objeto[j].getCustoSobra();
//
//            if(periodo!=P)
//                FO = cte * CPLEX_GETD(custoObjetosEntrada) - CPLEX_GETD(totalValorSobras);
//            else
//                FO = cte * CPLEX_GETD(custoObjetos) - CPLEX_GETD(valorSobras);
//
//            cerr << setw(0) << setprecision(4) << fixed;
//            cerr << setw(17) << (periodo!=P?CPLEX_GETD(custoObjetosEntrada):CPLEX_GETD(custoObjetos))
//                 << setw(16) << (periodo!=P?CPLEX_GETD(totalValorSobras):CPLEX_GETD(valorSobras))
//                 << setw(15) << FO << setw(15) << (MyClock() - inicio) / CLOCKS_PER_SEC
//                 << setw(12) << fabs(cplexGap) * 100 << "%"
//                 << endl;
//            cerr << setw(0);

        cerr << setw(0) << setprecision(4) << fixed;
        cerr << setw(17) << CPLEX_GETD(custoObjetos)
             << setw(16) << CPLEX_GETD(valorSobras)
             << setw(15) << fObjValue << setw(15) << (MyClock() - inicio) / CLOCKS_PER_SEC
             << setw(12) << fabs(cplexGap) * 100 << "%"
             << endl;
        cerr << setw(0);
#endif
        for (j = 1; j <= m; ++j) cte += objeto[j].getCustoSobra();

        if (periodo != P)
            FO = cte * CPLEX_GET(custoObjetosEntrada) - CPLEX_GET(totalValorSobras);
        else
            FO = cte * CPLEX_GET(custoObjetos) - CPLEX_GET(valorSobras);

        cout << setw(0) << setprecision(4) << fixed;
        cout << setw(3) << periodo;
        cout << setw(17) << (periodo != P ? CPLEX_GET(custoObjetosEntrada) : CPLEX_GET(custoObjetos))
             << setw(16) << (periodo != P ? CPLEX_GET(totalValorSobras) : CPLEX_GET(valorSobras))
             << setw(15) << FO << setw(15) << (MyClock() - inicio) / CLOCKS_PER_SEC
             << setw(12) << fabs(cplexGap) * 100 << "%"
             << setw(10) << numObjetos
             << setw(10) << numSobrasValidas
             << endl;
        cout << setw(0);
        if (periodo == P)
            cout << "Valor sobras remanescentes: " << CPLEX_GET(valorSobras) << endl;


        for(j = 1; j<=m; ++j) {
            if(CPLEX_GET(u[j]) == 1) {
                objeto[j].utilizado = true;
                ulong idx = objeto[j].codigo.objeto;//objeto[j].codigo%100;
                instancia[periodo].objeto[idx].utilizado = true;
//                cout << "Objeto: " << objeto[j].id << endl;
                for(i = 1; i<=n; ++i) {
                    if(CPLEX_GET(v[i][j]) == 1) {
                        ulong per = objeto[j].id.periodo;//objeto[j].id / 100;
                        ulong obj = objeto[j].id.objeto;//objeto[j].id % 100;
                        Item it = item[i];
                        it.x = CPLEX_GETD(x[i]) - item[i].w / 2.0;
                        it.y = CPLEX_GETD(y[i]) - item[i].h / 2.0;
                        it.x += objeto[j].x;
                        it.y += objeto[j].y;
//                        cout << setprecision(2) << fixed << i << ": (" << it.x << ", " << it.y  << ")\n";
                        listaItensDoObjeto[per][obj].push_back(it);
                    }
                }
            }
        }

        if(periodo < P) {
            for (j = 1; j <= m; ++j) {
                ulong per = objeto[j].id.periodo;
                ulong obj = objeto[j].id.objeto;
                /*********** Verifica se a sobra pode ser adicionada ************/
                if(periodo - objeto[j].id.periodo >= xi)
                    continue;
                if (CPLEX_GET(u[j]) == 1){
                    j1 = 2 * j - 1;
                    j2 = 2 * j;
                    ulong WsobraSuperior = CPLEX_GET(W_bar[1][j1]);
                    ulong HsobraSuperior = CPLEX_GET(H_bar[1][j1]);
                    ulong WsobraDireita  = CPLEX_GET(W_bar[1][j2]);
                    ulong HsobraDireita  = CPLEX_GET(H_bar[1][j2]);
                    ulong r_j = CPLEX_GET(r[j]);
                    ulong t_j = CPLEX_GET(t[j]);
                    double x_j, y_j;
                    if(WsobraSuperior*HsobraSuperior > 0){
                        x_j = objeto[j].x + 0;
                        y_j = objeto[j].y + (objeto[j].h - t_j);
                        Objeto novoObjeto(WsobraSuperior, HsobraSuperior, objeto[j].c, x_j, y_j, objeto[j].id, objeto[j].c_bar, true, objeto[j].codigo, objeto[j].codigoObjetoOrigem);
                        instancia[periodo + 1].adicionaObjetoSobraSuperior(novoObjeto);
                        listaItensDoObjeto[per][obj].push_back(Item(WsobraSuperior, HsobraSuperior, x_j, y_j, periodo, true));
                    }
                    if(WsobraDireita*HsobraDireita > 0){
                        x_j = objeto[j].x + (objeto[j].w - r_j);
                        y_j = objeto[j].y + 0;
                        Objeto novoObjeto(WsobraDireita, HsobraDireita, objeto[j].c, x_j, y_j, objeto[j].id, objeto[j].c_bar, true, objeto[j].codigo, objeto[j].codigoObjetoOrigem);
                        instancia[periodo + 1].adicionaObjetoSobraDireita(novoObjeto);
                        listaItensDoObjeto[per][obj].push_back(Item(WsobraDireita, HsobraDireita, x_j, y_j, periodo, true));
                    }
                }
                // Se o objeto é originado de uma sobra aproveitável, então a sobra é passada para o próximo período
                // Neste caso, o objeto é colocado na posição equivalente a sobra da direita
                else if (objeto[j].id.periodo != periodo) {
                    instancia[periodo + 1].adicionaObjetoSobraDireita(objeto[j]);
                    listaItensDoObjeto[per][obj].push_back(Item(objeto[j].w, objeto[j].h, objeto[j].x, objeto[j].y, periodo, true));
                }
            }
        }
//        cout << "Area das sobras aproveitaveis do período " << periodo << endl;
        //Atualiza a área da sobra aproveitável de cada objeto (de *entrada*) do período
        for (j = 1; j <= m; ++j) {
            j1 = 2 * j - 1;
            j2 = 2 * j;
            ulong idx = objeto[j].codigo.objeto;
            if (/*objeto[j].id / 100 == periodo &&*/ CPLEX_GET(u[j]) == 1){
                instancia[periodo].objeto[idx].areaSobraSuperior = CPLEX_GET(gamma[j1]);
                instancia[periodo].objeto[idx].areaSobraDireita = CPLEX_GET(gamma[j2]);
//                cout << "Objeto utilizado " << instancia[periodo].objeto[idx].id << ": "
//                     <<instancia[periodo].objeto[idx].areaSobraDireita << " + " << instancia[periodo].objeto[idx].areaSobraSuperior << endl;
            }
            else if(objeto[j].id.periodo != periodo){
                instancia[periodo].objeto[idx].areaSobraDireita = objeto[j].w*objeto[j].h;
                instancia[periodo].objeto[idx].areaSobraSuperior = 0;
//                cout << "Objeto nao-utilizado " << instancia[periodo].objeto[idx].id << ": "
//                     << instancia[periodo].objeto[idx].areaSobraDireita << " + " << instancia[periodo].objeto[idx].areaSobraSuperior << endl;
            }
        }

        //***Atualiza a área utilizada de objetos (originados de sobras?)***//
        //Contabilia a área dos itens cortados no objeto
        for(j = 1; j<=m; ++j) {
            if(CPLEX_GET(u[j]) == 1) {
                for(i = 1; i<=n; ++i) {
                    if((CPLEX_GET(v[i][j]) == 1) ){
                        //instancia[per].objeto[obj].areaSobraUtilizada += (item[i].w*item[i].h);
                        ulong idx = objeto[j].codigo.objeto;//objeto[j].codigo%100;
                        instancia[periodo].objeto[idx].areaUtilizada += (item[i].w*item[i].h);
                    }
                }
            }
        }

        if(periodo == P) {
            for (j = 1; j <= m; ++j) {
                ulong per = objeto[j].id.periodo;
                ulong obj = objeto[j].id.objeto;
                /*********** Verifica se a sobra pode ser adicionada ************/
                if(periodo - objeto[j].id.periodo >= xi)
                    continue;
                if ((CPLEX_GET(u[j]) == 1)) {
                    ulong r_j = CPLEX_GET(r[j]);
                    ulong t_j = CPLEX_GET(t[j]);
                    j1 = 2 * j - 1;
                    j2 = 2 * j;
                    ulong WsobraSuperior = CPLEX_GET(W_bar[1][j1]);
                    ulong HsobraSuperior = CPLEX_GET(H_bar[1][j1]);
                    ulong WsobraDireita = CPLEX_GET(W_bar[1][j2]);
                    ulong HsobraDireita = CPLEX_GET(H_bar[1][j2]);
                    double x_j, y_j;
                    ulong w, h;
                    //Contabilia a área das sobra aproveitáveis remanescentes dos objetos dos perídos anteriores à P
                    if (WsobraSuperior * HsobraSuperior > 0) {
                        x_j = objeto[j].x + 0;
                        y_j = objeto[j].y + (objeto[j].h - t_j);
                        listaItensDoObjeto[per][obj].push_back(
                                Item(WsobraSuperior, HsobraSuperior, x_j, y_j, periodo, true));
                    }
                    if (WsobraDireita * HsobraDireita > 0) {
                        x_j = objeto[j].x + (objeto[j].w - r_j);
                        y_j = objeto[j].y + 0;
                        listaItensDoObjeto[per][obj].push_back(
                                Item(WsobraDireita, HsobraDireita, x_j, y_j, periodo, true));
                    }
                    instancia[per].objeto[obj].areaSobraRemanescente +=
                            WsobraSuperior * HsobraSuperior + WsobraDireita * HsobraDireita;
                }
                else if((objeto[j].id.periodo != P)) { //Objeto originado de sobra mas não utilizado
                    listaItensDoObjeto[per][obj].push_back(Item(objeto[j].w, objeto[j].h, objeto[j].x, objeto[j].y, periodo, true));
                    instancia[per].objeto[obj].areaSobraRemanescente += objeto[j].w*objeto[j].h;
                }
            }
        }
#ifdef IMPRIMIR_INFO
        cout << endl;
#endif
        env.end();
        return finished;
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
    return false;
}

void setParametros(IloCplex &cplex, bool solucaoInteira){
    if(setLimiteMemoriaCplex) {
        cplex.setParam(IloCplex::WorkDir, ".");
        cplex.setParam(IloCplex::WorkMem, maxMemoriaCplex);
        cplex.setParam(IloCplex::NodeFileInd, 3);
    }
    if(setLimiteThreadsCplex)
        cplex.setParam(IloCplex::Threads, numThreadsCplex);

    cplex.setParam(IloCplex::TiLim, tempoLimiteCplex);
    cplex.setParam(IloCplex::ClockType, 2); //1 CPU time. 2 Wall clock time (total physical time elapsed).
    //cplex.setParam(IloCplex::EpOpt, 1e-9);
    if(solucaoInteira)
        cplex.setParam(IloCplex::EpAGap, 0.999999); //default: 1e-06
    else
        cplex.setParam(IloCplex::EpAGap, 1e-7);  //default: 1e-06
    cplex.setParam(IloCplex::EpGap, 0.0); //default: 1e-04
    cplex.setParam(IloCplex::PreInd, preSolveCplex);
    //cplex.setParam(IloCplex::NumericalEmphasis, true);
    //cplex.setParam(IloCplex::AdvInd, 2);
    //cplex.setParam(IloCplex::Symmetry, 0);

    /* Value	Meaning
     *   0	    Balance optimality and feasibility; default
     *   1	    Emphasize feasibility over optimality
     *   2	    Emphasize optimality over feasibility
     *   3	    Emphasize moving best bound
     *   4	    Emphasize finding hidden feasible solutions
     */
    //cplex.setParam(IloCplex::MIPEmphasis, 1);
}

bool cabeAlgumItemNoObjeto(Objeto& obj, vetor(Item)& item){
    ulong numItens = (ulong) item.size() - 1;
    for (ulong i = 1; i <= numItens; ++i) {
        if (item[i].w <= obj.w && item[i].h <= obj.h)
            return true;
    }
    return false;
}


ostream& operator<<(ostream &out, const TuplaCor &t){
    out << t.nome;
    return out;
}

