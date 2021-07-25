//
// Created by Oberlan Romão on 05/09/15.
//

#ifndef OBJETO_H
#define OBJETO_H
#include <iostream>

using namespace std;

#include "Delta.h"
#include "Util.h"

class Objeto {
public:
    ulong w, h;
    double c, x, y;
    Par id;
    double c_bar; //Custo por unidade de área da sobra
    bool ehSobra;
    Par codigo; //codigo do objeto de entrada
    Par codigoObjetoOrigem; //codigo do objeto ao qual o objeto foi gerado (no caso de sobra), se o for um objeto de entrada codigo == codigoObjetoOrigem
    double areaSobraDireita;
    double areaSobraSuperior;
    double areaSobraDireitaUtilizada; //Área utilizada na sobra da direita
    double areaSobraSuperiorUtilizada; //Área utilizada na sobra superior
    double areaUtilizada; //Área utilizada mais a área utilizada das sobras
    int areaSobraRemanescente;
    bool utilizado;
    bool deltaMudou;
    Delta delta; //Economia associada a área da sobra de cada objeto
    Objeto(ulong _w=0, ulong _h=0, double _c=0, double _x=0, double _y=0, Par _id={0,0}, double _c_bar=0, bool sobra=false, Par cod={0,0}, Par codObjOrigem={0,0}):
            w(_w),h(_h),c(_c),x(_x),y(_y),id(_id),c_bar(_c_bar),ehSobra(sobra),codigo(cod), codigoObjetoOrigem(codObjOrigem)
    {
        areaSobraDireita = w*h;
        areaSobraSuperior = w*h;
        areaUtilizada = 0.0;
        areaSobraDireitaUtilizada = 0.0;
        areaSobraSuperiorUtilizada = 0.0;
        utilizado = false;
        delta = Delta(valorInicialDelta, valorInicialDelta, valorInicialDelta, valorInicialDelta);//Delta(1.0, 1.0, 1.0, 1.0);
        deltaMudou = false;
        areaSobraRemanescente = 0;
    }
    bool operator<(const Objeto& obj) const {
        if(w < obj.w) return true;
        else if(w == obj.w) return (h < obj.h);
        else return false;
    }
    inline double getCusto() const {
        return w*h*c;
    }
    inline double getCustoPorUnidadeArea() const {
        return c;
    }
    inline double getCustoSobra() const {
        return w*h*c_bar;
    }
    inline double getCustoSobraPorUnidadeArea() const {
        return c_bar;
    }
    inline void print() const {
//        cout << w << " x " << h << "\t(" << x << ", " << y << ")\n";
        cout << codigo << "\t" << codigoObjetoOrigem << "\t" << c << "\t" << c_bar << "\t" << w << " x " << h << endl;
    }
    inline bool ehObjetoValido(){
        return w*h != 0;
    }
    inline void atualizaAreaUtilizada(){
        areaUtilizada += areaSobraDireitaUtilizada + areaSobraSuperiorUtilizada;
    }
    inline void reset(){
        w = 0;
        h = 0;
        x = 0;
        y = 0;
        areaSobraDireita = 0;
        areaSobraSuperior = 0;
        areaUtilizada = 0.0;
        areaSobraDireitaUtilizada = 0.0;
        areaSobraSuperiorUtilizada = 0.0;
        utilizado = false;
        deltaMudou = false;
        areaSobraRemanescente = 0;
    }
    inline void atualizaDeltas(){
        atualizaDeltaSobraDireita();
        atualizaDeltaSobraSuperior();
        atualizaDeltaSobras();
        atualizaDeltaObjeto();
    }
    inline ulong getArea() const{
        return w*h;
    }
    inline double getValorSobraRemanescente(){
        return c*areaSobraRemanescente;
    }
private:
    inline bool atualizaDeltaSobraDireita(){
        double tx;
        //cout << id << " " << codigo << " " << codigoObjetoOrigem << " D(";
        //cout << (int)areaSobraDireita << " " << (int)areaSobraDireitaUtilizada << " ";
        if(!DOUBLE_EQUALS(areaSobraDireita, 0.0))
            tx = areaSobraDireitaUtilizada/areaSobraDireita;
        else
#ifdef AREANULA_DELTAZERO
        tx = 0.0;
#else
        tx = delta.sobraDireita;
#endif

#ifndef PODE_ZERAR_DELTAS
        if(!DOUBLE_EQUALS(tx, delta.sobraDireita) && !DOUBLE_EQUALS(tx,0.0)) {
#else
        if(!DOUBLE_EQUALS(tx, delta.sobraDireita)) {
#endif
            //delta.sobraDireita = (DOUBLE_EQUALS(tx,0.0)?0.01:tx);
            //delta.sobraDireita = (DOUBLE_EQUALS(tx,0.0)?delta.sobraDireita:tx);
            //delta.sobraDireita = tx;
            delta.sobraDireita = delta.sobraDireita + (tx - delta.sobraDireita)*(ALPHA);//pow(valorSigma, iteracao);
            deltaMudou = true;
        }
        else
            deltaMudou = false;
        //cout << tx << ")" << endl;
        return deltaMudou;
    }
    inline bool atualizaDeltaSobraSuperior(){
        double tx;
        //cout << id << " " << codigo << " " << codigoObjetoOrigem << " T(";
        //cout << (int)areaSobraSuperior << " " << (int)areaSobraSuperiorUtilizada << " ";
        if(!DOUBLE_EQUALS(areaSobraSuperior, 0.0))
            tx = areaSobraSuperiorUtilizada/areaSobraSuperior;
        else
#ifdef AREANULA_DELTAZERO
        tx = 0.0;
#else
        tx = delta.sobraSuperior;
#endif

#ifndef PODE_ZERAR_DELTAS
        if(!DOUBLE_EQUALS(tx, delta.sobraSuperior) && !DOUBLE_EQUALS(tx,0.0)) {
#else
        if(!DOUBLE_EQUALS(tx, delta.sobraSuperior)) {
#endif
        //    delta.sobraSuperior = (DOUBLE_EQUALS(tx,0.0)?0.01:tx);
            //delta.sobraSuperior = (DOUBLE_EQUALS(tx,0.0)?delta.sobraSuperior:tx);
            //delta.sobraSuperior = tx;
            delta.sobraSuperior = delta.sobraSuperior + (tx - delta.sobraSuperior)*(ALPHA);//pow(valorSigma, iteracao);
            deltaMudou = true;
        }
        else
            deltaMudou = false;
        //cout << tx << ")" << endl;
        return deltaMudou;
    }
    inline bool atualizaDeltaSobras(){
        double tx;
        //cout << id << " " << codigo << " " << codigoObjetoOrigem << " T(";
        //cout << (int)areaSobraSuperior+areaSobraDireita << " " << (int)areaSobraSuperiorUtilizada+areaSobraDireita << " ";
        if(DOUBLE_EQUALS(areaSobraDireita, 0.0) && DOUBLE_EQUALS(areaSobraSuperior, 0.0))
#ifdef AREANULA_DELTAZERO
        tx = 0.0;
#else
        tx = delta.sobras;
#endif
        else
            tx = (areaSobraDireitaUtilizada+areaSobraSuperiorUtilizada)/(areaSobraDireita+areaSobraSuperior);

#ifndef PODE_ZERAR_DELTAS
        if(!DOUBLE_EQUALS(tx, delta.sobras) && !DOUBLE_EQUALS(tx,0.0)) {
#else
        if(!DOUBLE_EQUALS(tx, delta.sobras)) {
#endif
            //delta.sobraSuperior = (DOUBLE_EQUALS(tx,0.0)?0.01:tx);
            //delta.sobras = (DOUBLE_EQUALS(tx,0.0)?delta.sobras:tx);
            //delta.sobras = tx;
            delta.sobras = delta.sobras + (tx - delta.sobras)*(ALPHA);//pow(valorSigma, iteracao);
            deltaMudou = true;
        }
        else
            deltaMudou = false;
        //cout << tx << ")" << endl;
        return deltaMudou;
    }
    inline bool atualizaDeltaObjeto(){
        if(ehSobra) {
            delta.objeto = 1.0;
            return false;
        }
        double tx;
        //cout << id << " " << codigo << " " << codigoObjetoOrigem << " O(";
        //cout << w*h << " " << (int)areaUtilizada << " ";
        if(w*h != 0)
            tx = (areaUtilizada/(w*h));//tx = (1.0-areaUtilizada/(w*h));
        else
#ifdef AREANULA_DELTAZERO
        tx = 0.0;
#else
        tx = delta.objeto;
#endif

        //ASSERT(tx<=1.0);
#ifndef PODE_ZERAR_DELTAS
        if(!DOUBLE_EQUALS(tx, delta.objeto) && !DOUBLE_EQUALS(tx,0.0)) {
#else
        if(!DOUBLE_EQUALS(tx, delta.objeto)) {
#endif
//            deltaMudou = !(DOUBLE_EQUALS(delta.objeto,0.01) && DOUBLE_EQUALS(tx,0.0));
//            //delta.objeto = (DOUBLE_EQUALS(tx,0.0)?0.01:tx);;
            //delta.objeto = (DOUBLE_EQUALS(tx,0.0)?delta.objeto:tx);
            //delta.objeto = tx;
            delta.objeto = delta.objeto + (tx - delta.objeto)*(ALPHA);//pow(valorSigma, iteracao);
            deltaMudou = true;
        }
        else
            deltaMudou = false;
        //cout << tx << ")" << endl;
        return deltaMudou;
    }

};

#endif //OBJETO_H
