//
// Created by Oberlan Romão on 05/09/15.
//

#ifndef DELTA_H
#define DELTA_H
#include "Util.h"

class Delta{
public:
    double objeto; //Delta do objeto, se o objeto é um dado de entrada, então obj = 1.0
    double sobraDireita; //Delta da sobra do lado direito do objeto
    double sobraSuperior; //Delta da sobra superior do objeto
    double sobras; //Delta das sobras do objeto
    Delta(double obj=valorInicialDelta, double r=valorInicialDelta, double t=valorInicialDelta, double s=valorInicialDelta):
            sobraDireita(r), sobraSuperior(t), objeto(obj), sobras(s){}
};
#endif //DELTA_H
