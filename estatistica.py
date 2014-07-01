# python module
# statistical quality indicators for air quality model performance evaluation
# carlagama

from math import *
from numpy import *
from decimal import *

__all__ = ['agreement', 'nmse', 'fac2', 'fractional_bias', 'boot_bias', 'anb',
           'rms', 'mg', 'vg', 'nsd', 'correlation']


def agreement(observed, predicted):
    aux = average(observed)
    average_observed = [aux] * len(predicted)
    numerador = sum(map(_aux_numerador, predicted, observed))
    denominador = sum(map(_aux_denominador, predicted, observed,
                          average_observed))
    return 1 - numerador / denominador

def _aux_numerador(x, y):
    return (x - y)**2

def _aux_denominador(x, y, average):
    modulo = (abs(x - average)) + (abs(y - average))
    return modulo**2

def nmse(observed, predicted):
    result = average(map(_nmse1, observed, predicted)) / (average (observed) *
                                                          average(predicted))
    return result

def _nmse1(observed, predicted):
    return (observed - predicted)**2

def fac2(observed, predicted):
    a = map(_fa, observed, predicted)
    result = Decimal(sum(a)) / len(observed)
    return result

def _fa(observed, predicted):
    if (predicted / observed <= 2.0 and predicted / observed >= 0.5):
        valor = 1
    else:
        valor = 0
    return valor

def fractional_bias(observed, predicted):   
    mean_observed = sum(observed) / len(observed)
    mean_predicted = sum(predicted) / len(predicted)
    return (mean_observed - mean_predicted) / (0.5 * (mean_observed +
                                                      mean_predicted))

def boot_bias(observed, predicted):   
    mean_observed = sum(observed) / len(observed)
    mean_predicted = sum(predicted) / len(predicted)
    return (mean_observed - mean_predicted)

def anb(observed, predicted):
    result = average(map(_anb_aux, observed, predicted))
    return result

def _anb_aux(observed, predicted):
    return (abs(observed - predicted))/observed
    
def rms(observed, predicted):
    N = len(observed)
    soma = sum(map(_nmse1, observed, predicted))
    return sqrt(soma/N)

def mg(observed, predicted):
    a = average(map(_ln, observed))
    b = average(map(_ln, predicted))
    return exp(a-b)

def _ln(x):
    return log(x)

def vg(observed, predicted):
    a = average(map(_ln, observed))
    b = average(map(_ln, predicted)) 
    return exp((a - b) ** 2)

def nsd(observed, predicted):
    return std(predicted) / std(observed)

def correlation(observed, predicted):
    aux_o = average(observed)
    aux_p = average(predicted)
    avr_observed = [aux_o] * len(observed)
    avr_predicted = [aux_p] * len(predicted)
    a = sum(map(_aux_r_num, observed, predicted, avr_observed, avr_predicted))
    b = std(observed)*std(predicted)*(len(observed)-1)
    return a / b

def _aux_r_num(observed, predicted, avr_observed, avr_predicted):
    return (observed - avr_observed) * (predicted - avr_predicted)

