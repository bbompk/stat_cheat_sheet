import numpy as np;
import math;
from math import sqrt;

def collect(a=5000) :
    data = list()
    i = 0
    while(True) :
        if(i >= a): break
        inp = input("Enter data: ")
        if(not inp.isnumeric()): break
        data.append(int(inp))
        i+=1
    return data

def skew(data, result) :
    ave = result["ave"]
    sd = result["sd"]
    sd1 = result["sd1"]
    N = result["n"]
    dt = np.array(data)
    sk = np.sum(((dt-ave)**3)/((sd**3)*N))
    sk1 = ((N**2)/((N-1)*(N-2)))*np.sum(((dt-ave)**3)/((sd1**3)*N))
    skw = {
        "pop" : sk,
        "samp" : sk1
    }
    return skw
    
def kurtosis(data, result) :
    ave = result["ave"]
    sd = result["sd"]
    sd1 = result["sd1"]
    N = result["n"]
    dt = np.array(data)
    kt = np.sum(((dt-ave)**4)/((sd**4)*N))-3
    kt1 = ((((N**2)*(N+1))/((N-1)*(N-2)*(N-3)))*np.sum(((dt-ave)**4)/((sd1**4)*N))) - ((3*((N-1)**2))/((N-2)*(N-3)))
    kurt = {
        "pop" : kt,
        "samp" : kt1
    }
    return kurt
    
def analyze(data) :
    n = len(data)
    if n==0 : return dict()
    dt = np.array(data)
    summ = np.sum(dt)
    sum2 = np.sum(dt*dt)
    ave = summ/n

    dts = (dt-ave)**2
    sums = np.sum(dts)

    var = sums/n
    var1 = sums/(n-1)
    sd = math.sqrt(var)
    sd1 = math.sqrt(var1)
    result = {
        "n" : n,
        "ave" : ave,
        "sum" : summ,
        "sum2" : sum2,
        "sd" : sd,
        "sd1" : sd1,
        "var" : var,
        "var1" : var1
    }
    return result

def analyze_x_minus_y(datax, datay) :
    n1 = len(datax)
    n2 = len(datay)
    if n1!=n2 : return dict()
    n = n1
    dtx = np.array(datax)
    dty = np.array(datay)
    dt = dtx-dty
    summ = np.sum(dt)
    sum2 = np.sum(dt*dt)
    ave = summ/n

    dts = (dt-ave)**2
    sums = np.sum(dts)

    var = sums/n
    var1 = sums/(n-1)
    sd = math.sqrt(var)
    sd1 = math.sqrt(var1)
    result = {
        "data" : list(dt),
        "n" : n,
        "ave" : ave,
        "sum" : summ,
        "sum2" : sum2,
        "sd" : sd,
        "sd1" : sd1,
        "var" : var,
        "var1" : var1
    }
    return result

def isnumber(string) :
    string = string.strip("-")
    if string.isnumeric() :
        return True
    if string.count(".") == 1:
        return True
    return False

def dash() :
    print("-"*40)
    