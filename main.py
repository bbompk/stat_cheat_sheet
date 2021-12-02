import os
import sys;

from stat_util import *
from stat_formula import *

def renovation() :
    print()
    print("*"*40)
    command  = input("Choose Command: ")
    if(command=="a") :
        print("--Analyze Data--")
        data = collect()
        result = analyze(data)
        skw = skew(data, result)
        kurt = kurtosis(data,result)
        print("-"*40)
        print("data:",data)
        print(f"analysis: {result}")
        print(f"skewness: {skw}")
        print(f"kurtosis: {kurt}")
    if(command=="b") :
        print("--Analyze Data X-Y--")
        print("--Begin Enter data x--")
        datax = collect()
        print("--Begin Enter data y--")
        datay= collect()
        print("--End Entering data--")
        resultxy = analyze_x_minus_y(datax, datay)
        if len(resultxy) == 0 :
            print("n is not equal")
            return True
        resultx = analyze(datax)
        skwx = skew(datax, resultx)
        kurtx = kurtosis(datax,resultx)
        resulty = analyze(datay)
        skwy = skew(datay, resulty)
        kurty = kurtosis(datay, resulty)
        dataxy = resultxy["data"]
        skwxy = skew(dataxy, resultxy)
        kurtxy = kurtosis(dataxy, resultxy)
        print("-"*40)
        print("data x:",datax)
        print(f"analysis x: {resultx}")
        print(f"skewness x: {skwx}")
        print(f"kurtosis x: {kurtx}")
        print("data y:",datay)
        print(f"analysis y: {resulty}")
        print(f"skewness y: {skwy}")
        print(f"kurtosis y: {kurty}")
        print("data xy:",dataxy)
        print(f"analysis xy: {resultxy}")
        print(f"skewness xy: {skwxy}")
        print(f"kurtosis xy: {kurtxy}")
    if(command=="f") :
        print("--Formula--")
        unit = input("Choose Unit: ")
        if(unit=="4") :
            unit4()
        if(unit=="5") :
            unit5()
        if(unit=="7") :
            unit7()
        if(unit=="8") :
            unit8()
        if(unit=="9") :
            unit9()
        
    if(command == "exit" or command == "q") : return False


    return True


if(__name__ == "__main__") :
    while(True) :
        if(not renovation()) : break