import numpy as np;
import math;
from math import sqrt;

from stat_util import *

def unit4() :
    print("--Unit 4 : Distribution Function--")
    f = input("enter function (normal/expo/gamma): ")
    if(f=="normal") :
        ave = float(input("ave"))
        sd = float(input("sd: "))
        comp = input("comp (P(? < X < ?)): ").split()
        if(len(comp)!=2): return
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-ave)/sd
            comp[1] = (comp[1]-ave)/sd
            print(f"P({comp[0]} < Z < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/sd
            print(f"P(Z < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/sd
            print(f"P(Z > {comp[0]})")
    if(f=="expo") :
        ave = float(input("ave (lambda): "))
        mu = 1/ave
        var = 1/(ave**2)
        sd = sqrt(var)
        comp = input("comp (P(? < X < ?)): ").split()
        if(len(comp)!=2): return
        print("-"*40) 
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            lhs = 1-math.exp((-1)*ave*comp[0])
            rhs = 1-math.exp((-1)*ave*comp[1])
            print(f"P({comp[0]} < Z < {comp[1]}) = F({comp[1]}) - F({comp[0]}) = {rhs} - {lhs} = {rhs-lhs}")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            rhs = 1-math.exp((-1)*ave*comp[1])
            print(f"P(Z < {comp[1]}) = F({comp[1]}) = {rhs}")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            rhs = 1-math.exp((-1)*ave*comp[1])
            print(f"P(Z > {comp[1]}) = 1 - F({comp[1]}) = {1-rhs}")
        print(f"ave: {mu}")
        print(f"sd: {sd}")
        print(f"var: {var}")
    if(f=="gamma") :
        ave =  float(input("ave success in 1 unit of time: "))
        r = int(input("number of success: "))
        mu = r/ave
        var = r/(ave**2)
        sd = sqrt(var)
        comp = input("comp (P(? < X < ?)): ").split()
        if(len(comp)!=2): return
        print("-"*40) 
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            lhs = np.sum(np.array([(math.exp((-1)*ave*comp[0])*((ave*comp[0])**k))/math.factorial(k) for k in range(r)]))
            rhs = np.sum(np.array([(math.exp((-1)*ave*comp[1])*((ave*comp[1])**k))/math.factorial(k) for k in range(r)]))
            print(f"P({comp[0]} < Z < {comp[1]}) = F({comp[1]}) - F({comp[0]}) = {lhs} - {rhs} = {lhs-rhs}")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            rhs = np.sum(np.array([(math.exp((-1)*ave*comp[1])*((ave*comp[1])**k))/math.factorial(k) for k in range(r)]))
            print(f"P(Z < {comp[1]}) = 1 - f({comp[1]}) = {1-rhs}")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            rhs = np.sum(np.array([(math.exp((-1)*ave*comp[1])*((ave*comp[1])**k))/math.factorial(k) for k in range(r)]))
            print(f"P(Z > {comp[1]}) = f({comp[1]}) = {rhs}")
        print(f"ave: {mu}")
        print(f"sd: {sd}")
        print(f"var: {var}")






def unit5() :
    print("--Unit 5 : Sampling Distribution--")
    i = int(input("enter variation: "))
    if(i==0) :
        n = int(input("n: "))
        ave = float(input("ave: "))
        sd = float(input("sd: "))
        comp = input("comp (P(? < X < ?)): ").split()
        if(len(comp)!=2): return
        print("-"*40) 
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-ave)/(sd/(sqrt(n)))
            comp[1] = (comp[1]-ave)/(sd/(sqrt(n)))
            print(f"P({comp[0]} < Z < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/(sd/(sqrt(n)))
            print(f"P(Z < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/(sd/(sqrt(n)))
            print(f"P(Z > {comp[0]})")
    if(i==1) :
        n = int(input("n: "))
        ave = float(input("ave (sample): "))
        sd = float(input("sd (sample): "))
        comp = input("comp (P(? < X < ?)): ").split()
        v = n-1
        if(len(comp)!=2): return 
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-ave)/(sd/(sqrt(n)))
            comp[1] = (comp[1]-ave)/(sd/(sqrt(n)))
            print(f"P({comp[0]} < T({v}) < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/(sd/(sqrt(n)))
            print(f"P(T({v}) < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/(sd/(sqrt(n)))
            print(f"P(T({v}) > {comp[0]})")
    if(i==2) :
        n1 = int(input("n 1: "))
        n2 = int(input("n 2: "))
        ave1 = float(input("ave 1: "))
        ave2 = float(input("ave 2: "))
        sd1 = float(input("sd 1: "))
        sd2 = float(input("sd 2: "))
        ave12 = ave1 - ave2
        var12 = ((sd1**2)/n1) + ((sd2**2)/n2)
        comp = input("comp (P(? < X1-X2 < ?)): ").split()
        if(len(comp)!=2): return
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-(ave12))/sqrt(var12)
            comp[1] = (comp[1]-(ave12))/sqrt(var12)
            print(f"P({comp[0]} < Z < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-(ave12))/sqrt(var12)
            print(f"P(Z < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-(ave12))/sqrt(var12)
            print(f"P(Z > {comp[1]})")
    if(i==3) :
        n1 = int(input("n 1: "))
        n2 = int(input("n 2: "))
        ave1 = float(input("ave 1: "))
        ave2 = float(input("ave 2: "))
        sd1 = float(input("sd 1 (sample): "))
        sd2 = float(input("sd 2 (sample): "))
        var_eq = input("Are var equal? (y/n): ")
        varEq = True if var_eq in ("y","Y") else False
        sp2 = (((n1-1)*(sd1**2))+((n2-1)*(sd2**2)))/(n1+n2-2)
        ave12 = ave1 - ave2
        var12 = ((sp2)/n1) + ((sp2)/n2)  if varEq else ((sd1**2)/n1) + ((sd2**2)/n2) 
        v = n1+n2-2 if varEq else ((var12)**2)/(((((sd1**2)/n1)**2)/(n1-1))+((((sd2**2)/n2)**2)/(n2-1)))
        comp = input("comp (P(? < X1-X2 < ?)): ").split()
        if(len(comp)!=2): return 
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-(ave12))/sqrt(var12)
            comp[1] = (comp[1]-(ave12))/sqrt(var12)
            print(f"P({comp[0]} < T({round(v,2)}) < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-(ave12))/sqrt(var12)
            print(f"P(T({round(v,2)}) < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-(ave12))/sqrt(var12)
            print(f"P(T({round(v,2)}) > {comp[1]})")
    if(i==4) :
        n = int(input("n: "))
        is_col = input("collect data? (y/n): ")
        col = True if is_col in ("y","Y") else False
        data =list()
        if(col) :
            print("--Begin Enter Data--")
            for j in range(n) :
                xi = float(input(f"Enter x{j+1}: "))
                yi = float(input(f"Enter y{j+1}: "))
                data.append(xi-yi)
            print("--End Enter Data--")
        dt = np.array(data)
        Dave = np.sum(dt)/n
        if(col) :
            ave = Dave
        else :
            ave = float(input("ave: "))
        if(col) :
            sd = sqrt((np.sum(dt**2)-(n*Dave*Dave))/(n-1))
        else :
            sd = float(input("sd (sample): "))
        v = n-1
        comp = input("comp (P(? < /D < ?)): ").split()
        if(len(comp)!=2): return 
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-ave)/(sd/sqrt(n))
            comp[1] = (comp[1]-ave)/(sd/sqrt(n))
            print(f"P({comp[0]} < T({v}) < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/(sd/sqrt(n))
            print(f"P(T({v}) < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/(sd/sqrt(n))
            print(f"P(T({v}) > {comp[1]})")
    if(i==5) :
        n = int(input("n: "))
        var = float(input("var (population): "))
        v = n-1
        comp = input("comp (P(? < S^2 < ?)): ").split()
        if(len(comp)!=2) : return
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = ((n-1)*comp[0])/var
            comp[1] = ((n-1)*comp[1])/var
            print(f"P({comp[0]} < X2({v}) < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = ((n-1)*comp[1])/var
            print(f"P(X2({v}) < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = ((n-1)*comp[1])/var
            print(f"P(X2({v}) > {comp[1]})")
    if(i==6) :
        n1 = int(input("n 1: "))
        n2 = int(input("n 2: "))
        varp1 = float(input("var 1 (population): "))
        varp2 = float(input("var 2 (population): "))
        df1 = n1-1
        df2 = n2-1
        comp = input("comp (? < (s1^2)/(s2^2) < ?): ").split()
        if len(comp)!=2 : return
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = comp[0]*(varp2/varp1)
            comp[1] = comp[1]*(varp2/varp1)
            print(f"P({comp[0]} < F({df1},{df2}) < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = comp[1]*(varp2/varp1)
            print(f"P(F({df1},{df2}) < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = comp[1]*(varp2/varp1)
            print(f"P(F({df1},{df2}) > {comp[1]})")
    if(i==7) :
        n = int(input("n: "))
        ave = float(input("ave probability to happen (population): "))
        comp = input("comp (? < P_hat < ?): ").split()
        if len(comp)!=2 : return
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-ave)/sqrt((ave*(1-ave))/n)
            comp[1] = (comp[1]-ave)/sqrt((ave*(1-ave))/n)
            print(f"P({comp[0]} < Z < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/sqrt((ave*(1-ave))/n)
            print(f"P(Z < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-ave)/sqrt((ave*(1-ave))/n)
            print(f"P(Z > {comp[1]})")
    if(i==8) :
        n1 = int(input("n 1: "))
        n2 = int(input("n 2: "))
        p1 = float(input("ave to happen for P1 (Population): "))
        p2 = float(input("ave to happen for P2 (Population): "))
        comp = input("comp (? < P_hat1 - P_hat2 < ?): ").split()
        if len(comp)!=2 : return
        print("-"*40)
        if(isnumber(comp[0])) :
            comp[0] = float(comp[0])
            comp[1] = float(comp[1])
            comp[0] = (comp[0]-(p1-p2))/sqrt(((p1*(1-p1))/n1)+((p2*(1-p2))/n2))
            comp[1] = (comp[1]-(p1-p2))/sqrt(((p1*(1-p1))/n1)+((p2*(1-p2))/n2))
            print(f"P({comp[0]} < Z < {comp[1]})")
        elif comp[0] == "<" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-(p1-p2))/sqrt(((p1*(1-p1))/n1)+((p2*(1-p2))/n2))
            print(f"P(Z < {comp[1]})")
        elif comp[0] == ">" :
            comp[1] = float(comp[1])
            comp[1] = (comp[1]-(p1-p2))/sqrt(((p1*(1-p1))/n1)+((p2*(1-p2))/n2))
            print(f"P(Z > {comp[1]})")




def unit7() :
    print("--Unit 7 : Test of Hypothesis--")
    i = int(input("enter variation: "))
    if(i==2) :
        yes = ("y","Y")
        comp = input("H1 (ave1-ave2 ? d): ").split()
        d = float(comp[1])
        alpha = float(input("alpha: "))
        if(len(comp)!=2) : return
        n1 = int(input("n1: "))
        n2 = int(input("n2: "))
        ave1 = float(input("ave 1: "))
        ave2 = float(input("ave 2: "))
        var1 = float(input("var 1: "))
        var2 = float(input("var 2: "))
        n_large = input("consider n large? (y/n): ")
        n_large = True if n_large in yes else False
        zt_alpha = 0
        v=-1
        if(n_large) :
            zt = ((ave1-ave2)-d)/sqrt((var1/n1)+(var2/n2))
            dash()
            print(f"Z = {zt}")
        else :
            know_sd = input("do you know sd? (y/n): ")
            know_sd = True if know_sd in yes else False
            if(know_sd) :
                zt = ((ave1-ave2)-d)/sqrt((var1/n1)+(var2/n2))
                dash()
                print(f"Z = {zt}")
            else:
                sd_eq = input("are sd equal? (y/n): ")
                sd_eq = True if sd_eq in yes else False
                if(sd_eq) :
                    v = n1+n2-2
                    sp2 = (((n1-1)*var1)+((n2-1)*var2))/(n1+n2-2)
                    zt = ((ave1-ave2)-d)/sqrt((sp2/n1)+(sp2/n2))
                    dash()
                    print(f"T = {zt}")
                else : 
                    v = (((var1/n1)+(var2/n2))**2)/((((var1/n1)**2)/(n1-1))+(((var2/n2)**2)/(n2-1)))
                    v = round(v,2)
                    zt = ((ave1-ave2)-d)/sqrt((var1/n1)+(var2/n2))
                    dash()
                    print(f"T = {zt}")
        if(comp[0]=="!=") :
            if v==-1:
                zt_alpha = float(input(f"z_{alpha/2}: "))
            else :
                zt_alpha = float(input(f"t_{alpha/2}_{v}: "))
            H1 = math.fabs(zt) > zt_alpha
            H0 = not H1
        elif(comp[0]==">") :
            if v==-1:
                zt_alpha = float(input(f"z_{alpha}: "))
            else :
                zt_alpha = float(input(f"t_{alpha}_{v}: "))
            H1 = zt > zt_alpha
            H0 = not H1
        elif(comp[0]=="<") :
            if v==-1:
                zt_alpha = float(input(f"z_{alpha}: "))
            else :
                zt_alpha = float(input(f"t_{alpha}_{v}: "))
            H1 = zt < (-1)*zt_alpha
            H0 = not H1
        print(f"H0: {H0}")
        print(f"H1: {H1}")
    if(i==3) :
        yes = ("y","Y")
        comp = input("H1 (ave1-ave2 ? d): ").split()
        d = float(comp[1])
        alpha = float(input("alpha: "))
        n = int(input("n: "))
        col = input("collect data? (y/n): ")
        col = True if col in yes else False
        data =list()
        if(col) :
            print("--Begin Enter Data--")
            for j in range(n) :
                xi = float(input(f"Enter x{j+1}: "))
                yi = float(input(f"Enter y{j+1}: "))
                data.append(xi-yi)
            print("--End Enter Data--")
            dt = np.array(data)
            ave_d = np.sum(dt)/n
            sd_d = sqrt((np.sum(dt**2)-(n*(ave_d**2)))/(n-1))
        else :
            ave_d = float(input("ave of d: "))
            sd_d = float(input("sd of d: "))
        T = (ave_d-d)/(sd_d/sqrt(n))
        v = n-1
        dash()
        print(f"ave of d: {ave_d}")
        print(f"sd of d: {sd_d}")
        print(f"T: {T}")
        if(comp[0]=="!=") :
            t_alpha = float(input(f"t_{alpha/2}_{v}: "))
            H1 = math.fabs(T) > t_alpha
            H0 = not H1
        elif(comp[0]==">") :
            t_alpha = float(input(f"t_{alpha}_{v}: "))
            H1 = T > t_alpha
            H0 = not H1
        elif(comp[0]=="<") :
            t_alpha = float(input(f"t_{alpha}_{v}: "))
            H1 = T < (-1)*t_alpha
            H0 = not H1
        print(f"H0: {H0}")
        print(f"H1: {H1}")     
    if(i==6) :
        yes = ("y","Y")
        comp = input("H1 (p ? p0): ").split()
        if(len(comp)!=2) : return 
        p0 = float(comp[1])
        alpha = input("alpha: ")
        n = int(input("n: "))
        p_hat = float(input("p_Hat (ratio that we got from sampling): "))
        Z = (p_hat-p0)/sqrt((p0*(1-p0))/n)
        dash()
        print(f"Z: {Z}")
        if(comp[0]=="!=") :
            z_alpha = float(input(f"z_{alpha/2}: "))
            H1 = math.fabs(Z) > z_alpha
            H0 = not H1
        elif(comp[0]==">") :
            z_alpha = float(input(f"z_{alpha}: "))
            H1 = Z > z_alpha
            H0 = not H1
        elif(comp[0]=="<") :
            z_alpha = float(input(f"z_{alpha}: "))
            H1 = Z < (-1)*z_alpha
            H0 = not H1
        print(f"H0: {H0}")
        print(f"H1: {H1}")     
    if(i==7):
        yes = ("y","Y")
        comp = input("H1 (p1-p2 ? d): ").split()
        if(len(comp)!=2) : return
        d = float(comp[1])
        alpha = float(input("alpha: "))
        n1 = int(input("n1: "))
        n2 = int(input("n2: "))
        a1 = int(input("a1: "))
        a2 = int(input("a2: "))
        p_hat1 = a1/n1
        p_hat2 = a2/n2
        if d==0 :
            pbar = (a1+a2)/(n1+n2)
            Z = ((p_hat1-p_hat2)-(d))/sqrt(pbar*(1-pbar)*((1/n1)+(1/n2)))
        else :
            Z = ((p_hat1-p_hat2)-(d))/sqrt(((p_hat1*(1-p_hat1))/n1)+((p_hat2*(1-p_hat2))/n2))
        dash()
        print(f"Z: {Z}")
        if(comp[0]=="!=") :
            z_alpha = float(input(f"z_{alpha/2}: "))
            H1 = math.fabs(Z) > z_alpha
            H0 = not H1
        elif(comp[0]==">") :
            z_alpha = float(input(f"z_{alpha}: "))
            H1 = Z > z_alpha
            H0 = not H1
        elif(comp[0]=="<") :
            z_alpha = float(input(f"z_{alpha}: "))
            H1 = Z < (-1)*z_alpha
            H0 = not H1
        print(f"H0: {H0}")
        print(f"H1: {H1}")     





def unit8() :
    print("-- Unit 8 : One-Way Analysis of Variance --")
    k = int(input("k: "))
    data = list()
    raw = list()
    N = 0
    for i in range(k) :
        print(f"--Begin Enter Data {i+1}--")
        dt = collect()
        N += len(dt)
        raw.append(dt)
        data.append(analyze(dt))
    print("--End Enter Data--")
    F_alpha = float(input("F_alpha_d.f-Between-G_d.f-Within-G : "))
    TT = np.sum(np.array([dt["sum"] for dt in data]))
    SST = np.sum(np.array([dt["sum2"] for dt in data])) - ((TT**2)/N)
    SSB = np.sum(np.array([(dt["sum"]**2)/dt["n"] for dt in data])) - ((TT**2)/N)
    SSW = SST-SSB
    dof1 = k-1
    dof2 = N-k
    MSB = SSB/(k-1)
    MSE = SSW/(N-k)
    F_Test = MSB/MSE
    H0 = not F_Test > F_alpha
    H1 = not H0
    print("-"*40)
    for i in range(k) :
        print(f"data {i+1}: {raw[i]}")
        print(f"analysis: {data[i]}")
    print(f"k: {k}")
    print(f"N: {N}")
    print(f"TT: {TT}")
    print(f"SST: {SST}")
    print(f"SSB: {SSB}")
    print(f"SSW: {SSW}")
    print(f"dof between-groups: {dof1}")
    print(f"dof within-groups: {dof2}")
    print(f"MSB: {MSB}")
    print(f"MSE: {MSE}")
    print(f"F_Test: {F_Test}")
    print(f"F_alpha_d.f-Between-G_d.f-Within-G: {F_alpha}")
    print(f"H0: {H0}")
    print(f"H1: {H1}")

def unit9() :
    print("-- Unit 9 : Linear Regression and Correlation --")
    n = int(input("n : "))
    print("--Begin Enter Data x-")
    datax = collect(n)
    print("--Begin Enter Data y--")
    datay = collect(n)
    print("--End Enter Data--")
    dtx = np.array(datax)
    dty = np.array(datay)
    resultx = analyze(datax)
    resulty = analyze(datay)
    Sxy = np.sum(dtx*dty)-((resultx["sum"]*resulty["sum"])/n)
    Sxx = resultx["sum2"]-((resultx["sum"]**2)/n)
    Syy = resulty["sum2"]-((resulty["sum"]**2)/n)
    r = Sxy/math.sqrt(Sxx*Syy)
    T = r*math.sqrt((n-2)/(1-(r**2)))
    b1 = Sxy/Sxx
    b0 = resulty["ave"]-(b1*resultx["ave"])
    dtyhat = (dtx*b1)+b0
    e2 = (dty-dtyhat)**2
    SEE = math.sqrt(np.sum(e2)/(n-2))
    R2 = (b1*Sxy)/Syy
    c = input("Check Correlation? (y/n): ")
    if(c == "y" or c == "Y") :
        t_alpha = float(input("t_alpha/2_n-2: "))
        h0 = not (math.fabs(T) > t_alpha)  
        h1 = not h0
    print("-"*40)
    print(f"data x: {datax}")
    print(resultx)
    print(f"data y: {datay}")
    print(resulty)
    print(f"Sxy: {Sxy}")
    print(f"Sxx: {Sxx}")
    print(f"Syy: {Syy}")
    print(f"r: {r}")
    print(f"T: {T}")
    if(c == "y" or c == "Y") :
        print(f"t_alpha/2_n-2: {t_alpha}")
        print(f"H0: {h0}")
        print(f"H1: {h1}")
    print(f"b0: {b0}")
    print(f"b1: {b1}")
    print(f"SEE: {SEE}")
    print(f"R2: {R2}")