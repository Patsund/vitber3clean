import numpy as np
from matplotlib import pyplot as plt
import time

alpha = 5e-5
m = 1e-2

def Vwater(t,X):
    T=24*3600
    factor = 2 * np.pi / T
    newX = [factor*(-X[1]), factor*X[0]]
    return newX

def eulerForEq1(X,t,Vwater,Xdot,h):
    x, y = X[0], X[1]
    VwaterVector = Vwater(t,X)
    x += h*Xdot[0]
    y += h*Xdot[1]
    newX = [x,y]
    dx, dy = Xdot[0],Xdot[1]
    dx += h*(alpha/m)*(VwaterVector[0]-dx)
    dy += h*(alpha/m)*(VwaterVector[1]-dy)
    return newX

def eulerForEq2(X,h):
    x, y=X[0], X[1]
    newX=[x,y]
    T=24*60*60
    factor=2*np.pi/T
    newX=[newX[0]+h*factor*(-y),newX[1]+h*factor*(x)]
    return newX

def fForEq2(X,t):
    dX=np.array(Vwater(t,X))
    return dX

def fForEq1(X,t):
    Vw = Vwater(t, [X[0], X[1]])
    dx, dy = X[2], X[3]
    dvx = alpha/m*(Vw[0]-X[2])
    dvy = alpha/m*(Vw[1]-X[3])
    return [dx,dy,dvx,dvy]

def rk2(X,f,h,t): #X er en array som kan inneholde enten x,y eller x,y,vx,vy avhengin av om hhv eq 2 eller eq 1 skal brukes
    K1 = np.array(f(X, t))
    K2 = np.array(f(X+h*K1, t+h))
    X_= X+h/2*(K1+K2)
    return X_

def errorAndTrajectoryForRK2(X,h,figname):
    h0=h
    timeNow, timeFinal = 0, 48*3600
    n=int(timeFinal/h)
    newX=np.array(X)
    xList,yList=[],[]
    xList.append(newX[0])
    yList.append(newX[1])
    for i in range(n+1):
        h=min(h,timeFinal-timeNow)
        newX=rk2(newX,fForEq2,h,timeNow)
        xList.append(newX[0])
        yList.append(newX[1])
        timeNow+=h
    X=np.array(X)
    newX=np.array(newX)
    vecError=X-newX
    if figname=="dummy":
        X=[]
        plt.figure("1b h="+str(h0)+" Trap")
        plt.title("Trajectory for particle using trapezoid, h = "+str(h0))
        plt.plot(xList,yList)
    elif figname=="errorOnly":
        pass
    else:
        plt.figure(figname)
        plt.title("Tracectory for particle using trapezoid, h = "+str(h0))
        plt.plot(xList,yList)
        plt.savefig("Oppgave1pdfer\opg"+figname)
        plt.close()
    error=np.sqrt(vecError[0]**2+vecError[1]**2)
    return error

def errorAndTrajectoryForEuler(X,h,figname):
    h0=h
    x,y=X[0],X[1]
    xList,yList=[],[]
    newX=[x,y] #oppretter en ny matrise for å unngå å endre på X. Må skje siden vi endrer på newX ila neste for-løkke
    xList.append(newX[0])
    yList.append(newX[1])
    timeFinal=2*24*60*60
    timeNow=0
    for i in range(int(2*24*60*60/h)+1):
        h=min(h,timeFinal-timeNow)
        timeNow+=h
        newX=eulerForEq2(newX,h)
        xList.append(newX[0])
        yList.append(newX[1])
    X=np.array(X)#gjør om til numpyarrays for å kunne gjøre vektoraritmetikk. Kunne i utgangspunktet tatt inn nparrays,
    newX=np.array(newX)#men funksjonen er nå kompatibel med både vanlige arrays og nparrays.
    vecError=X-newX #finner vektoren fra X til newX
    if figname=="dummy":
        plt.figure("1a h="+str(h0)+" Eul")
        plt.title("Trajectory for particle using Euler, h = "+str(h0))
        plt.plot(xList,yList)
    elif figname=="errorOnly":
        pass
    else:
        plt.figure(figname)
        plt.title("Trajectory for particle using Euler, h = "+str(h0))
        plt.plot(xList,yList)
        plt.savefig("Oppgave1pdfer\opg"+figname)
        plt.close()
    error=np.sqrt(vecError[0]**2+vecError[1]**2)#finner avstanden mellom startpunktet og sluttpunktet
    return error

def errorFor1c(h,analyticEndpoint):
    L = 1.0E+02
    numberOfTimeSteps = int(2 * 24 * 60 * 60 / h) + 1
    X = np.array([L, 0, 0, 0])
    timeNow=0
    timeFinal=48*3600
    for j in range(numberOfTimeSteps):
        h=min(h,timeFinal-timeNow)
        X=rk2(X,fForEq1,h,timeNow)
        timeNow+=h
    error = np.linalg.norm(np.array([X[0], X[1]]) - analyticEndpoint)
    return error

def analyticSolution(L,alpha,m):
    k = alpha/m
    w = 2*np.pi/(24*60*60)*k
    A = np.array([
    [ 0, 0, 1, 0],
    [ 0, 0, 0, 1],
    [ 0,-w,-k, 0],
    [ w, 0, 0,-k]
    ])
    lams , V = np.linalg.eig(A)
    X0 = np.array([L,0.,0.,0.])
    C = np.linalg.solve(V,X0)
    t = 2*24*3600
    X = V.dot(C*np.exp(lams*t)) # pun
    X = [X[0].real, X[1].real]
    #print("Analytic position after 48 hours = ", X[0].real, X[1].real)
    return X

def task1a(savefig=False):
    print("**Oppgave 1a**")
    X=[100,0]
    h0=60*60 #sekund
    thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    while thiserror>10:#finner høyeste tidssteg som gir feil mindre enn 10 meter
        h0-=1
        thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    timestepEuler=h0
    print("Det minste tidssteget som kreves for å ikke få en global feil større enn 10 m med Euler er h =",timestepEuler)
    if savefig:
        X=[100,0]
        h=100
        errorAndTrajectoryForEuler(X,h,"1ahlik100.pdf")
        times=np.logspace(1,4,10)
        errorArr=[]
        for h in times:
            X=[100,0]
            errorArr.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1aerror")
        plt.title("Error for Euler")
        plt.loglog(times,errorArr,"ro")
        plt.axvline(x=timestepEuler)
        plt.savefig("Oppgave1pdfer\opg1aerror.pdf")
    else:
        X=[100,0]
        h=100
        errorAndTrajectoryForEuler(X,h,"dummy")
        times=np.logspace(1,4,10)
        errorArr=[]
        for h in times:
            X=[100,0]
            errorArr.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1aerror")
        plt.title("Error for Euler")
        plt.loglog(times,errorArr,"ro")
        plt.axvline(x=timestepEuler)
        #plt.show()

def task1b(savefig=False):
    print("**Oppgave 1b**")
    X=[100,0]
    h0=60*60 #sekund
    thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    while thiserror>10:#finner høyeste tidssteg som gir feil mindre enn 10 meter
        h0-=1
        thiserror = errorAndTrajectoryForEuler(X,h0,"errorOnly")
    timestepEuler=h0
    h0=60*60 #sekund
    thiserror = errorAndTrajectoryForRK2(X,h0,"errorOnly")
    while thiserror>10:
        h0-=1
        thiserror = errorAndTrajectoryForRK2(X,h0,"errorOnly")
    timestepTrap=h0
    print("Det minste tidssteget som kreves for å ikke få en global feil større enn 10 m med trapesmetoden er h =",timestepTrap)
    if savefig:
        X=[100,0]
        h=100
        errorAndTrajectoryForRK2(X,h,"1ahlik100.pdf")
        times=np.logspace(1,4,10)
        errorArr1=[]
        errorArr2=[]
        for h in times:
            X=[100,0]
            errorArr1.append(errorAndTrajectoryForRK2(X,h,"errorOnly"))
            errorArr2.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1berror")
        plt.title("Error for Euler and trapezoid with biggest error smaller than 10 marked with a line")
        plt.loglog(times,errorArr1,"bo")
        plt.loglog(times,errorArr2,"ro")
        plt.axvline(x=timestepEuler, color="r")
        plt.axvline(x=timestepTrap, color="b")
        #plt.legend(bbox_to_anchor=(0., 1.1, 1., .052), loc="best", ncol=2, mode="expand", borderaxespad=0.)
        plt.savefig("Oppgave1pdfer\opg1berror.pdf")
    else:
        X=[100,0]
        h=100
        errorAndTrajectoryForRK2(X,h,"dummy")
        times=np.logspace(1,4,10)
        errorArr1=[]
        errorArr2=[]
        for h in times:
            X=[100,0]
            errorArr1.append(errorAndTrajectoryForRK2(X,h,"errorOnly"))
            errorArr2.append(errorAndTrajectoryForEuler(X,h,"errorOnly"))
        plt.figure("1berror")
        plt.title("Error for Euler and trapezoid with biggest error smaller than 10 marked with a line")
        plt.loglog(times,errorArr1,"bo")
        plt.loglog(times,errorArr2,"ro")
        plt.axvline(x=timestepEuler, color="r")
        plt.axvline(x=timestepTrap, color="b")
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #plt.show()
    X=(100,0)
    n=10
    endTimeEuler,endTimeRK2=0,0
    for i in range(n):
        startTimeEuler=time.time()
        errorEuler=errorAndTrajectoryForEuler(X,timestepEuler,"errorOnly")
        endTimeEuler+=time.time()-startTimeEuler
        startTimeRK2=time.time()
        errorETM=errorAndTrajectoryForRK2(X,timestepTrap,"errorOnly")
        endTimeRK2+=time.time()-startTimeRK2
    endTimeEuler/=n
    endTimeRK2/=n
    print("Euler bruker",endTimeEuler,"sekunder med error < 10 m. Trapesmetoden bruker",endTimeRK2,"sekunder med error < 10 m")

def task1c(savefig=False):
    print("**Oppgave 1c**")
    if savefig:
        L = 1.0E+02
        differentTimesteps = 10
        analyticEndpoint = np.array(analyticSolution(L, alpha, m))
        timestepArray = np.logspace(2, 3, num=differentTimesteps)
        errorArray = np.zeros(differentTimesteps)
        for i in range(differentTimesteps):
            h = timestepArray[i]
            numberOfTimeSteps = int(2 * 24 * 60 * 60 / h) + 1
            coordinateArray = np.zeros((numberOfTimeSteps,2))
            X = np.array([L, 0, 0, 0])
            timeNow=0
            timeFinal=48*3600
            for j in range(numberOfTimeSteps):
                h=min(h,timeFinal-timeNow)
                X=rk2(X,fForEq1,h,timeNow)
                timeNow+=h
                coordinateArray[j] = np.array([X[0], X[1]])
            xValueArray = [k[0] for k in coordinateArray]
            yValueArray = [k[1] for k in coordinateArray]
            errorArray[i] = np.linalg.norm(np.array([X[0], X[1]]) - analyticEndpoint)
            if i == 3:
                plt.figure("1c"+str(i))
                plt.title("tidssteg " + str(h))
                plt.plot(xValueArray, yValueArray, 'ro', markersize=1)
                plt.plot([analyticEndpoint[0]], [analyticEndpoint[1]], 'ko', markersize=4,label="analytisk ende")
                plt.legend(loc="best")
        plt.figure()
        plt.loglog(timestepArray,errorArray)
        plt.xlabel("timestep / sekunder")
        plt.ylabel("error / meter")
        plt.title("Avvik fra analytisk løsning")
        plt.savefig("Oppgave1pdfer\opg1cplot.pdf")
        h0=350#sekund
        thiserror = errorFor1c(h0,analyticEndpoint)
        while thiserror>10:
            h0-=0.5
            thiserror = errorFor1c(h0,analyticEndpoint)
        timestep1c=h0
        print("Det minste tidssteget som kreves for å ikke få en global feil større enn 10 m med trapesmetoden er h =",timestep1c)
    else:
        L = 1.0E+02
        differentTimesteps = 10
        analyticEndpoint = np.array(analyticSolution(L, alpha, m))
        timestepArray = np.logspace(2, 3, num=differentTimesteps)
        errorArray = np.zeros(differentTimesteps)
        for i in range(differentTimesteps):
            h = timestepArray[i]
            numberOfTimeSteps = int(2 * 24 * 60 * 60 / h) + 1
            coordinateArray = np.zeros((numberOfTimeSteps,2))
            X = np.array([L, 0, 0, 0])
            timeNow=0
            timeFinal=48*3600
            for j in range(numberOfTimeSteps):
                h=min(h,timeFinal-timeNow)
                X=rk2(X,fForEq1,h,timeNow)
                timeNow+=h
                coordinateArray[j] = np.array([X[0], X[1]])
            xValueArray = [k[0] for k in coordinateArray]
            yValueArray = [k[1] for k in coordinateArray]
            errorArray[i] = np.linalg.norm(np.array([X[0], X[1]]) - analyticEndpoint)
            if i == 3:
                plt.figure("1c"+str(i))
                plt.title("tidssteg " + str(h))
                plt.plot(xValueArray, yValueArray, 'ro', markersize=1)
                plt.plot([analyticEndpoint[0]], [analyticEndpoint[1]], 'ko', markersize=4, label="analytisk ende")
                plt.legend(loc="best")
        plt.figure("1c avvik")
        plt.loglog(timestepArray,errorArray)
        plt.xlabel("timestep / sekunder")
        plt.ylabel("error / meter")
        plt.title("Avvik fra analytisk løsning")
        h0=420#sekund
        thiserror = errorFor1c(h0,analyticEndpoint)
        while thiserror>10:
            h0-=0.5
            thiserror = errorFor1c(h0,analyticEndpoint)
        timestep1c=h0
        print("Det minste tidssteget som kreves for å ikke få en global feil større enn 10 m med trapesmetoden er h =",timestep1c)
        #plt.show()

def task1d(savefig=False):
    print("**Oppgave 1d**")
    if savefig:
        L = 1.0E+02
        timeNow = 0
        totalTime = 48* 3600
        deviationLimit = 0.1
        # X = [x, y, dx, dy]
        X = np.array([L, 0, 0, 0])
        #Vi ønsker å logge både tidssteg, avvik for steg som godkjennes samt coordinater
        timeStepArray = []
        deviationArray = []
        coordinateArray = []
        analyticEndpoint = analyticSolution(L, alpha, m)
        while (timeNow < totalTime):
            if timeNow==0:
                h=1000
                Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
                Xtrap = rk2(X,fForEq1,h,timeNow)
                thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
                while thisDeviation > deviationLimit:
                    #Først finjusteres h slik at det første steget akkurat har kommer innenfor nøyaktighetskravet
                    h-=1
                    Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
                    Xtrap = rk2(X, fForEq1, h, timeNow)
                    thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
            h = min(h, totalTime - timeNow)
            Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
            Xtrap = rk2(X, fForEq1, h, timeNow)
            thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
            while (thisDeviation > deviationLimit):
                h = h/2
                Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
                Xtrap = rk2(X,fForEq1,h,timeNow)
                thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
            timeStepArray.append([timeNow, h])
            deviationArray.append(thisDeviation)
            timeNow += h
            if (thisDeviation < deviationLimit/10):
                h *= 2
            coordinateArray.append(Xtrap)
            X = Xtrap
        xValueArray = [k[0] for k in coordinateArray]
        yValueArray = [k[1] for k in coordinateArray]
        timeStepArray_xValues = [k[0] for k in timeStepArray]
        timeStepArray_yValues = [k[1] for k in timeStepArray]
        plt.figure("1dTidssteg")
        plt.title("Tidssteg")
        plt.plot(timeStepArray_xValues, timeStepArray_yValues, 'go', markersize=1)
        plt.xlabel("tid / sekunder")
        plt.ylabel("tidssteg / sekunder")
        plt.savefig("Oppgave1pdfer\opg1dtidssteg.pdf")
        plt.figure("1dAvvik")
        plt.title("Avvik mellom Euler og Trapes")
        plt.plot(timeStepArray_xValues, deviationArray, 'ko', markersize=1)
        plt.xlabel("tid / sekunder")
        plt.ylabel("avvik / meter")
        plt.savefig("Oppgave1pdfer\opg1dAvvik.pdf")
        plt.figure("1dBane")
        plt.title("Bane med variabelt tidssteg")
        plt.plot(xValueArray, yValueArray, 'ro', markersize=0.4)
        plt.plot([analyticEndpoint[0]], [analyticEndpoint[1]], 'ko', markersize=4, label="analytisk ende")
        plt.legend(loc="best")
        plt.savefig("Oppgave1pdfer\opg1dBane.pdf")
    else:
        L = 1.0E+02
        timeNow = 0
        totalTime = 48* 3600
        deviationLimit = 0.1
        # X = [x, y, dx, dy]
        X = np.array([L, 0, 0, 0])
        #Vi ønsker å logge både tidssteg, avvik for steg som godkjennes samt coordinater
        timeStepArray = []
        deviationArray = []
        coordinateArray = []
        analyticEndpoint = analyticSolution(L, alpha, m)
        while (timeNow < totalTime):
            if timeNow==0:
                h=1000
                Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
                Xtrap = rk2(X,fForEq1,h,timeNow)
                thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
                while thisDeviation > deviationLimit:
                    #Først finjusteres h slik at det første steget akkurat har kommer innenfor nøyaktighetskravet
                    h-=1
                    Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
                    Xtrap = rk2(X, fForEq1, h, timeNow)
                    thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
            h = min(h, totalTime - timeNow)
            Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
            Xtrap = rk2(X, fForEq1, h, timeNow)
            thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
            while (thisDeviation > deviationLimit):
                h = h/2
                Xeul = eulerForEq1(X[:2], timeNow, Vwater, X[2:], h)
                Xtrap = rk2(X,fForEq1,h,timeNow)
                thisDeviation = np.linalg.norm(Xtrap[:2] - Xeul)
            timeStepArray.append([timeNow, h])
            deviationArray.append(thisDeviation)
            timeNow += h
            if (thisDeviation < deviationLimit/10):
                h *= 2
            coordinateArray.append(Xtrap)
            X = Xtrap
        xValueArray = [k[0] for k in coordinateArray]
        yValueArray = [k[1] for k in coordinateArray]
        timeStepArray_xValues = [k[0] for k in timeStepArray]
        timeStepArray_yValues = [k[1] for k in timeStepArray]
        plt.figure("1dtidssteg")
        plt.title("Tidssteg")
        plt.plot(timeStepArray_xValues, timeStepArray_yValues, 'go', markersize=1)
        plt.xlabel("tid / sekunder")
        plt.ylabel("tidssteg / sekunder")
        plt.figure("1davvik")
        plt.title("Avvik mellom Euler og Trapes")
        plt.plot(timeStepArray_xValues, deviationArray, 'ko', markersize=1)
        plt.xlabel("tid / sekunder")
        plt.ylabel("avvik / meter")
        plt.figure("1dbane")
        plt.title("Bane med variabelt tidssteg")
        plt.plot(xValueArray, yValueArray, 'ro', markersize=0.4)
        plt.plot([analyticEndpoint[0]], [analyticEndpoint[1]], 'ko', markersize=4, label="analytisk ende")
        plt.legend(loc="best")
        #plt.show()


def oppgave1(saving=False):
    if saving:
         task1a(True)
         task1b(True)
         task1c(True)
         task1d(True)
    else:
        task1a()
        task1b()
        task1c()
        task1d()
        plt.show()

