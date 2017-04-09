import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import xarray as xr
from scipy.interpolate import RectBivariateSpline
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyproj
import random
import time
plt.style.use('bmh')
startTime=time.time()

def oppgave3(datapath, saving=False):
    opg3Time=time.time()
    global d, f
    d  = xr.open_dataset(datapath)
    f  = Interpolator(dataset = d)
    if saving:
        #print("**Oppgave 3a med separate plot**")
        #task3a(separate=True, savefig=True)
        print("**Oppgave 3a med samlet plot**")
        task3a(savefig=True)
        print("**Oppgave 3b**")
        task3b(savefig=True)
        print("brukte",time.time()-opg3Time,"sekunder på oppgave 3")
    else:
        print("**Oppgave 3a**")
        task3a()
        print("**Oppgave 3b**")
        task3b()
        print("brukte",time.time()-opg3Time,"sekunder på oppgave 3")
        plt.show()

#Denne funksjonen gir oss vannhastigheten med bakgrunn i dataene
class Interpolator():
    def __init__(self, dataset):
        self.dataset = dataset

    def get_interpolators(self, X, it):
        # Add a buffer of cells around the extent of the particle butt
        buf  = 3
        # Find extent of particle butt in terms of indices
        imax = np.searchsorted(self.dataset.X, np.amax(X[0,:])) + buf
        imin = np.searchsorted(self.dataset.X, np.amin(X[0,:])) - buf
        jmax = np.searchsorted(self.dataset.Y, np.amax(X[1,:])) + buf
        jmin = np.searchsorted(self.dataset.Y, np.amin(X[1,:])) - buf
        # Take out subset of array, to pass to
        # interpolation object
        # Fill NaN values (land cells) with 0, otherwise
        # interpolation won't work
        u    = self.dataset.u[it, 0, jmin:jmax, imin:imax].fillna(0.0)
        v    = self.dataset.v[it, 0, jmin:jmax, imin:imax].fillna(0.0)
        # RectBivariateSpline essentially returns a function,
        # which can be called to get value at arbitrary position
        # kx and ky sets order of spline interpolation along either direction (must be 1 <= kx <= 5)
        # transpose arrays to switch order of coordinates
        fu   = RectBivariateSpline(self.dataset.X[imin:imax], self.dataset.Y[jmin:jmax], u.T)#, kx = 3, ky = 3)
        fv   = RectBivariateSpline(self.dataset.X[imin:imax], self.dataset.Y[jmin:jmax], v.T)#, kx = 3, ky = 3)
        return fu, fv

    def get_time_index(self, t):
        # Get index of largest timestamp smaller than (or equal to) t
        return np.searchsorted(self.dataset.time, t, side='right') - 1

    def __call__(self, X, t):
        # get index of current time in dataset
        it = self.get_time_index(t)
        # get interpolating functions,
        # covering the extent of the particle
        fu, fv = self.get_interpolators(X, it)
        # Evaluate velocity at position(x[:], y[:])
        dx = fu(X[0,:], X[1,:], grid = False)
        dy = fv(X[0,:], X[1,:], grid = False)
        return np.array([dx, dy])


t = np.datetime64('2017-02-01T12:00:00')
X = np.array([-3000000, -1200000]).reshape(2, 1)  #reshape (2,Np)
#print("Hastighet i dette punktet ved denne tiden:\n",f(X, t))

t0   = np.datetime64('2017-02-01T12:00:00')
h  = np.timedelta64(3600, 's')
#Step forward
i = 2
t = t0 + i*h
#Get number of seconds in h
h_seconds = h / np.timedelta64(1, 's')

def Vwater(t,X):
    return f(X,t)

def rk2(X,t,h,Vwater):
    h_seconds = h / np.timedelta64(1, 's')
    Xvel=np.array(Vwater(t,X)).reshape(X.shape) #finner hastighet i nåværende punkt
    Xnext=X+h_seconds*Xvel #finner approksimert neste punkt
    XvelNext=np.array(Vwater(t+h,Xnext)).reshape(X.shape) #finner hastighet i approksimert neste punkt
    Xfinal=X+h_seconds/2*(Xvel+XvelNext)
    return Xfinal

def particleTrajectory(X0, time_final, h, time_initial, velocityField, integrator):
    numberOfTimeSteps = int((time_final - time_initial) / h)
    X = np.zeros((numberOfTimeSteps + 1, *X0.shape))
    # X lagrer alle partiklenes posisjon gjennom hele tidsforløpet
    X[0, :] = X0
    time_now = time_initial
    for step in range(numberOfTimeSteps):  # Tok bort +1 pga hele timer
        X[step + 1, :] = integrator(X[step, :], time_now, h, velocityField)
        time_now += h  # numpy håndterer time64-opplegg
        # Denne bør virkelig returnere koordinater på formen (2,1) - og det gjør den
    return X

def randomX0Array(lowendY,highendY,lowendX,highendX,numberOfParticles):
    finalArray=[]
    for i in range(numberOfParticles):
        finalArray.append(random.randint(lowendY,highendY))
    for i in range(numberOfParticles):
        finalArray.append(random.randint(lowendX,highendX))
    return finalArray

def task3a(separate=False, savefig=False):
    if savefig:
        if separate:
            numberOfParticles = 10000
            initArray=randomX0Array(-3.01e6,-2.99e6,-1.21e6,-1.19e6,numberOfParticles)
            X0 = np.array(initArray).reshape(2, numberOfParticles)  # reshape (2,Np)
            t0 = np.datetime64('2017-02-01T12:00:00')
            tEnd = np.datetime64('2017-02-11T12:00:00')
            h = np.timedelta64(3600, 's')
            trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
            trajectories = np.hsplit(trajectories, len(trajectories[0]))
            xArray = trajectories[0]
            yArray = trajectories[1]
            plt.figure()
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
            ax.add_feature(land_10m)
            ax.coastlines(resolution='10m')
            ax.set_extent((0, 6, 58.5, 62.5))
            p1 = pyproj.Proj(d.projection_stere.proj4)
            p2 = pyproj.Proj(proj='latlong')
            plt.title("Partikkelens bane")
            colors=["r.","b.","g.","c.","k.","y."]
            for index in range(0,6): #endret fra len(xArrayParticleSplit)
                plt.figure("3a dag"+str(index*2))
                ax = plt.axes(projection=ccrs.NorthPolarStereo())
                land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
                ax.add_feature(land_10m)
                ax.coastlines(resolution='10m')
                ax.set_extent((0, 6, 58.5, 62.5))
                lons, lats = pyproj.transform(p1, p2, xArray[index*2*24], yArray[index*2*24])
                ax.plot(lons, lats, colors[index], transform=ccrs.PlateCarree(), zorder=2)
                print("lagrer figur")
                figTime=time.time()
                plt.savefig("Oppgave3pdfer\opg3adag"+str(index*2)+"plot.png")
                plt.close()
                print("figur lagret")
        else:
            numberOfParticles = 10000
            initArray=randomX0Array(-3.01e6,-2.99e6,-1.21e6,-1.19e6,numberOfParticles)
            X0 = np.array(initArray).reshape(2, numberOfParticles)  # reshape (2,Np)
            t0 = np.datetime64('2017-02-01T12:00:00')
            tEnd = np.datetime64('2017-02-11T12:00:00')
            h = np.timedelta64(3600, 's')
            trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
            trajectories = np.hsplit(trajectories, len(trajectories[0]))
            xArray = trajectories[0]
            yArray = trajectories[1]
            plt.figure()
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
            ax.add_feature(land_10m)
            ax.coastlines(resolution='10m')
            ax.set_extent((0, 6, 58.5, 62.5))
            p1 = pyproj.Proj(d.projection_stere.proj4)
            p2 = pyproj.Proj(proj='latlong')
            plt.title("Partikkelens bane")
            colors=["r.","b.","g.","c.","k.","y."]
            for index in range(0,6):
                lons, lats = pyproj.transform(p1, p2, xArray[index*2*24], yArray[index*2*24])
                ax.plot(lons, lats, colors[index], transform=ccrs.PlateCarree(), zorder=2)
            day0 = mpatches.Patch(color='r', label='Day 0')
            day2 = mpatches.Patch(color='b', label='Day 2')
            day4 = mpatches.Patch(color='g', label='Day 4')
            day6 = mpatches.Patch(color='c', label='Day 6')
            day8 = mpatches.Patch(color='k', label='Day 8')
            day10 = mpatches.Patch(color='y', label='Day 10')
            ax.legend(handles=[day0,day2,day4,day6,day8,day10],bbox_to_anchor=(1.05,1), loc = 2, borderaxespad=0.)
            print("lagrer figur")
            plt.savefig("Oppgave3pdfer\opg3atotalfigur3a.png")
            plt.close()
            print("figur lagret")
            #endTime=time.time()
    else:
        if separate:
            numberOfParticles = 10000
            initArray=randomX0Array(-3.01e6,-2.99e6,-1.21e6,-1.19e6,numberOfParticles)
            X0 = np.array(initArray).reshape(2, numberOfParticles)  # reshape (2,Np)
            t0 = np.datetime64('2017-02-01T12:00:00')
            tEnd = np.datetime64('2017-02-11T12:00:00')
            h = np.timedelta64(3600, 's')
            trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
            trajectories = np.hsplit(trajectories, len(trajectories[0]))
            xArray = trajectories[0]
            yArray = trajectories[1]
            plt.figure()
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
            ax.add_feature(land_10m)
            ax.coastlines(resolution='10m')
            ax.set_extent((0, 6, 58.5, 62.5))
            p1 = pyproj.Proj(d.projection_stere.proj4)
            p2 = pyproj.Proj(proj='latlong')
            plt.title("Partikkelens bane")
            colors=["r.","b.","g.","c.","k.","y."]
            for index in range(0,6): #endret fra len(xArrayParticleSplit)
                plt.figure("3a dag"+str(index*2))
                ax = plt.axes(projection=ccrs.NorthPolarStereo())
                land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
                ax.add_feature(land_10m)
                ax.coastlines(resolution='10m')
                ax.set_extent((0, 6, 58.5, 62.5))
                lons, lats = pyproj.transform(p1, p2, xArray[index*2*24], yArray[index*2*24])
                ax.plot(lons, lats, colors[index], transform=ccrs.PlateCarree(), zorder=2)
            #plt.show()
        else:
            numberOfParticles = 10000
            initArray=randomX0Array(-3.01e6,-2.99e6,-1.21e6,-1.19e6,numberOfParticles)
            X0 = np.array(initArray).reshape(2, numberOfParticles)  # reshape (2,Np)
            t0 = np.datetime64('2017-02-01T12:00:00')
            tEnd = np.datetime64('2017-02-11T12:00:00')
            h = np.timedelta64(3600, 's')
            trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
            trajectories = np.hsplit(trajectories, len(trajectories[0]))
            xArray = trajectories[0]
            yArray = trajectories[1]
            plt.figure("3a samlet")
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
            ax.add_feature(land_10m)
            ax.coastlines(resolution='10m')
            ax.set_extent((0, 6, 58.5, 62.5))
            p1 = pyproj.Proj(d.projection_stere.proj4)
            p2 = pyproj.Proj(proj='latlong')
            plt.title("Partikkelens bane")
            colors=["r.","b.","g.","c.","k.","y."]
            for index in range(0,6):
                lons, lats = pyproj.transform(p1, p2, xArray[index*2*24], yArray[index*2*24])
                ax.plot(lons, lats, colors[index], transform=ccrs.PlateCarree(), zorder=2)
            day0 = mpatches.Patch(color='r', label='Day 0')
            day2 = mpatches.Patch(color='b', label='Day 2')
            day4 = mpatches.Patch(color='g', label='Day 4')
            day6 = mpatches.Patch(color='c', label='Day 6')
            day8 = mpatches.Patch(color='k', label='Day 8')
            day10 = mpatches.Patch(color='y', label='Day 10')
            ax.legend(handles=[day0,day2,day4,day6,day8,day10],bbox_to_anchor=(1.05,1), loc = 2, borderaxespad=0.)
            #plt.show()
            #endTime=time.time()

def frequencyCounter(X,Nx,Ny,numberOfParticles):
    cellSize = 800
    x = np.array((  X[0, :]+3.01e6)/cellSize, dtype= int)
    y = np.array((X[1, :]+1.3e6) / cellSize, dtype= int)
    frequencyGrid = np.zeros((Ny, Nx))
    for i in range(numberOfParticles):
        frequencyGrid[y[i], x[i]] += 1
    return frequencyGrid

def task3b(savefig=False):
    if savefig:
        Nx, Ny = 600, 300 #antall ruter i hhv x- og y-retning
        x = -3010000 + 800 * np.arange(Nx) #deler inn i ruter
        y = -1300000 + 800 * np.arange(Ny)
        x, y = np.meshgrid(x, y) #setter opp i et rutenett
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
        ax.add_feature(land_10m)
        ax.coastlines(resolution='10m')
        p1 = pyproj.Proj(d.projection_stere.proj4)
        p2 = pyproj.Proj(proj='latlong')
        lons, lats = pyproj.transform(p1, p2, x, y)
        numberOfParticles=100000
        X0 = np.array(randomX0Array(-3.01e6, -2.99e6, -1.21e6, -1.19e6, numberOfParticles)).reshape(2, numberOfParticles)
        print("laget X0",time.time()-startTime)
        t0 = np.datetime64('2017-02-01T12:00:00')
        tEnd = np.datetime64('2017-02-11T12:00:00')
        h = np.timedelta64(3600, 's')
        trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
        xArray = trajectories[:, 0, :]
        yArray = trajectories[:, 1, :]
        dataGrid_X = d.X.values
        dataGrid_Y = d.Y.values
        ax.set_extent((0, 6, 58.5, 62.5))
        colormap = ['Reds', 'Oranges', 'Greens', 'Blues', 'Purples', 'PuRd'] #['autumn', 'cool', 'copper' , 'summer', 'winter', 'autumn']
        for day in range(0, 6):
            plt.figure("3b dag"+str(day*2))
            plt.title("3b dag"+str(day*2))
            x = -3010000 + 800 * np.arange(Nx)
            y = -1300000 + 800 * np.arange(Ny)
            x, y = np.meshgrid(x, y)
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
            ax.add_feature(land_10m)
            ax.coastlines(resolution='10m')
            ax.set_extent((0, 6, 58.5, 62.5))
            p1 = pyproj.Proj(d.projection_stere.proj4)
            p2 = pyproj.Proj(proj='latlong')
            concentration, coordinates_X, coordinates_Y = np.histogram2d(xArray[day*2*24], yArray[day*2*24], bins=(dataGrid_X, dataGrid_Y))
            concentration = np.ma.masked_array(concentration, mask=concentration == 0)
            coordinates_X, coordinates_Y = np.meshgrid(coordinates_X, coordinates_Y)
            lons, lats = pyproj.transform(p1, p2, coordinates_X, coordinates_Y)
            ax.pcolormesh(lons, lats, concentration.T, transform=ccrs.PlateCarree(), zorder=2, cmap='gist_heat_r')
            print("lagrer figur")
            plt.savefig("Oppgave3pdfer\oppgave3b"+str(2*day)+".png")
            plt.close()
            print("figur lagret")
        #plt.show()
    else:
        Nx, Ny = 600, 300
        x = -3010000 + 800 * np.arange(Nx)
        y = -1300000 + 800 * np.arange(Ny)
        x, y = np.meshgrid(x, y)
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
        ax.add_feature(land_10m)
        ax.coastlines(resolution='10m')
        p1 = pyproj.Proj(d.projection_stere.proj4)
        p2 = pyproj.Proj(proj='latlong')
        lons, lats = pyproj.transform(p1, p2, x, y)
        numberOfParticles=100000
        X0 = np.array(randomX0Array(-3.01e6, -2.99e6, -1.21e6, -1.19e6, numberOfParticles)).reshape(2, numberOfParticles)
        t0 = np.datetime64('2017-02-01T12:00:00')
        tEnd = np.datetime64('2017-02-11T12:00:00')
        h = np.timedelta64(3600, 's')
        trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
        xArray = trajectories[:, 0, :]
        yArray = trajectories[:, 1, :]
        dataGrid_X = d.X.values
        dataGrid_Y = d.Y.values
        ax.set_extent((0, 7, 58, 64))
        colormap = ['Reds', 'Oranges', 'Greens', 'Blues', 'Purples', 'PuRd'] #['autumn', 'cool', 'copper' , 'summer', 'winter', 'autumn']
        for day in range(0, 6):
            plt.figure("3b dag"+str(day*2))
            plt.title("3b dag"+str(day*2))
            x = -3010000 + 800 * np.arange(Nx)
            y = -1300000 + 800 * np.arange(Ny)
            x, y = np.meshgrid(x, y)
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
            ax.add_feature(land_10m)
            ax.coastlines(resolution='10m')
            ax.set_extent((0, 6, 58.5, 62.5))
            p1 = pyproj.Proj(d.projection_stere.proj4)
            p2 = pyproj.Proj(proj='latlong')
            concentration, coordinates_X, coordinates_Y = np.histogram2d(xArray[day*2*24], yArray[day*2*24], bins=(dataGrid_X, dataGrid_Y))
            concentration = np.ma.masked_array(concentration, mask=concentration == 0)
            coordinates_X, coordinates_Y = np.meshgrid(coordinates_X, coordinates_Y)
            lons, lats = pyproj.transform(p1, p2, coordinates_X, coordinates_Y)
            ax.pcolormesh(lons, lats, concentration.T, transform=ccrs.PlateCarree(), zorder=2, cmap='gist_heat_r')
        #plt.show()