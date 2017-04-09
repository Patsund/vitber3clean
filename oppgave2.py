import math
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
from scipy.interpolate import RectBivariateSpline
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyproj
import time
plt.style.use('bmh')


def oppgave2(datapath, savefig=False):
    global d, f
    d  = xr.open_dataset(datapath)
    f  = Interpolator(dataset = d)
    if savefig:
        print("**Oppgave 2a**")
        task2a(True)
        print("**Oppgave 2b**")
        task2b(True)
    else:
        print("**Oppgave 2a**")
        task2a()
        print("**Oppgave 2b**")
        task2b()
        plt.show()

class Interpolator():
    def __init__(self, dataset):
        self.dataset = dataset

    def get_interpolators(self, X, it):
        # Add a buffer of cells around the extent of the particle cloud
        buf  = 3
        # Find extent of particle cloud in terms of indices
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

def Vwater(t,X):
    return f(X,t)

def rk2(X, t, h, Vwater):
    h_seconds = h / np.timedelta64(1, 's')
    Xvel = np.array(Vwater(t, X)).reshape(X.shape)           #finner hastighet i nåværende punkt
    Xnext = X+h_seconds*Xvel                                #finner approksimert neste punkt
    XvelNext = np.array(Vwater(t+h, Xnext)).reshape(X.shape) #finner hastighet i approksimert neste punkt
    Xfinal = X+h_seconds/2*(Xvel+XvelNext)
    return Xfinal

def particleTrajectory(X0, time_final, h, time_initial, velocityField, integrator):
    numberOfTimeSteps = int((time_final - time_initial) / h)
    X = np.zeros((numberOfTimeSteps + 1, *X0.shape))
    X[0, :] = X0
    time_now = time_initial
    for step in range(numberOfTimeSteps):
        time_now += h
        X[step + 1, :] = integrator(X[step, :], time_now, h, velocityField)
    return X

def task2a(savefig=False):
    if savefig:
        clockStart = time.time()
        numberOfParticles = 1
        X0 = np.array([-3e6, -1.2e6]).reshape(2, numberOfParticles)
        startTimes = [np.datetime64('2017-02-01T12:00:00'), np.datetime64('2017-02-05T12:00:00'), np.datetime64('2017-02-07T12:00:00')]
        endTimes = [np.datetime64('2017-02-11T12:00:00'), np.datetime64('2017-02-15T12:00:00'), np.datetime64('2017-02-17T12:00:00')]
        dateList = [1, 5, 7]
        plt.figure("2a")
        for day in range(3):
            t0, tEnd = startTimes[day], endTimes[day]
            h = np.timedelta64(3600, 's')
            trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
            xArray = trajectories[:, 0, :]
            yArray = trajectories[:, 1, :]
            print("Tid brukt:", time.time() - clockStart)
            print("Plotter bane nr", day + 1)
            for index in range(numberOfParticles):
                plt.plot(xArray[:, index], yArray[:, index], label=(str(dateList[day])+". feb"))
            plt.legend(bbox_to_anchor=(0., 1.1, 1., .052), loc="best", ncol=3, mode="expand", borderaxespad=0.)
            #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title("Bane for ulike startdager")
        print("Tid brukt:", time.time() - clockStart)
        print("Plott lagret")
        plt.savefig("Oppgave2pdfer\oppgave2a.pdf")
    else:
        clockStart = time.time()
        numberOfParticles = 1
        X0 = np.array([-3e6, -1.2e6]).reshape(2, numberOfParticles)
        startTimes = [np.datetime64('2017-02-01T12:00:00'), np.datetime64('2017-02-05T12:00:00'), np.datetime64('2017-02-07T12:00:00')]
        endTimes = [np.datetime64('2017-02-11T12:00:00'), np.datetime64('2017-02-15T12:00:00'), np.datetime64('2017-02-17T12:00:00')]
        dateList = [1, 5, 7]
        plt.figure("2a")
        for day in range(3):
            t0, tEnd = startTimes[day], endTimes[day]
            h = np.timedelta64(3600, 's')
            trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
            xArray = trajectories[:, 0, :]
            yArray = trajectories[:, 1, :]
            print("Tid brukt:", time.time() - clockStart)
            print("Plotter bane nr", day + 1)
            for index in range(numberOfParticles):
                plt.plot(xArray[:, index], yArray[:, index], label=(str(dateList[day])+". feb"))
            plt.legend(bbox_to_anchor=(0., 1.1, 1., .052), loc="best", ncol=3, mode="expand", borderaxespad=0.)
            #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title("Bane for ulike startdager")
        print("Tid brukt:", time.time() - clockStart)
        print("2aplot klart")

def task2b(savefig=False):
    if savefig:
        startTime = time.time()
        numberOfParticles = 3
        X0 = np.array([-3e6, -3e6, -3e6, -1.2e6, -1.23e6, -1.26e6]).reshape(2, numberOfParticles)
        t0 = np.datetime64('2017-02-01T12:00:00')
        tEnd = np.datetime64('2017-02-11T12:00:00')
        h = np.timedelta64(3600, 's')
        trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
        xArray = trajectories[:, 0, :]
        yArray = trajectories[:, 1, :]
        plt.figure("2b")
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
        ax.add_feature(land_10m)
        ax.coastlines(resolution='10m')
        p1 = pyproj.Proj(d.projection_stere.proj4)
        p2 = pyproj.Proj(proj='latlong')
        plt.title("Ulike startposisjoner 1. feb")
        ax.set_extent((0.7, 5.7, 58.8, 61.4))
        for index in range(numberOfParticles):
            ##############################################################################################
            lons, lats = pyproj.transform(p1, p2, xArray[:,index], yArray[:,index])
            ax.plot(lons, lats,".", transform=ccrs.PlateCarree(), zorder=2)
        plt.savefig("Oppgave2pdfer\oppgave2b.pdf")
        print("oppgave2b.pdf ferdiglagret")
        endTime=time.time()
        print("tid brukt på 2b",endTime-startTime)
    else:
        startTime = time.time()
        numberOfParticles = 3
        X0 = np.array([-3e6, -3e6, -3e6, -1.2e6, -1.23e6, -1.26e6]).reshape(2, numberOfParticles)
        t0 = np.datetime64('2017-02-01T12:00:00')
        tEnd = np.datetime64('2017-02-11T12:00:00')
        h = np.timedelta64(3600, 's')
        trajectories = particleTrajectory(X0, tEnd, h, t0, Vwater, rk2)
        xArray = trajectories[:, 0, :]
        yArray = trajectories[:, 1, :]
        plt.figure("2b")
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        land_10m = cfeature.NaturalEarthFeature('physical','land','10m',color='#dddddd')
        ax.add_feature(land_10m)
        ax.coastlines(resolution='10m')
        p1 = pyproj.Proj(d.projection_stere.proj4)
        p2 = pyproj.Proj(proj='latlong')
        plt.title("Ulike startposisjoner 1. feb")
        ax.set_extent((0.7, 5.7, 58.8, 61.4))
        for index in range(numberOfParticles):
            ##############################################################################################
            lons, lats = pyproj.transform(p1, p2, xArray[:,index], yArray[:,index])
            ax.plot(lons, lats,".", transform=ccrs.PlateCarree(), zorder=2)
        print("2bplot klart")
        endTime=time.time()
        print("tid brukt på 2b",endTime-startTime)
