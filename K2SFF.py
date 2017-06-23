import os
import matplotlib.pyplot as pl
import numpy as np
from lmfit import minimize, Parameters

# pipeline files
from auxiliaries import *
from numpy import mean, cov, cumsum, dot, linalg, size, flipud
from scipy.integrate import quad

def example():
    x = np.linspace(0,20,500)
    m = -2
    c = 5
    y =  m*x  + 5 + 10*np.random.random(len(x))

    print np.std(10*np.random.random(len(x)))

    X = np.array(zip(x,y))
    Cov = np.cov(x,y)
    eig = np.linalg.eig(Cov)


    y1,x1 = eig[1][0]
    y2,x2 =  eig[1][0]
    y2 = -y2



    print "x1, y1::", x1, y1

    theta = np.arctan(y1/x1)

    if x1>0 and y1>0:
        theta = np.pi/2 - theta
    print "Angle:",np.rad2deg(theta)
    print "Angle:",np.rad2deg(np.arctan(-2))

    RotMatrix = [np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]
    NewPoints = np.dot(X,RotMatrix)


    length = (max(x)- min(x)) / 5



    pl.figure()
    pl.plot(x,y,"k.")
    pl.arrow(np.median(x), np.median(y),x1*length,y1*length,fc="k", ec="k",head_width=0.5*length, head_length=0.75*length )
    pl.arrow(np.median(x), np.median(y),x2*length,y2*length,fc="k", ec="k",head_width=0.5*length, head_length=0.75*length )
    pl.show()

    print np.std(NewPoints[:,0]-np.median(NewPoints[:,0]))
    NewPoints[:,0]-np.median(NewPoints[:,0])
    pl.figure(figsize=(15,6))
    pl.subplot(121)
    pl.plot(NewPoints[:,0]-np.median(NewPoints[:,0]),NewPoints[:,1] - np.median(NewPoints[:,1]),"k.")
    pl.axis('equal')
    pl.title('Transformed coordinates')
    pl.subplot(122)
    pl.plot(x-np.median(x),y-np.median(y),"bo")
    pl.axis('equal')
    pl.title('Original coordinates')
    pl.show()


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield list(l[i:i+n])

def median_filter(time,data,binsize=100):
  # do a running median filter dividing all data points by the median of their immediate surroundings
  i = 0
  data_filtered = []
  while i < len(time):
	  bin_begin = max(0, int(i - binsize/2))
	  bin_end = int(min(len(time),(i+binsize/2)))
	  the_bin = data[bin_begin:bin_end]
	  the_bin = sorted(the_bin)
	  median = np.median(the_bin) #[len(the_bin)/2]
	  data_filtered.append(data[i]/median)
	  i = i + 1
  return data_filtered

def sff_residual(params,time,data,s,X, Y, robust=True):
  #
  # residual function used for calculating a fit to centroid (and time), borrowed from reducing data for Spitzer
  #

  # unpack all parameters (note: some may be fixed rather than variable)
  S1 = params['S1'].value
  S2 = params['S2'].value
  S3 = params['S3'].value
  X1 = params['X1'].value
  X2 = params['X2'].value
  X3 = params['X3'].value
  Y1 = params['Y1'].value
  Y2 = params['Y2'].value
  Y3 = params['Y3'].value
  T0 = params['T0'].value
  T1 = params['T1'].value
  T2 = params['T2'].value
  T3 = params['T3'].value
  T4 = params['T4'].value
  TsinAmp = params['TsinAmp'].value
  TsinOff = params['TsinOff'].value

  mean_s = np.mean(s)
  mean_X = np.mean(X)
  mean_Y = np.mean(Y)
  time0 = time[0] - 1.#np.array(1994.0) #time[0]-1.

  model = (T0 + TsinAmp*np.sin((time-time0)+TsinOff) + T1*(time-time0) + T2*((time-time0)**2.) + T3*((time-time0)**3.) + T4*((time-time0)**4.) + S1*(s-mean_s) + S2*((s-mean_s)**2) + S3*((s-mean_s)**3) + X1*(X-mean_X) + X2*((X-mean_X)**2) + X3*((X-mean_X)**3)+ Y1*(Y-mean_Y) + Y2*((Y-mean_Y)**2) + Y3*((Y-mean_Y)**3))

  residual = np.array(data-model)

  if robust:
    # calculate residual in a robust way
    residual = residual[np.abs(residual) < np.mean(residual) + 3.*np.std(residual)]
    #if len(residual2) >= 15:
    #    residual = residual2
  return residual


def sff_fit(time,flux,Xc,Yc,starname='',outputpath='',chunksize=200,  PolyDeg=5):
  #
  # Fit a polynomial to the data and return corrected data
  #
  outputfolder = os.path.join(outputpath,str(starname))

  flux = np.array(flux) / np.median(flux)

  params = Parameters() # fitting parameters, set to vary=false to fix
  params.add('X1', value = 0.,vary=True)
  params.add('X2', value = 0.,vary=True)
  params.add('X3', value = 0.,vary=True)
  params.add('Y1', value = 0.,vary=True)
  params.add('Y2', value = 0.,vary=True)
  params.add('Y3', value = 0.,vary=True)
  params.add('S1', value = 0.,vary=True)
  params.add('S2', value = 0.,vary=True)
  params.add('S3', value = 0.,vary=True)
  params.add('T0', value = 0.,vary=True)
  params.add('T1', value = 0.,vary=True)
  params.add('T2', value = 0.,vary=True) #
  params.add('T3', value = 0.,vary=True) #
  params.add('T4', value = 0.,vary=False)
  params.add('TsinAmp', value = 0.,vary=False)
  params.add('TsinOff', value = 0.,vary=False)
  # first divide data in different chunks

  print chunksize
  time_chunks = list(chunks(time,chunksize))
  flux_chunks = list(chunks(flux,chunksize))
  Xc_chunks = list(chunks(Xc,chunksize))
  Yc_chunks = list(chunks(Yc,chunksize))


  #if the last chunk is really small
  if len(time_chunks[-1])<chunksize/2.5:
      time_chunks[-2].extend(time_chunks[-1])
      time_chunks.pop()
      flux_chunks[-2].extend(flux_chunks[-1])
      flux_chunks.pop()
      Xc_chunks[-2].extend(Xc_chunks[-1])
      Xc_chunks.pop()
      Yc_chunks[-2].extend(Yc_chunks[-1])
      Yc_chunks.pop()

  N_iter = 1
  counter = 1

  #convert the centroid movement to arclength motion
  CorrectedTime = []
  CorrectedFlux = []
  TotalArcLength = []
  for t,X,Y,Flux  in zip(time_chunks,Xc_chunks,Yc_chunks,flux_chunks):

      X = np.array(X)
      Y = np.array(Y)
      Flux = np.array(Flux)/np.median(Flux)

      #Calculate the drift in the field
      #TODO check how having 1.5 changes the graphs
      X = X - np.median(X)+1.5
      Y = Y - np.median(Y)+1.5

      Xc_Yc =  np.array(zip(X,Y))
      Cov = np.cov(X,Y)
      eig = np.linalg.eig(Cov)

      y1,x1 = eig[1][0]
      theta = np.arctan(y1/x1)

      if x1>0 and y1>0:
          theta = 2*np.pi - theta

      RotMatrix = [np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]

      TempHolder = np.dot(Xc_Yc,RotMatrix)
      X_Transformed, Y_Transformed = TempHolder[:,0], TempHolder[:,1]
      X_Transformed = X_Transformed - np.median(X_Transformed)
      Y_Transformed = Y_Transformed - np.median(Y_Transformed)
      del TempHolder

      params = np.polyfit(X_Transformed, Y_Transformed,PolyDeg)


      #calculate the arc length instead of the star drift
      Mult = np.arange(len(params)-1,0,-1)


      NewParams = params[:-1]*Mult

      def Integrand(x,params):
        return np.sqrt(1+(np.polyval(params,x))**2)

      #for XVal in x:

      ArcLength = []
      for XVal in X_Transformed:
          tempArc,_ = quad(Integrand,0,XVal,args=(NewParams),epsabs=1e-5)
          ArcLength.append(tempArc)

      TotalArcLength.extend(ArcLength)

      ##TODO remove the thruster events


      #calculate ds/dt at every point





      ParamFlux = np.polyfit(ArcLength,Flux,PolyDeg)
      FluxPredicted = np.polyval(ParamFlux,ArcLength)
      CorrectedTime.extend(t)
      CorrectedFlux.extend(np.array(Flux)/FluxPredicted)

      '''
      #TO DO A printing flag
      pl.figure()
      pl.plot(X,Y,"k.-")
      pl.plot(X_Transformed, Y_Transformed, "b.")
      pl.plot(np.linspace(min(X_Transformed), max(X_Transformed),1000), np.polyval(params,np.linspace(min(X_Transformed), max(X_Transformed),1000)),"r-")
      pl.axvline(x=0,color="black",lw=2)
      pl.axhline(y=0,color="black",lw=2)
      pl.title('Drift in the frame:'+str(counter))
      pl.show()

      counter+=1

      pl.figure()
      pl.plot(ArcLength,Flux, "ko")
      pl.plot(np.linspace(min(ArcLength), max(ArcLength),1000),np.polyval(ParamFlux, np.linspace(min(ArcLength), max(ArcLength),1000)),"r-")
      pl.axvline(x=0,color="green",lw=2)
      pl.axhline(y=1,color="green",lw=2)
      pl.title("Arclength vs Flux")
      pl.show()

      pl.figure()
      pl.plot(t,Flux, "ko")
      pl.plot(t,FluxPredicted,"bo")
      pl.show()
      '''
  pl.plot()
  pl.plot()
  Diff = np.diff(TotalArcLength)
  pl.figure()
  pl.plot(Diff,"k.", MarkerSize=2)
  pl.show()
  return [CorrectedTime,CorrectedFlux]
