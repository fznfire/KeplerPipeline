import numpy as np
import batman
import matplotlib.pyplot as plt


def InjectSignal(t,Flux):
    NoiseLevel = 1e-2

    t = np.arange(0,90,0.0204)
    print t
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = 0.0                       #time of inferior conjunction
    params.per = 10                       #orbital period
    params.rp = 0.05                      #planet radius (in units of stellar radii)
    params.a = 10                         #semi-major axis (in units of stellar radii)
    params.inc = 90                       #orbital inclination (in degrees)
    params.ecc = 0.01                     #eccentricity
    params.w = 45                         #longitude of periastron (in degrees)
    params.limb_dark = "uniform"

    #ld_options = ["uniform", "linear", "quadratic", "nonlinear"]
    #ld_coefficients = [[], [0.3], [0.1, 0.3], [0.5, 0.1, 0.1, -0.1]]

    params.u = []
    m = batman.TransitModel(params, t)
    ModelFlux = m.light_curve(params)


    plt.figure(figsize=(20,5))
    plt.plot(t,ModelFlux,"ko",MarkerSize=2)
    plt.xlabel("Flux(F)")
    plt.ylabel("Days")
    plt.show()
