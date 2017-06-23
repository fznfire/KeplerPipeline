import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model, decomposition

x = np.linspace(0,10,100)
y = x**2 + np.random.random(len(x))

pca = decomposition.PCA()

Regression = linear_model.LogisticRegression()
Regression.fit(x,y)
plt.figure()
plt.plot(x,y,'k.')
plt.show()
#Does it affect a sinusoidal graph

#Does it take off long term variation in data
