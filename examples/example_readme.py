import numpy as np
from TVDCondat2013 import TVD
from matplotlib import pyplot as plt


# Generating n segments of a noisy signal
size_segment = 400
n_segments = 5
C = np.random.rand(size_segment) + np.random.randint(-5,5)
for d in range(n_segments-1):
	A = np.random.rand(size_segment) + np.random.randint(-5,5)
	C = np.concatenate((C,A))

X = np.arange(n_segments*size_segment)

plt.plot(X,C,color= 'r', lw = 0.2, zorder = 1, label= "Raw signal")

# Setting the regulation parameters
lambda_TVD = [0.1,1,10,100,1000]

# Denoising with different lambda

for l in lambda_TVD:
	denoised = TVD(C,l)
	plt.plot(X,denoised, lw = 1, zorder = 2, label= "TVD l = %s"%(l))



plt.xlabel("X")
plt.ylabel("Signal")

plt.legend()

plt.savefig("Example.png", dpi = 500)
