import matplotlib
matplotlib.use('Agg')
import numpy as np
from TVDCondat2013 import TVD, D_TVD_R
from matplotlib import pyplot as plt


# Generating noisy segmented data

n_segment = 5
size_segment = 200
noise_emphasiizer = 500
A = noise_emphasiizer * np.random.rand(size_segment)
X = np.arange(1,1+size_segment)

a = np.random.randint(-4,4)
for i in range(size_segment):
	A[i] = a*X[i]+A[i]

current_run = 1
lA = A[-1]
for i in range(n_segment-1):
	a = np.random.randint(-4,4)
	tA = noise_emphasiizer * np.random.rand(size_segment)
	tX = np.arange(current_run*size_segment,(current_run+1)*size_segment)
	for i in range(size_segment):
		tA[i] = a*tX[i]+tA[i]+lA
	A = np.concatenate((A,tA))
	X = np.concatenate((X,tX))
	current_run += 1
	lA = tA[-1]

plt.plot(X,A,color= 'r', lw = 0.4, zorder = 1, label= "Raw signal")

# Setting the regulation parameters
lambda_TVD = [100,250,500,5000,10000]

# Denoising with different lambda

for l in lambda_TVD:
	print('lambda: ', l)
	denoised = D_TVD_R(A,l)
	plt.plot(X,denoised, lw = 1, zorder = 2, label= r"TVD $\lambda$ = %s"%(l))



plt.xlabel("X")
plt.ylabel("Signal")

plt.legend()

plt.savefig("Example_curve.png", dpi = 500)
