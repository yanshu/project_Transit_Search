include("functions.jl")
using PyCall
@pyimport kplr
@pyimport matplotlib.pyplot as pl

TIME = zeros(1000)
for i=1:1000
	TIME[i] += i
end
FLUX = randn(1000)*0.1
for i=50:75
	FLUX[i] += -1.+rand()*0.1
end
FLUX[50] -= 0.5
FLUX[75] -= 0.5
for i=450:475
	FLUX[i] += -1.+rand()*0.1
end
FLUX[450] -= 0.5
FLUX[475] -= 0.5
for i=950:975
	FLUX[i] += -1.+rand()*0.1
end
FLUX[850] -= 0.5
FLUX[875] -= 0.5

#plot FLUX vs TIME
#plotEverything(TIME,FLUX,"97.01")

# use Aigrian & Irwin's method to search for the period, epoch and duration that gives the biggest log10(Q) value (see more about this in README)
# the input parameters are (TIME, FLUX, length_p,f_min, f_max), where length_f is the length of the trial frequency array,
# f_min and f_max are the minimum and maximum trial frequency. f_min must >= 1/total time span.

(length_f,f_min, f_max) = (100,1/403.,1/397.)
transit_detection!(TIME,FLUX,length_f,f_min, f_max);

