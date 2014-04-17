include("para_functions.jl")		# this is faster than functions_2.jl
@everywhere using Distributions
#using PyCall
#@pyimport matplotlib.pyplot as pl

println("Use ",nprocs(), " processors")

function getSimulatedArray(n::Int64,p::Float64,d::Float64)
        TIME = zeros(n)
        for i=1:n
                TIME[i] += i
        end
        FLUX = randn(n)*0.1
        e = [50+p*i for i=0:ifloor((n-50)/p)]           # all the starting time of transits
        for e_idx = 1:length(e)
                for i=e[e_idx]:e[e_idx]+d
                        FLUX[i] += -1.+rand()*0.1
                end
        end
        return (TIME,FLUX)
end

# generate a simulated TIME and FLUX array with period = 400, duration = 25, epoch = 50
(TIME,FLUX) = getSimulatedArray(1000,400.,25.)

# use Aigrian & Irwin's method to search for the period, epoch and duration that gives the biggest log10(Q) value (see more about this in README)
# the input parameters are (TIME, FLUX, length_p,f_min, f_max), where length_f is the length of the trial frequency array,
# f_min and f_max are the minimum and maximum trial frequency. f_min must >= 1/total time span.

(length_f,f_min, f_max) = (100,1/403.,1/397.)
tic()
transit_detection!(TIME,FLUX,length_f,f_min, f_max)
println("Total time used in transit_detection!():",toq())

