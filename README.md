
# Transit Search

This project searches for transiting planets candidates in Kepler data using Aigran & Irwin's transit search algorithm[1]. The basic idea is: loop through different trial periods, trial transit durations and trial transit epochs, find out which (p, d, e) maximizes the Q value, which is basically the square of the signal to noise ratio (SNR). In our code, we use log10(Q) instead of just Q.  

The core code for transit search is the function: transit_detection!(), whose input parameters are: (TIME,FLUX,length_f,f_min, f_max). TIME and FLUX are the preprocessed time array and flux array, length_f is the number of trial frequencies, f_min, f_max are the minimum and maximum trial frequencies. The output of transit_detection!() is ((best_p,best_d,best_e,best_logQ,logQ_p), i.e. the best estimated period, transit duration and epoch, the corresponding logQ value and the array that stores each logQ value for each trial period.

The main steps of transit_detection!() is: (1) For a given trial period value, do phase folding; (2) For each trial duration and trial epoch, identify the in-transit data points in this phased flux, and calculate the logQ for it; (3) Compare all logQ values, get the maximum logQ, then the corresponding (p, e, d) is the best estimated result.

### Phase folding

Purpose: given the input time and flux array, fold the array using period p. After the folding, the time range is [0,p]. It has two steps: (1) Fold all data points into a time range[0,p]; (2) Merge all data points that're within one small time window (set as t_avg_step), since after the first step, the density of the data points will increase, combine data points that are too close to each other as one new data point can accelerate the algorithm.

### Simulated data test

test_simulated_data.jl uses a simulated data set for a simple test. 
The data shows has a period = 400, duration = 25 and epoch = 50. See the plot SimulatedDataTest.png for this simulated data set.
Run the code with trial frequency range [1/403, 1/397] and number of trial frequency = 100, the result gives: (best_p, best_d, best_e) = (399.5579641376555,26.0,49.0), it is close to the true value.
Run Time:It uses 147.44s to run on a laptop with a 2.9 GHz Intel Core i7 and 8 GB memroy. It uses 160.0s to run on Tesla with 1 core.

### One planet test

transit_search_test.jl is a small, crude test using true kepler data of a confirmed planet.
It uses only 4 usable segments of one lightcurve of a known planet (KOI 97.01), its true period is 4.88548917 days, duration is 4.6296 hours. Run transit_search_test.jl, with input parameters set to: length_f=20,f_min=1/5.,f_max=1/4.7, avg_t_step = 0.003 (this is the time step for trial duration and trial epoch, normally it is determined by the average time gap between two consecutive data points,for our light curve data, it should be about 0.0007 days. But for the purpose of a quick test, we use 0.003 days). 
The best estimated results are: period = 4.890738813735692 days, duration = 4.536 hours
The true values are: period = 4.88548917 days, duration = 5.2295 hours
Run Time: 2225.0 s on Tesla with 1 core.



[1]Aigran & Irwin 2004 (http://arxiv.org/abs/astro-ph/0401393)
