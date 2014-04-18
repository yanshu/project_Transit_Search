
# Transit Search

This project searches for transiting planets candidates in Kepler data using Aigran & Irwin's transit search algorithm[1]. The basic idea is: loop through different trial periods, trial transit durations and trial transit epochs, find out which (p, d, e) maximizes the Q value, which is basically the square of the signal to noise ratio (SNR). In our code, we use log10(Q) instead of just Q.  

The core code for transit search is the function: transit_detection!(), whose input parameters are: (TIME,FLUX,length_f,f_min, f_max). TIME and FLUX are the preprocessed time array and flux array, length_f is the number of trial frequencies, f_min, f_max are the minimum and maximum trial frequencies. The output of transit_detection!() is ((best_p,best_d,best_e,best_logQ,logQ_p), i.e. the best estimated period, transit duration and epoch, the corresponding logQ value and the array that stores each logQ value for each trial period.

The main steps of transit_detection!() is: (1) For a given trial period value, do phase folding; (2) For each trial duration and trial epoch, identify the in-transit data points in this phased flux, and calculate the logQ for it; (3) Compare all logQ values, get the maximum logQ, then the corresponding (p, e, d) is the best estimated result.

### Phase folding

Purpose: given the input time and flux array, fold the array using period p. After the folding, the time range is [0,p]. It has two steps: (1) Fold all data points into a time range[0,p]; (2) Merge all data points that're within one small time window (set as t_avg_step), since after the first step, the density of the data points will increase, combine data points that are too close to each other as one new data point can accelerate the algorithm.

### Test without parallelization

#### Simulated data test

test_simulated_data.jl uses a simulated data set for a simple test. (To run this test, remember to comment line36 in functions.jl and uncomment line33; because we need the avg_t_step parameter to be set as the default definition)

The simulated data has a period = 400, duration = 25 and epoch = 50. See the plot SimulatedDataTest.png for this simulated data set.
Run the code with trial frequency range [1/403, 1/397] and number of trial frequency = 100, the result gives: (best_p, best_d, best_e) = (399.5579641376555,26.0,49.0), it is close to the true value.

Run Time: It uses 147.44s to run on a laptop with a 2.9 GHz Intel Core i7 and 8 GB memroy. It uses 160.0s to run on Tesla with 1 core. Most of the process time is spent in the get_in_transit_index() function:

(1) number of total input data points: 1000:
Time used in phase folding: 0.09459777300000001
Time used in get_in_transit_index: 190.12969029804088
Time used in calculate_logQ: 8.025332680003013
Total time used in transit_detection!(): 211.536398995 seconds
Result: (399.5579641376555,26.0,49.0)

(2) number of total input data points: 10000:
Time used in phase folding: 0.217275445
Time used in get_in_transit_index: 166.29636909701702
Time used in calculate_logQ: 7.568610243004216
Total time used in transit_detection!():186.867834682
Result: (399.9775,26.0,49.0)

(3) number of total input data points: 100000:
Time used in phase folding: 1.6532172079999992
Time used in get_in_transit_index: 159.94216397200384
Time used in calculate_logQ: 7.574739042001622
Total time used in transit_detection!():180.057972289
Result: (399.9775,26.0,51.0)

We see: (1) Most of the process time is spent in the get_in_transit_index() function;
(2) as the total size of input flux array increases, the total process time increases, and the time spent in phase folding increase almost linearly with input array size, while other functions spent more or less similar time.

So for parallelization, we need to work on the get_in_transit_index!() function.

#### One planet test

transit_search_test.jl is a small, crude test using true kepler data of a confirmed planet.	(To run it, remember to comment line33 in functions.jl and uncomment line36, because here we want t_avg_step to be 0.003 days)

It uses only 4 usable segments of one lightcurve of a known planet (KOI 97.01), its true period is 4.88548917 days, duration is 4.6296 hours. Run transit_search_test.jl, with input parameters set to: length_f=20,f_min=1/5.,f_max=1/4.7, avg_t_step = 0.003 (this is the time step for trial duration and trial epoch, normally it is determined by the average time gap between two consecutive data points,for our light curve data, it should be about 0.0007 days. But for the purpose of a quick test, we use 0.003 days). 

The best estimated results are: period = 4.890738813735692 days, duration = 4.536 hours
The true values are: period = 4.88548917 days, duration = 5.2295 hours

Run Time: 2225.0 s on Tesla with 1 core.


### Parallelization using distributed arrays

para_functions.jl shows the core functions needed for parallelization of the original code. Most of the functions in it are the same with functions.jl , the most important change is in transit_detection!() funcntion, where the 3-layer for loops is modified and a new function get_in_transit_index_and_calc_logQ! is created to make it easy to parallelize:
```
for i=1:length_p            # loop through all trial periods
        for j=1:max_length_d_and_e      # loop through all trial durations for this period
            (best_p,best_d,best_e,best_logQ) = get_in_transit_index_and_calc_logQ!(foldedTIME[i], foldedFLUX[i],p,d,e,i,j,best_logQ,best_p,best_d,best_e,avg_t_step,noise)
        end
    end
```

Run 
```
julia -p n para_simulated_data.jl 
```
to test the parallelization for simulated data.


[1]Aigran & Irwin 2004 (http://arxiv.org/abs/astro-ph/0401393)
