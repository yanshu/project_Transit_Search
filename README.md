
# Transit Search

This project searches for transiting planets candidates in Kepler data using Aigran & Irwin's transit search algorithm[1]. The basic idea is: loop through different trial periods, trial transit durations and trial transit epochs, find out which (p, d, e) maximizes the Q value, which is basically the square of the signal to noise ratio (SNR). In our code, we use log10(Q) instead of just Q.  

functions.jl and para_functions.jl contain all the functions for this project. The former is the serial vesion, and the latter is the parallel version. The core function is: transit_detection!(), whose input parameters are: (TIME,FLUX,length_f,f_min, f_max). TIME and FLUX are the preprocessed time array and flux array, length_f is the number of trial frequencies, f_min, f_max are the minimum and maximum trial frequencies. The output of transit_detection!() is ((best_p,best_d,best_e,best_logQ,logQ_p), i.e. the best estimated period, transit duration and epoch, the corresponding logQ value and the array that stores each logQ value for each trial period.

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

#### Simulated data test

para_functions.jl shows the core functions needed for parallelization of the original code. Most of the functions in it are the same with functions.jl , the most important change is in transit_detection!() funcntion, where the 3-layer for loops is modified and a new function get_in_transit_index_and_calc_logQ! is created to make it easy to parallelize.

The index [i,j] are distributed over several processors.

Run 
```
julia -p n para_simulated_data.jl 
```
to test the parallelization for simulated data. (n is the number of processors)

Run time vs number of processors: (Telsa)

(1) n = 1: Total time used in transit_detection!():235.974000392

(2) n = 2: Total time used in transit_detection!():123.928522869

(3) n = 3: Total time used in transit_detection!():128.286416456

(4) n = 4: Total time used in transit_detection!():99.218811207

(5) n = 5: Total time used in transit_detection!():107.634166859 

(6) n = 6: Total time used in transit_detection!():71.39507179

(7) n = 7: Total time used in transit_detection!():71.575478587

(8) n = 8: Total time used in transit_detection!():57.220656825

(9) n = 9: Total time used in transit_detection!():57.918904511

(10) n = 10: Total time used in transit_detection!():50.570660036

(11) n = 11: Total time used in transit_detection!():48.488129586

We see the increase of processors makes the process time significantly decrease. 
The reason why p=2 & 3, p= 4 & 5, p= 6 & 7 , p = 8 & 9, p = 10 & 11 spent similar time is: The parallelization is using distributed arrays: distribute a 2d array of indices (i,j) into different processors, the default setup is even distribution. But for our transit search, we don't need all indices - for a given period value p[i], the number of the trial durations is int(p[i]/t_avg_step),so we only need the corresponding duration array index to go from 1 to int(p[i]/t_avg_step), so at a large period array index i, we need a larger duration array index j. It's similar to a triangular matrix, so when this 2d array is evenly distributed, this triangular array is not evenly distributed.

To fix this:

Instead of using a 2d distributed array, just use a 1d distributed array to store the index (i,j), i.e. an array of tuples, see the code here:
```
oned_idx = Any[]
for i=1:length_p
    length_d = int(p[i]/avg_t_step)
    for j = 1:length_d
        append!(oned_idx,[(i,j)])
    end
end
idx = distribute(oned_idx)
```
After this modification, the run time for different processors are:

(1) n = 1: Total time used in transit_detection!():215.752591133

(2) n = 2: Total time used in transit_detection!():123.095294461

(3) n = 3: Total time used in transit_detection!():90.721392053

(4) n = 4: Total time used in transit_detection!():75.883206517

(5) n = 5: Total time used in transit_detection!():68.422990043

(6) n = 6: Total time used in transit_detection!():64.156633272

(7) n = 7: Total time used in transit_detection!():54.97174648

(8) n = 8: Total time used in transit_detection!():51.917523137

(9) n = 9: Total time used in transit_detection!():46.779422036

(10) n = 10: Total time used in transit_detection!():49.136618533

(11) n = 11: Total time used in transit_detection!():44.801815168

(12) n = 12: Total time used in transit_detection!():42.965138617

Compared to the previous results, this shows much better parallization.


#### One planet test

Run
```
julia -p n para_one_planet.jl
```
to do the parallelization for the real data test.

Run time vs number of processors:

(1)

(2)

(3)

(4)

(5) n = 5, run time: 

(6) n = 6, run time: 605.452938523 seconds




 
[1]Aigran & Irwin 2004 (http://arxiv.org/abs/astro-ph/0401393)
