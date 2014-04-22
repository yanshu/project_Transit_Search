
# Transit Search

This is a Astro585 project that prepares some simple code to search for transiting planets candidates in Kepler data using transit method.

### Introduction

Kepler is a space observatory lauched by NASA in 2009, its purpose is to search for planets. It uses a photometer that records the brightness of stars. [1] A planet is defined as an astronomical object that travels arount a star. When  planet goes across the face of a star, the brightness of the star would decrease, this is called a transit, the duration of the decrease is the time the planet spends on moving across the surface, the period of a transit is equal to the orbital period of the planet, and the epoch of the transit is the start time of the first transit detected by Kepler. 

So we can use the kepler data to search for planet candidates by searching for transits: given the light curves (more specifically the "sap_flux" and "time" data from the light curves), find out the period, duration, and epoch of the transit, and the corresponding signal-to-noise ratio(SNR).

### Search algorithm

There are different transit search algorithms. This project uses Aigran & Irwin's transit search algorithm[2]. The basic idea is: (1) divide the data points in the flux vs time plot into two bins, one in-transit bin and one out-of-transit bin, (2) loop over different trial periods (p), trial transit durations (d) and trial transit epochs (e), find out the (p, d, e) that maximizes the Q value for the data points in the in-transit bin. See the definition of Q in the Aigran & Irwin paper, it is basically the square of the signal to noise ratio (SNR) [2].

###Steps

The main steps of this transit search project:

(1) Given the raw time array and flux array, strip out bad data points; ( Good data points selection)

(2) Generate arrays of trial periods, durations, and epochs;

(3) For each trial periods, do phase folding;

(4) Use the phase folded flux, use the Aigran & Irwin search algorithm, find out the best estimated p, d, e.

functions.jl and para_functions.jl contain all the functions for this project. The former is the serial vesion, and the latter is the parallel version. The core function is: transit_detection!(), whose input parameters are: (TIME,FLUX,length_f,f_min, f_max). TIME and FLUX are the preprocessed time array and flux array, length_f is the number of trial frequencies, f_min, f_max are the minimum and maximum trial frequencies. The output of transit_detection!() is ((best_p,best_d,best_e,best_logQ,logQ_p), i.e. the best estimated period, transit duration and epoch, the corresponding logQ value and the array that stores each logQ value for each trial period.

The main steps of transit_detection!() is: (1) For a given trial period value, do phase folding; (2) For each trial duration and trial epoch, identify the in-transit data points in this phased flux, and calculate the logQ for it; (3) Compare all logQ values, get the maximum logQ, then the corresponding (p, e, d) is the best estimated result.

#### Good data points selection

The defintion for bad data points is: quality flag!=0, sap_flux or pdf_flux = NaNs, Infs, first day of each quarter.

#### Phase folding

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

(1) n = 1: Total time used in transit_detection!(): 235.974000392

(2) n = 2: Total time used in transit_detection!(): 123.928522869

(3) n = 3: Total time used in transit_detection!(): 128.286416456

(4) n = 4: Total time used in transit_detection!(): 99.218811207

(5) n = 5: Total time used in transit_detection!(): 107.634166859 

(6) n = 6: Total time used in transit_detection!(): 71.39507179

(7) n = 7: Total time used in transit_detection!(): 71.575478587

(8) n = 8: Total time used in transit_detection!(): 57.220656825

(9) n = 9: Total time used in transit_detection!(): 57.918904511

(10) n = 10: Total time used in transit_detection!(): 50.570660036

(11) n = 11: Total time used in transit_detection!(): 48.488129586

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

(1) n = 1: Total time used in transit_detection!(): 215.752591133

(2) n = 2: Total time used in transit_detection!(): 123.095294461

(3) n = 3: Total time used in transit_detection!(): 90.721392053

(4) n = 4: Total time used in transit_detection!(): 75.883206517

(5) n = 5: Total time used in transit_detection!(): 68.422990043

(6) n = 6: Total time used in transit_detection!(): 64.156633272

(7) n = 7: Total time used in transit_detection!(): 54.97174648

(8) n = 8: Total time used in transit_detection!(): 51.917523137

(9) n = 9: Total time used in transit_detection!(): 46.779422036

(10) n = 10: Total time used in transit_detection!(): 49.136618533

(11) n = 11: Total time used in transit_detection!(): 44.801815168

(12) n = 12: Total time used in transit_detection!(): 42.965138617

Compared to the previous results, this shows much better parallization.


#### One planet test

Run
```
julia -p n para_one_planet.jl
```
to do the parallelization for the same one plant data test. Vary n from n= 1, 2, ... 12, the result is always the same: period = 4.890738813735692 days, duration = 4.536 hours. The run time vs number of processors:

(1) n = 1, run time: 2418.599992104 seconds

(2) n = 2, run time: 1429.745141357 seconds

(3) n = 3, run time: 999.511613143 seconds

(4) n = 4, run time: 779.576668775 seconds

(5) n = 5, run time: 677.916969961 seconds

(6) n = 6, run time: 605.452938523 seconds

(7) n = 7, run time: 576.667202719 seconds

(8) n = 8, run time: 518.513130904 seconds

(9) n = 9, run time: 501.187991569 seconds 

(10) n = 10, run time: 458.875286608 seconds

(11) n = 11, run time: 468.668311456 seconds

(12) n = 12, run time: 437.128084147 seconds

We see a significant decrease of process time as number of processors increase. This means the parallelized code works well.



[1] http://en.wikipedia.org/wiki/Kepler_%28spacecraft%29 
[2]Aigran & Irwin 2004 (http://arxiv.org/abs/astro-ph/0401393)
