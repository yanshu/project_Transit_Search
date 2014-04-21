# This script first gets the kplr light curve data of one confirmed planet, then does good data points selection, then identify
# usable light curve segments using these good data points, and in the end runs the transit_detection!() function on the segments,
# the output is the best estimated period and duration for the input light curve. The script name para_one_planet indicates
# there is parallelization here, the parallelization exists in the transit_detection!() function only (see it in para_functions.jl), 
# other parts are not parallelized.

include("para_functions.jl")
using PyCall
@pyimport kplr
#@pyimport matplotlib.pyplot as pl

#############################################################################################################################################
##### get kepler lightcurve data

println("Getting lightcurves of a planet candidate : KOI 97.01")

client = kplr.API()
koi = client[:koi](97.01)		# in this test, we use the transiting planet candidate:  KOI ID =97.01 
println("True period = ",koi[:koi_period], " days, true transit duration = ", koi[:koi_duration], " hours")	#true period and duration
lightcurves = koi[:get_light_curves]()  # get lightcurves
n_curves = length(lightcurves)		# the number of lightcurves

#############################################################################################################################################
##### select good data points

println("Entering good data points selection...")

TIME = Float64[]	# the array to store the time information that'll be used as inputs for transit detection 
FLUX = Float64[]	# the array to store the flux information that'll be used as inputs for transit detection 

for i = 5:5        # select only the 5th lightcurve for test (can be changed to loop over all light curves: for i=1:n_curves)
    curve = lightcurves[i][:open]()     	# open lightcurve[5]	
    data = curve[2][:data]              	# get data from curve 
    sap_quality = convert(Array{Float64,1},data[:field]("SAP_QUALITY"))  # get the sap_quality information
    sap_flux = convert(Array{Float64,1},data[:field]("SAP_FLUX"))	# get the flux information, convert it into an array sap_flux[]
    time = convert(Array{Float64,1},data[:field]("TIME"))	# get the time information, convert it into an array time[]
    timecorr = convert(Array{Float64,1},data[:field]("TIMECORR"))	# get the time correction information, it'll be used for good data selection
    first_day_quarter = [54964.,55002.,55092.,55184.,55277.,55371.,55462.,55567.,55644.,55739.,55833.,55930.,56014.,56203.,56304.]
					 # the first day of each quarter (format: MJD time), any data points in these days is considered bad points
    
    ### convert time to MJD format in order to compare with fisrt_day_quarter[]
    bjd = time + 2454833.0                 # get BJD time
    time_slice_correction = (0.25+0.62*(5-1))/(86400)	# time_slice_correction is used to convert BJD time to JD time
    jd = bjd - timecorr + time_slice_correction  # get JD time
    mjd = jd - 2400000.5                    # get MJD time
    
    ### get rid of bad points, store the indices of good data points
    good_points_index = Int64[]		# an array to store the indices of good data points
    n_time = length(time)
    for i = 1:n_time		# loop over all data points, store good data points' indices into good_points_index[]
        if (is(sap_flux[i],NaN) || is(sap_flux[i],Inf) || is(time[i],NaN)|| is(time[i],Inf) || sap_quality[i]!=0 || in(floor(mjd[i]),first_day_quarter))
            # it is a bad data point (bad points definition: time = NaN or Inf, sap_quality !=0, sap_flux=NaN or Inf or it is on the first day of quarters
        else			# it is a good data point, add index i to the good_points_index array
            append!(good_points_index,[i])
        end
    end
    println("Finished data points selection")

#############################################################################################################################################
##### identify usable light curve segments using the good data points we just got [ see segments' definition in the code]
    
    println("Identifying and combining usable light curve segments ... ")
    ### identify usable light curve segments [ see function getSegmentIndex() in para_functions.jl]
    N_good = 500	# to define one segment, first, we require each segment to have at least a minimum number of good points: N_good
    N_bad = 3		# second, we define segment breaks whenever several bad points are in a row, the number of bad points = N_bad
   			# where bad points definition is: time = NaN, Inf or the first day of each quarter, sap_quality!=0, sap_flux=NaN,Inf
    segmentIndex = getSegmentIndex(N_good, N_bad,good_points_index)       # the getSegmentIndex() function returns an array of arrays, 
    			#segmentIndex is an array of arrays, its length is the number of usable segments, its each element is an array of the 
			#indices of good data points in one segment

    ### normalize each segment and combine all segments, get the TIME[] and FLUX[] array
    nseg = length(segmentIndex)
       for i= 1:4      # for this test we only use 4 segments out of 11 usable segments(can be changed to loop over all segments: for i= 1:nseg)
	median_value = median(sap_flux[segmentIndex[i]])
        sap_flux[segmentIndex[i]]= (sap_flux[segmentIndex[i]]-median_value)/median_value		# normalize each segment
        append!(TIME,time[segmentIndex[i]])
        append!(FLUX,sap_flux[segmentIndex[i]])		# combine good segments, get TIME[] and FLUX[] arrays, the inputs for transit_detection!()
    end
end
println("Finished combining usable segments, length of TIME = ", length(TIME))

#plot FLUX vs TIME, if you need
#plotEverything(TIME,FLUX,"97.01")

#############################################################################################################################################
# use Aigrian & Irwin's method to search for the period, epoch and duration that gives the biggest log10(Q) value (see more about this in README)
# see transit_detection!() in para_function.jl, its input parameters are (TIME, FLUX, length_p,f_min, f_max), where length_f is the length of the 
# trial frequency array, f_min and f_max are the minimum and maximum trial frequency. f_min must >= 1/total time span.
# the output parameter is (best_period, best_duration, best_epoch, best_logQ).

(length_f,f_min, f_max) = (20,1/5.,1/4.7)	# here for a quick test, use 20 as the number of trial frequencies, f_min = 1/5. and f_max = 1/4.7,
						# since we use a confirmed planet and we know the true frequency is 1/4.88.  But normally, we wouldn't
						# know anything about its true frequency in advance, so he trial frequency would have a large range
						# and a big size, note: f_min must >= 1/total time span, it will return an error if this condition is
						# not satisfied.
tic()
result = transit_detection!(TIME,FLUX,length_f,f_min, f_max,func);
println("Time used to run transit_detection:",toq())

println("The best estimated period and duration is: period = ",result[1], " days, duration = ", result[2]*24, " hours") 
println("The true period and duration is : period = ",koi[:koi_period], " days, duration = ", koi[:koi_duration], " hours")


