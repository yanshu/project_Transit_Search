include("functions.jl")
using PyCall
@pyimport kplr
@pyimport matplotlib.pyplot as pl


# download lightcurve data
client = kplr.API()
println("Downloading lightcurves of planet candidate : KOI 97.01")
koi = client[:koi](97.01)		# in this test, we use the transiting planet candidate:  KOI ID =97.01 
println("True period = ",koi[:koi_period], " days, true transit duration = ", koi[:koi_duration], " hours")	# the true period and duration
lightcurves = koi[:get_light_curves]()
n_curves = length(lightcurves)

TIME = Float64[]	# the array to store the final time information for transit detection 
FLUX = Float64[]	# the array to store the final flux information for transit detection 

println("Entering good data points selection...")

for i = 5:5        # select only the 5th lightcurve for test (can be changed to loop through all light curves: for i=1:n_curves)
    curve = lightcurves[i][:open]()     	# open lightcurve[5]	
    data = curve[2][:data]              	# get data from curve 
    sap_quality = convert(Array{Float64,1},data[:field]("SAP_QUALITY"))  # get the sap_quality information
    sap_flux = convert(Array{Float64,1},data[:field]("SAP_FLUX"))	# get the flux information, convert it into an array sap_flux[]
    time = convert(Array{Float64,1},data[:field]("TIME"))	# get the time information, convert it into an array time[]
    timecorr = convert(Array{Float64,1},data[:field]("TIMECORR"))	# get the time correction information, it'll be used for good data selection
    
    # get rid of bad points, store the indices of good data points into the array good_points_index
    N_good = 500	# to define one segment, first, we require segments to have at least a minimum number of good points: N_good
    N_bad = 3		# second, we define segment breaks whenever several bad points are in a row, number of bad points: N_bad
   			# where bad points definition is:  time = NaN or Inf or the first day of each quarter, sap_quality !=0, sap_flux=NaN or Inf
    first_day_quarter = [54964.,55002.,55092.,55184.,55277.,55371.,55462.,55567.,55644.,55739.,55833.,55930.,56014.,56203.,56304.] # the first day of each quarter (format: MJD time), any data points in these days is considered bad data
    # convert time to MJD format in order to compare with fisrt_day_quarter[]
    bjd = time + 2454833.0                 # get BJD time
    time_slice_correction = (0.25+0.62*(5-1))/(86400)	# time_slice_correction is used to convert BJD time to JD time
    jd = bjd - timecorr + time_slice_correction  # get JD time
    mjd = jd - 2400000.5                    # get MJD time
    
    good_points_index = Int64[]		# an array to store the indices of good data points
    n_time = length(time)
    for i = 1:n_time
        if (is(sap_flux[i],NaN) || is(sap_flux[i],Inf) || is(time[i],NaN)|| is(time[i],Inf) || sap_quality[i]!=0 || in(floor(mjd[i]),first_day_quarter))
            # it is a bad data point (bad points definition: time = NaN or Inf, sap_quality !=0, sap_flux=NaN or Inf or it is on the first day of quarters
        else			# it is a good data point, add index i to the good_points_index array
            append!(good_points_index,[i])
        end
    end

   # identify usable light curve segments, using the good_points_index 
    segmentIndex = getSegmentIndex(N_good, N_bad,good_points_index)       # segmentIndex is an array of arrays, its each element is an array of the indices of good data points in one segment

    # normalize each segment and combine all segments together into one array
    nseg = length(segmentIndex)
       for i= 1:4      # for this test we only use 4 segments out of 11 usable segments(can be changed to loop through all segments: for i= 1:nseg)
	median_value = median(sap_flux[segmentIndex[i]])
        sap_flux[segmentIndex[i]]= (sap_flux[segmentIndex[i]]-median_value)/median_value		# normalize each segment
        append!(TIME,time[segmentIndex[i]])
        append!(FLUX,sap_flux[segmentIndex[i]])
    end
end
println("Finished data points selection, number of good data points = ", length(TIME))

#plot FLUX vs TIME
#plotEverything(TIME,FLUX,"97.01")

# use Aigrian & Irwin's method to search for the period, epoch and duration that gives the biggest log10(Q) value (see more about this in README)
# the input parameters are (TIME, FLUX, length_p,f_min, f_max), where length_f is the length of the trial frequency array,
# f_min and f_max are the minimum and maximum trial frequency. f_min must >= 1/total time span.
# the output paramete is (best_period, best_duration, best_epoch, best_logQ, logQ_p). (logQ_p is an array of logQ values for each trial period)

(length_f,f_min, f_max) = (20,1/5.,1/4.7)

result = transit_detection!(TIME,FLUX,length_f,f_min, f_max);

println("The best estimated period and duration for this planet is: period = ",result[1], " days, duration = ", result[2]*24, " hours") 
println("The true period and duration is : period = ",koi[:koi_period], " days, duration = ", koi[:koi_duration], " hours")

