using PyCall
#using PyPlot
@pyimport kplr
#@pyimport matplotlib.pyplot as pl
#@pyimport bls

first_day_quarter = [54964.,55002.,55092.,55184.,55277.,55371.,55462.,55567.,
    55644.,55739.,55833.,55930.,56014.,56203.,56304.,]
time_slice_correction = (0.25+0.62*(5-1))/(86400)
N_good = 500
N_bad = 3

function getData()
    client = kplr.API()
    koi = client[:koi](97.01)
    println("period = ",koi[:koi_period])
    lightcurves = koi[:get_light_curves]()
    TIME = Float64[]
    FLUX = Float64[]

    n_curves = length(lightcurves)
    print ("length of lightcurves = ", n_curves, "\n")
    for i = 5:5        # loop over all light curves
        curve = lightcurves[i][:open]()     # get one curve, i.e. open one fits file
        data = curve[2][:data]              # get the data in curve[2]
        sap_quality = convert(Array{Float64,1},data[:field]("SAP_QUALITY"))
        sap_flux = convert(Array{Float64,1},data[:field]("SAP_FLUX"))
        time = convert(Array{Float64,1},data[:field]("TIME"))
        timecorr = convert(Array{Float64,1},data[:field]("TIMECORR"))
        bjd = time + 2454833.0                 # get BJD time
        jd = bjd - timecorr + time_slice_correction  # get JD time
        mjd = jd - 2400000.5                    # get MJD time
        #("mjd time =", mjd)
        good_index = Int64[]
        n_time = length(time)
        for i = 1:n_time
            if (is(sap_flux[i],NaN) || is(sap_flux[i],Inf) || is(time[i],NaN)|| is(time[i],Inf) || sap_quality[i]!=0 || in(floor(mjd[i]),first_day_quarter))
                # ith element is a bad data point (time = NaN or Inf, sap_quality !=0, sap_flux=NaN or Inf)
                # or it is on the first day of quarters, do nothing
            #elseif ()
                # ith element is near a known transits candidates , do nothing
            else
                # ith element is considered to be a good data point, add index i to the good_index array
                append!(good_index,[i])
            end
        end

        # identify usable light curve segments
        segmentIndex = getSegmentIndex(N_good, N_bad,good_index)       # segmentIndex is an array of arrays, its each element is an array of the indexes of good data points that can be considered as one segment

        # combine all segments together into one array
        nseg = length(segmentIndex)
        println("nseg = ", nseg)
        for i= 1:nseg
            # normalize each segment
            sap_flux[segmentIndex[i]]= normalize(sap_flux[segmentIndex[i]])
            # estimate each segment's noise level
            noise = estimateNoise(sap_flux[segmentIndex[i]])
	    println("Noise = ", noise)
            #sap_flux[segmentIndex[i]] -= noise

            append!(TIME,time[segmentIndex[i]])
            append!(FLUX,sap_flux[segmentIndex[i]])
        end
    end
    println("length of TIME = ", length(TIME))

    #plot everything
   #plotEverything(TIME,FLUX,"97.01")

   # use Aigrian & Irwin's method to search for the period, epoch and duration that gives the biggest SNR (or Q value)
   # input parameters are (TIME, FLUX, length_p,f_min, f_max), whre length_f is the length of the trial frequency 
    # array, f_min and f_max are the minimum and maximum trial frequency. f_min must >= 1/total time span.
   length_f=20
   f_min=1/5.
   f_max=1/4.6
   transit_detection!(TIME,FLUX,length_f,f_min, f_max);
end

function getSegmentIndex(Ngood::Int64, Nbad::Int64, index::Array)
    # Ngood: minimum number of good points in a row
    # Nbad:  number of bad points in a row makes the segment break
    n = length(index)
    segment = Any[]
    push!(segment,[index[1]])
    for i = 1:n-1
        Nseg = length(segment)         # the current length of segment[]
        if (Nseg == 0)
            push!(segment,[index[i+1]])
        else
            if(index[i+1] - index[i] < Nbad)
                append!(segment[Nseg],[index[i+1]])
            elseif(length(segment[Nseg]) < Ngood)
                pop!(segment)
            else
                push!(segment,[index[i+1]])
            end
        end
    end
    if(length(segment[end]) < Ngood)        # check the last element
        pop!(segment)
        end
    return segment
end

function normalize(x::Array)
    (x-median(x))/median(x)
end

function estimateNoise(FLUX::Array)
    median(abs(FLUX - median(FLUX)))
end

function transit_detection!(TIME::Array,FLUX::Array,length_f::Int64,trial_f_min::Float64, trial_f_max::Float64)
    @assert length(TIME) == length(FLUX)
    @assert typeof(length_f) == Int64
    tot = TIME[end]-TIME[1]        # total time span
    @assert trial_f_min >= 1/tot
    println("total time span = ",tot)

    # get the avage time gap between consecutive data points avg_t_step, it will be used as the trial duration step
    n_time_gaps = 0
    n = length(TIME)
    sum_time_gaps = 0
    for i=2:n
        if(TIME[i]-TIME[i-1]<0.5)     # significant time gaps are not counted
            sum_time_gaps += TIME[i] - TIME[i-1]
            n_time_gaps += 1
        end
    end
    #avg_t_step = sum_time_gaps/n_time_gaps
    avg_t_step = 0.06
    println("avg_t_step = ", avg_t_step)

    # get the trial period array from trial frequency
    p=zeros(Float64,length_f)
    length_p = length_f
    df = (trial_f_max - trial_f_min)/length_f   # trial frequency step
    println("df = ", df, "\n")
    for i=1:length_f
        p[i]= 1./(trial_f_max - (i-1)*df)
    end
	println("p = ",p)

    # get the trial period array p, trial duration array d, and trial epoch array e
    # calculate Q(p,d,e) for each p,d,e
    best_logQ =	-1000 
    (best_p,best_d,best_e) = (0.,0.,0.)
    logQ_p=zeros(length_p)
    for i=1:length_p
        # at a given trial period, do phase-folding
        trial_p = p[i]
        foldedTIME=[]       # the array to store folded time
        foldedFLUX =[]      # the array to store folded flux
        (foldedTIME,foldedFLUX) = phase_folding(TIME,FLUX,trial_p)
        # at this given trial period, get the trial duration array, the default minimum trial duration = avg_t_step;
        # max duration = trial period
	println("trial_p = ", trial_p)
        length_d = int(trial_p/avg_t_step)
        println("\n","length_d = ", length_d)
        d = zeros(length_d)
        for j=1:length_d
            d[j] = avg_t_step*j      # trial duration step = avg_t_step
            trial_d = d[j]
            # at a given trial duration, get the trial epoch array, the default minimum trial epoch = 0;
            # maximum trial epoch = tot - current trial duration, trial epoch step = avg_t_step,too
            length_e = length_d - j
	    #println("length_e = ", length_e)
            e = zeros(length_e)
            I_index = []                               # the array to store in-transit data points' index
            current_idx = 1                             # the index used in the searching for transit times
            for k=1:length_e
                e[k]= avg_t_step*(k-1)
                trial_e = e[k]
                I_index = get_in_transit_index!(foldedTIME,trial_e,trial_d,current_idx)
                logQ = calculate_logQ(I_index,foldedFLUX,trial_p,trial_d,trial_e)
		if(logQ>best_logQ)
			best_logQ = logQ
			(best_p,best_d,best_e) = (trial_p,trial_d,trial_e)
		end
            end
        end
	logQ_p[i] = best_logQ
	println("at trial_p = ", trial_p, ", best_logQ = ", best_logQ,"\n", "(best_p,best_d,best_e) = ", (best_p,best_d,best_e))
     end
     #plot(p,logQ_p)
     #fig_name = string("/gpfs/home/fxh140/Astro585/project/logQ_p_", fileName , ".pdf")
     #savefig(fig_name)
     #plotEverything(p,logQ_p,"5th_segment")
     return (best_p,best_d,best_e)
end

function calculate_logQ(I_index::Array,foldedFLUX::Array,trial_p::Float64,trial_d::Float64,trial_e::Float64)
    length_I_index = length(I_index)
    intransitFLUX = foldedFLUX[I_index]
    corr_flux = intransitFLUX - mean(intransitFLUX)    # mean corrected flux
    MAD = median(abs(corr_flux- median(corr_flux)))
    logQ = 2*log(abs(mean(corr_flux))) + log(length_I_index) - 2*log(MAD)
#	println("e =",trial_e)
    #println("logQ = ", logQ)
    return logQ
end

function get_in_transit_index!(foldedTIME::Array,trial_e::Float64,trial_d::Float64, current_idx::Int64)
    @assert current_idx < length(foldedTIME)
    # search for in-transit data points in foldedTIME[]
    trial_transit_end = trial_e + trial_d
    # use I_index to store the index of in-transit data points
    I_index = Int64[]
    while(foldedTIME[current_idx] < trial_e)
        current_idx+=1
    end     # after this while loop, foldTIME[current_idx] >= trial_e
    if(current_idx ==1 || foldedTIME[current_idx] >=  trial_e && foldedTIME[current_idx-1] < trial_e)
        while(foldedTIME[current_idx] < trial_transit_end)
            # step forward and get the index for the transit
            append!(I_index,[current_idx])
            current_idx += 1
        end # after this while loop, foldedTIME[current_idx] >= trial_transit_end
    else
        while(foldedTIME[current_idx] >=  trial_e)
            current_idx -= 1
        end     # after this while loop, foldedTIME[current_idx] < trial_e, so we need to add 1 to it.
        current_idx +=1             # then step forward again
        while(foldedTIME[current_idx] < trial_transit_end)
            append!(I_index,[current_idx])
            current_idx += 1
        end
    end
    # when this function finishes, current_idx = 1+ the index corresponding the last in-transit point
    return I_index
end

function phase_folding(TIME::Array, FLUX::Array, period::Float64)
    n = length(TIME)
    t0 = TIME[1]
    phase=zeros(Float64,n)
    for i=2:n
        t = TIME[i]
        phase[i] = mod((t-t0),period)   #phase's range is [t0,period]
    end
    p=sortperm(phase)
    foldedTIME = phase[p]   # sort the TIME array
    foldedFLUX = FLUX[p]    # sort the FLUX array accordingly
    return (foldedTIME, foldedFLUX)
    end

function plotEverything(TIME::Array, FLUX::Array,fileName::String)
    #pl.scatter(TIME,FLUX)
	plot(TIME,FLUX)
    x_min = minimum(TIME)
    x_max = maximum(TIME)
    x_min -= (x_max - x_min)*0.05
    x_max += (x_max - x_min)*0.05
    y_min = minimum(FLUX)
    y_max = maximum(FLUX)
    y_min -= (y_max - y_min)*0.05
    y_max += (y_max - y_min)*0.05
    pl.xlim([x_min,x_max])
    pl.ylim([y_min,y_max])
    fig_name = string("/gpfs/home/fxh140/Astro585/project/logQ_p_", fileName , ".pdf")
    #fig_name = string("/Users/feifeihuang/project/plots/all_segment", fileName , ".pdf")
    pl.savefig(fig_name)
    pl.show()
    return
end

getData()
