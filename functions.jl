# This script contains five functions:
# four core functions for transit detection:
# (1) transit_detection!()
# (2) phase_folding!()
# (3) calculate_logQ!()
# (4) get_in_transit_index!()
# and two other functions for getting usable segment indices and plotting flux vs time:
# (5) plot_Everything()
# (6) getSetmentIndex()

function transit_detection!(TIME::Array,FLUX::Array,length_f::Int64,trial_f_min::Float64,trial_f_max::Float64)
    println("Entering transit_detection function")
    @assert length(TIME) == length(FLUX)
    @assert typeof(length_f) == Int64
    tot = TIME[end]-TIME[1]        # total time span, the trial_f_min must be larger than 1/tot
    @assert trial_f_min >= 1/tot

    # get the noise for FLUX array
    median_value = median(FLUX)
    noise = median(abs(FLUX - median_value))	# noise is defined as MAD, median absolute deviation
    
    # get the avage time gap between consecutive data points avg_t_step, it will be used as the trial duration and epoch step
    tic()
    n = length(TIME)
    sum_time_gaps = 0
    n_time_gaps = 0
    for i=2:n
        if(TIME[i]-TIME[i-1]<5)     # significant time gaps are not counted
            sum_time_gaps += TIME[i] - TIME[i-1]
            n_time_gaps += 1
        end
    end
    avg_t_step = sum_time_gaps/n_time_gaps	# use the normal definition for simulated_data test
    println("Time used in getting avg_t_step:")
    toc()
    #avg_t_step = 0.003 		# for the one_planet_test, use 0.003, a larger step to speed things up
    println("	trial duration and trial epoch step: avg_t_step = ", avg_t_step)

    # get the trial arrays p[], d[] and e[]
    tic()
    p=zeros(Float64,length_f)
    length_p = length_f
    df = (trial_f_max - trial_f_min)/length_f   # trial frequency step
    for i=1:length_p
        p[i]= 1./(trial_f_max - (i-1)*df)
    end
    max_length_d_and_e = int(p[end]/avg_t_step)
    d=zeros(Float64,max_length_d_and_e)
    e=zeros(Float64,max_length_d_and_e)
    for i=1:max_length_d_and_e	 # p[end]/avg_t_step is the largest possible length for the duration array, even though we only need a portion of them, because given a duration value d[i], the length of e we need is only length_d-i+1
        d[i] = avg_t_step*i      # trial duration and epoch step = avg_t_step
	e[i] = avg_t_step*(i-1)		# e[1] = 0.
    end
    println("Time used in getting trial arrays p, e and d:")
    toc()

    # calculate log10(Q)(p,d,e) for each trial period, trial duration, trial epoch

    tic()
    best_logQ =	-10000 		# best_logQ is used to store the largest logQ for all (p,d,e), initialized to be a very small value
    (best_p,best_d,best_e) = (0.,0.,0.)		# best_p, best_d, best_e are used to store the period, duration, epoch corresponding the best_logQ
    logQ_p=zeros(length_p)	# an array to store the logQ value for each trial period
    time_phase_folding = 0.	# record the time for doing phase folding
    time_get_in_transit_index = 0.	# record the time for get_in_transit_index
    time_calculate_logQ = 0.	# record the time for calculate_logQ
    for i=1:length_p		# loop through all trial periods
        # at a given trial period, do phase-folding, foldedTIME[] is the array to store folded time , foldedFLUX[] is the folded flux
        trial_p = p[i]
	tic()
        (foldedTIME,foldedFLUX) = phase_folding!(TIME,FLUX,trial_p,avg_t_step)
	time_phase_folding += toq()
        # at this given trial period, get the trial duration array, trial duration step = avg_t_step, the maximum trial duration = trial period
        length_d = int(trial_p/avg_t_step)
    	local_best_logQ = -10000 		# best logQ for the current trial period, initialized to be a very small value
	(local_best_p,local_best_d,local_best_e) = (0.,0.,0.)		# best logQ for the current trial period
        for j=1:length_d	# loop through all trial durations
            trial_d = d[j]
            # at this given trial duration, get the trial epoch array, the minimum trial epoch = 0, trial epoch step = avg_t_step,too, so length of e is:
            length_e = length_d-j+1
            #current_idx = 1                     # the index used in the searching for transit times
            for k=1:length_e
                trial_e = e[k]
		tic()
                I_index = get_in_transit_index!(foldedTIME,trial_e,trial_d,1)		# get the indices of all in-transit data points
		time_get_in_transit_index += toq()
		if length(I_index)==0			# if there is no in-transit data points, move to the next epoch
			current_idx+=1
			continue
		else				# if there is in-transit data points, calculate the logQ value for these points
			tic()
                	logQ = calculate_logQ!(I_index,foldedFLUX,trial_p,trial_d,trial_e,noise)
			time_calculate_logQ += toq()
			if(logQ>local_best_logQ)
				local_best_logQ = logQ
			        (local_best_p,local_best_d,local_best_e) = (trial_p,trial_d,trial_e)
			end
		end
            end
	end
	logQ_p[i] = local_best_logQ
	if(local_best_logQ>best_logQ)
                best_logQ = local_best_logQ
		(best_p,best_d,best_e) = (local_best_p,local_best_d,local_best_e)
    	end
    	println("  for trial_p = ", trial_p, ", best_logQ = ", local_best_logQ, ", (best_p,best_d,best_e) = ", (local_best_p,local_best_d,local_best_e),"\n")
    end
    
    println("Finishing transit detection, of all trial periods,durations and epochs, the best one is : (best_p, best_d, best_e) = ",(best_p, best_d, best_e), "\n")
    println("Time used in phase folding: ", time_phase_folding)
    #println("The time used in running the 3-level for loops", toq())
    println("Time used in get_in_transit_index: ", time_get_in_transit_index)
    println("Time used in calculate_logQ: ", time_calculate_logQ)

    return (best_p,best_d,best_e,best_logQ,logQ_p)
end

function calculate_logQ!(I_index::Array,foldedFLUX::Array,trial_p::Float64,trial_d::Float64,trial_e::Float64,noise::Float64)
    @assert length(I_index)>0
    length_I_index = length(I_index)
    intransitFLUX = foldedFLUX[I_index]			# flux array that is in transit
    logQ = 2*log10(abs(mean(intransitFLUX))) + log10(length_I_index) - 2*log10(noise)		# calulate log10(Q) 
    return logQ
end

function get_in_transit_index!(foldedTIME::Array,trial_e::Float64,trial_d::Float64, current_idx::Int64 =1)
    if current_idx > length(foldedTIME)
        return Int64[]
    end
    # search for in-transit data points in foldedTIME[]
    trial_transit_end = trial_e + trial_d
    # use I_index to store the index of in-transit data points
    I_index = Int64[]
    while(foldedTIME[current_idx] < trial_e)
        current_idx+=1
    end     						# after this while loop, foldTIME[current_idx] >= trial_e
    if(current_idx ==1 || foldedTIME[current_idx-1] < trial_e)      # it is the 1st element after the trial_e, i.e. the start of the transit, so
        while(current_idx <= length(foldedTIME) && foldedTIME[current_idx] < trial_transit_end)			 # we step forward and record the indices for the transit
            append!(I_index,[current_idx])
            current_idx += 1
        end             # after this while loop, foldedTIME[current_idx] >= trial_transit_end or foldedTIME[current_idx] is the last array element
    else        # it has some distance to the 1st element after the trial_e, so move backward to get to that element
        while(current_idx >0 && foldedTIME[current_idx] >=  trial_e)
            current_idx -= 1
        end     # after this while loop, foldedTIME[current_idx] < trial_e, so we need to add 1 to it
        current_idx +=1             # now, current_idx is the 1st element after trial_e, so we can step forward again, and add indices to I_index
        while(current_idx <= length(foldedTIME) && foldedTIME[current_idx] < trial_transit_end)
            append!(I_index,[current_idx])
            current_idx += 1
        end
    end
    # when this function finishes, current_idx = 1+ the index of the last in-transit point
    return I_index
end

function phase_folding!(TIME::Array, FLUX::Array, period::Float64, avg_t_step::Float64)
    @assert length(TIME) == length(FLUX)
    println("  Entering phase_folding function")
    n_time = length(TIME)
    println("  before phase folding, length of FLUX array = ", n_time)
    t0 = TIME[1]
    phase=zeros(Float64,n_time)
    for i=2:n_time
        t = TIME[i]
        phase[i] = mod((t-t0),period)   #phase's range is [0,period]
    end
    p=sortperm(phase)
    foldedTIME = phase[p]   # sort the TIME array
    foldedFLUX = FLUX[p]    # sort the FLUX array accordingly
	
    # now the foldedFLUX array has the same data points as the input FLUX array, a lot of data points will be very close in time,
    # we want to combine data points that are very near to each other as one new data point
    # this time window is set to: avg_t_step
    mergedLength = 1+ifloor(foldedTIME[end]/avg_t_step)  # number of bins for this new merged flux array
    mergedTIME =zeros(mergedLength)
    mergedFLUX =zeros(mergedLength)
    j = 1
    good_index = Int64[]	# some of bins may have empty data points, i.e. there isn't any data points inside that bin, that's why we need an array good_index[] to record the non zero flux's indices
    for index = 1:mergedLength
        length_one_merged = 0		# the number of merged data points that are within the time window
        while j <= n_time && index == 1+ifloor(foldedTIME[j]/avg_t_step)
                mergedTIME[index] += foldedTIME[j]
                mergedFLUX[index] += foldedFLUX[j]
                j+=1
                length_one_merged += 1
        end
	if length_one_merged >0			# when this bin is not empty, store the index in good_index[]
       		append!(good_index,[index])
		mergedTIME[index] /= length_one_merged		# the average time of the merged point
        	mergedFLUX[index] /= length_one_merged		# the average flux value of the merged point
	end
    end
    goodMergedTIME = Float64[]		# the final flux and time array it returns
    goodMergedFLUX = Float64[] 
    for i=1:length(good_index)
	append!(goodMergedTIME,[mergedTIME[good_index[i]]])
	append!(goodMergedFLUX,[mergedFLUX[good_index[i]]])
    end
    println("  Exiting phase_folding, after phase folding, length of FLUX array = ", length(goodMergedFLUX))
    return (goodMergedTIME, goodMergedFLUX)
end

function getSegmentIndex(Ngood::Int64, Nbad::Int64,index::Array)
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

function plotEverything(TIME::Array, FLUX::Array,fileName::String)
    pl.scatter(TIME,FLUX)
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
    fig_name = string("logQ_p_", fileName , ".png")
    pl.savefig(fig_name)
    pl.show()
    return
end
