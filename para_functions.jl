# This script contains five functions:
# four core functions for transit detection:
# (1) transit_detection!()
# (2) phase_folding!()
# (3) func()
# (4) get_in_transit_index_and_calc_logQ!()
# (5) calculate_logQ!()
# (6) get_in_transit_index!()
# and one function for getting usable segment indices:
# (7) getSetmentIndex()

#Most of the functions here are the same with serialized functions. The biggest change is in transit_detection!() where 
#the original three-layer for loop is modified and two new functions are added to make it easy to parallelize. 
#The parallelization method is using distributed arrays. The index [i,j] is distributed over different processors.

function transit_detection!(TIME::Array,FLUX::Array,length_f::Int64,trial_f_min::Float64,trial_f_max::Float64,func::Function)
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
    println("Time used in getting avg_t_step:",toq())
    #avg_t_step = 0.003 		# for the one_planet_test, use 0.003, a larger step to speed things up
    println("   trial duration and trial epoch step: avg_t_step = ", avg_t_step)

    # get the trial arrays p[], d[] and e[]
    tic()
    p=zeros(Float64,length_f)
    length_p = length_f
    df = (trial_f_max - trial_f_min)/length_f   # trial frequency step
    for i=1:length_p
        p[i]= 1./(trial_f_max - (i-1)*df)
    end
    max_length_d = max_length_e = int(p[end]/avg_t_step)
    d=zeros(Float64,max_length_d)
    e=zeros(Float64,max_length_e)
    for i=1:max_length_d	 # p[end]/avg_t_step is the largest possible length for the duration array, even though we only need a portion of them, because given a duration value d[i], the length of e we need is only length_d-i+1
        d[i] = avg_t_step*i      # trial duration and epoch step = avg_t_step
	e[i] = avg_t_step*(i-1)		# e[1] = 0.
    end
    println("Time used in getting trial arrays p, e and d: ",toq())

    # calculate log10(Q)(p,d,e) for each trial period, trial duration, trial epoch
    logQ_all = Array(Any,length_p,max_length_d)
    (best_logQ,best_p,best_d,best_e) = (-10000.,0.,0.,0.)		# best_logQ: the largest logQ for all (p,d,e), initialized to be a very small value, best_p, best_d, best_e are the corresponding period, duration and epoch 
    
    time_phase_folding = 0.	# record the time for doing phase folding
    time_calc_logQ = 0.		# time for calc logQ

    tic()
    println("  Entering phase folding ...")
    foldedTIME = Array(Array,length_p)
    foldedFLUX = Array(Array,length_p)
    for i=1:length_p		# loop through all trial periods
        (foldedTIME[i],foldedFLUX[i]) = phase_folding!(TIME,FLUX,p[i],avg_t_step)
    end
    time_phase_folding =  toq()
    println("  Finished phase folding.")

    tic()
    println("  Entering calculation of logQ ... ")
    idx = distribute([ (i,j) for i in 1:length_p, j in 1:max_length_d])
    proclist = procs(idx)
    refs = [@spawnat proclist[proc] func(localpart(idx),foldedTIME,foldedFLUX,p,d,e,avg_t_step,noise) for proc in 1:length(proclist)]
    result = Any[]
    for proc in 1:length(proclist)
 	result = append!(result,fetch(refs[proc]))
        println("size of array from worker ",proclist[proc]," = ",length(fetch(refs[proc])))
    end
    println("   Finishing calculation of logQ")
    time_calc_logQ = toq()

    tic()
    for n=1:length(result)
	if(result[n][4]>best_logQ)
       		(best_p,best_d,best_e,best_logQ) = result[n] 
       	end
    end
    println("   Time spent on finding the maximum logQ = ", toq())
    println("Finishing transit detection, of all trial periods,durations and epochs, the best one is : (best_p, best_d, best_e) = ",(best_p, best_d, best_e), "\n")
    println("Time used in phase folding: ", time_phase_folding)
    println("Time used in calculation of logQ: ", time_calc_logQ)
    return (best_p,best_d,best_e,best_logQ)
end

@everywhere function get_in_transit_index_and_calc_logQ!(foldedTIME::Array,foldedFLUX::Array,p::Array,d::Array,e::Array,i::Int64,j::Int64,avg_t_step::Float64,noise::Float64)
	trial_p = p[i]
	trial_d = d[j]
	j_max = int(trial_p/avg_t_step)
	k_max = j_max-j+1
	(local_best_logQ,local_best_e) = (-10000.,e[1]) 
	if(j <= j_max)
		current_idx = 1                     # the index used in the searching for transit times
        	for k=1:k_max
        	    trial_e = e[k]
        	    I_index = get_in_transit_index!(foldedTIME,trial_e,trial_d,current_idx)         # get the indices of all in-transit data points
        	    if length(I_index)==0                   # if there is no in-transit data points, move to the next epoch
        	            current_idx+=1
        	            continue
        	    else                            # if there is in-transit data points, calculate the logQ value for these points
        	            logQ = calculate_logQ!(I_index,foldedFLUX,trial_p,trial_d,trial_e,noise)
			    if(logQ>local_best_logQ)                                   
				(local_best_e,local_best_logQ) = (trial_e,logQ)
                            end
        	    end
        	end
	end
	return (trial_p, trial_d, local_best_e,local_best_logQ)
end

@everywhere function calculate_logQ!(I_index::Array,foldedFLUX::Array,trial_p::Float64,trial_d::Float64,trial_e::Float64,noise::Float64)
    @assert length(I_index)>0
    length_I_index = length(I_index)
    intransitFLUX = foldedFLUX[I_index]			# flux array that is in transit
    logQ = 2*log10(abs(mean(intransitFLUX))) + log10(length_I_index) - 2*log10(noise)		# calulate log10(Q) 
    return logQ
end

@everywhere function get_in_transit_index!(foldedTIME::Array,trial_e::Float64,trial_d::Float64, current_idx::Int64)
    # search for in-transit data points in foldedTIME[]
    trial_transit_end = trial_e + trial_d
    # use I_index to store the index of in-transit data points
    I_index = Int64[]
    if (foldedTIME[end] < trial_e || foldedTIME[1] >= trial_transit_end)
        return Int64[]
    else
        I_index = Int64[]                           # use I_index to store the indices of in-transit data points
        while(current_idx < length(foldedTIME) && foldedTIME[current_idx] < trial_e)
            current_idx+=1               # current_idx is before the start of the trial transit,move forward
        end                              # until foldedTIME[current_idx] >= trial_e
        if(current_idx ==1 || foldedTIME[current_idx-1] < trial_e)          # it is the 1st element after the start of the transit, so:
            while(current_idx < length(foldedTIME) && foldedTIME[current_idx] < trial_transit_end) # step forward and record the indices
                append!(I_index,[current_idx])
                current_idx += 1
            end             # after this while loop, foldedTIME[current_idx] >= trial_transit_end
        else                # there're several points to the start of the trial transit, so move backward until it is before the start of the transit
            while(current_idx >0 && foldedTIME[current_idx] >=  trial_e)
                current_idx -= 1
            end              # after this while loop, foldedTIME[current_idx] < trial_e, so we need to add 1 to current_idx
            current_idx +=1             # now, current_idx is the 1st element after trial_e, so step forward and add indices to I_index
            while(current_idx < length(foldedTIME) && foldedTIME[current_idx] < trial_transit_end)
                append!(I_index,[current_idx])
                current_idx += 1
            end
    	end
    end
    # when this function finishes, current_idx = 1+ the index of the last in-transit point
    return I_index
end

@everywhere function func(idx::Array,foldedTIME::Array,foldedFLUX::Array,p::Array,d::Array,e::Array,avg_t_step::Float64,noise::Float64)
	result = Array((Float64,Float64,Float64,Float64),length(idx))
	println("entered func()")
	for l=1:length(idx)
		i = idx[l][1]
		j = idx[l][2]
		result[l] = get_in_transit_index_and_calc_logQ!(foldedTIME[i],foldedFLUX[i],p,d,e,i,j,avg_t_step,noise)
	end
	return result
end


function phase_folding!(TIME::Array, FLUX::Array, period::Float64, avg_t_step::Float64)
    @assert length(TIME) == length(FLUX)
    #println("  Entering phase_folding function")
    n_time = length(TIME)
    #println("  before phase folding, length of FLUX array = ", n_time)
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
    goodMergedTIME = Float64[]          # the final flux and time array it returns
    goodMergedFLUX = Float64[]
    for i=1:length(good_index)
        append!(goodMergedTIME,[mergedTIME[good_index[i]]])
        append!(goodMergedFLUX,[mergedFLUX[good_index[i]]])
    end
    #println("  Exiting phase_folding, after phase folding, length of FLUX array = ", length(goodMergedFLUX))
    return (goodMergedTIME, goodMergedFLUX)
end

function getSimulatedArray(n::Int64,p::Float64,d::Float64)
	TIME = zeros(n)
	for i=1:n
		TIME[i] += i
	end
	FLUX = randn(n)*0.1
	e = [50+p*i for i=0:ifloor((n-50)/p)]		# all the starting time of transits 
	for e_idx = 1:length(e)
		for i=e[e_idx]:e[e_idx]+d
			FLUX[i] += -1.+rand()*0.1
		end
	end
	return (TIME,FLUX)
end

