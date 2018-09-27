#!/usr/bin/python

####################################################################
# runSimulations.py
# Jaline Gerardin 2010
#
# Given a parameter file, the function runOneGroup simulates circuit
# response to inputs of various duration. Temporal dose responses and
# kintic filtering metrics (temporal ultrasensitivity score and
# trigger time) are calculated and outputed to file.
#
####################################################################

import commands
import os

HOMEDIR = ''

# runOneGroup reads in parameters from file and outputs a temporal
# dose response and kinetic filtering metrics for each parameter set
# that ran successfully.
#
# runOneGroup's 'program' argument is the code that simulates circuit
# response to input.
#
# required format for parameter file: 
# line 1: basal input, change in input, node logic for 3 nodes,
# topology ID number, starting parameter ID number
# subsequent lines: 26 kcat's and Km's for each line. Link regulation
# type (activator, inhibitor, absent) is coded in kcat values
# (positive, negative, 0 respectively)
#
# Temporal dose responses are outputed to rundir and kinetic filtering
# metrics are outputed to directory stem_[node A logic][node B
# logic][node C logic]
def runOneGroup(program, paramfile, rundir, stem):

    # open parameter file and read in settings in first line
    finparam = open(paramfile)
    settings = finparam.readline()
    [base_input, change, Aand, Band, Cand, topnum, circuit_index] = [n for n in settings.split()]
    circuit_index = int(circuit_index)

    # file for storing temporal dose response curves
    rawoutfile = rundir + 'raw/data_' + str(topnum)
    # file for storing temporal dose reponse metrics
    metdir = os.path.join(HOMEDIR, stem + '_' + Aand + Band + Cand)
    outfile = os.path.join(metdir, 'data_%d' % topnum)

    # start indexing parameters at circuit_index
    param = circuit_index

    # input durations to test
    intimes = [0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 800, 1000, 2000, 3000, 5000, 6000, 8000, 10000, 20000, 50000]

    # for each parameter set in parameter file
    for thisparam in finparam.readlines() :
        print 'running top ' + str(topnum) + ' param ' + str(param)

        data = [] # stores maximum output amplitudes
        times = [] # stores input durations
        params = ' '.join(thisparam.split())

        # try to initialize the circuit
        initialization = runOneTime(program, intimes[0], base_input, change, Aand, Band, Cand, params, '')
        # if circuit did not reach steady state
        if(initialization < 0) : 
            param += 1
            continue # continue to next parameter set

        # apply shortest input duration
        diff0 = runOneTime(program, intimes[0], base_input, change, Aand, Band, Cand, params, initialization)
        # if error in simulation
        if diff0 < 0 : 
            param += 1
            continue # continue to next parameter set

        # apply longest input duration
        difflast = runOneTime(program, intimes[len(intimes)-1], base_input, change, Aand, Band, Cand, params, initialization)
        # if amplitude is too small or amplitude is identical for shortest and longest inputs
        if(difflast < 10**-30 or abs(difflast - diff0) < 10**-30) : 
            param += 1
            continue # continue to next parameter set

        # store shortest input duration and output amplitude
        times.append(intimes[0])
        data.append(diff0)

        # check that all input duration simulations finished without error
        allok = 1

        # test the remaining input durations
        for i in range(1, len(intimes)-1) :
            next = runOneTime(program, intimes[i], base_input, change, Aand, Band, Cand, params, initialization)
            if next < 0 : # if error in simulations, break out of loop
                allok = 0
                break
            # store output amplitude and input duration
            data.append(next)
            times.append(intimes[i])
            # starting with 6th input duration, check if output amplitude is equivalent to amplitude for longest input duration
            if i > 4 and abs(next - difflast) < 10**-15 : # if so, skip remaining input durations and break out of loop
                break
        
        # if simulation error occurred, continue to next parameter set
        if allok == 0 :
            param += 1
            continue

        # if all input durations were run, append data for longest duration
        if len(times) == len(intimes) - 1 :
            data.append(difflast)
            times.append(intimes[len(intimes)-1])

        # record maximum output amplitudes to file
        with open(rawoutfile, 'a') as fout :
            fout.write(str(topnum) + '\t' + str(param) + '\t')
            for i in range(len(data)):
                fout.write(str(data[i])+'\t')
            fout.write('\n')
        
        # calculate temporal dose response metrics and write to file if calculations were error-free
        ultradata = getultra(times, data)
        if len(str(ultradata)) > 0 and ultradata != -1:
            with open(outfile, 'a') as fout :
                fout.write(str(topnum) + '\t' + str(param) + '\t' + str(ultradata) + '\t' + str(max(data) - min(data)) + '\n')

        # increment parameter tracker
        param += 1
    finparam.close()


# run simulation of one parameter set, one input duration
def runOneTime(program, intime, base_input, change, Aand, Band, Cand, params, initialization) :

    line = commands.getoutput('%(program)s %(intime)s %(base_input)s %(change)s %(Aand)s %(Band)s %(Cand)s %(params)s %(initialization)s' % vars())

    try :
        data = [float(n) for n in line.split()]
    except ValueError:
        print line
        return -1

    if(data[0] == -1) : # if circuit could not reach steady state during initialization, return error
        return -1
    if initialization == '' : # if in initialization phase, return initial steady state
        return line
    if data[10] == -1 or len(data) < 24 : # if something went wrong during simulation, return error
        return -1

    return data[14] # otherwise return maximum output amplitude


# calculate temporal dose response metrics
def getultra(times, data) :

    # program used to calculate temporal dose response metrics
    program = 'doseResponseMetrics'

    # calculate metrics
    numpoints = len(times)
    mytimes = ' '.join(['%.2f'%x for x in times])
    mydata = ' '.join(['%.15f'%x for x in data])
    TUmetrics = commands.getoutput(HOMEDIR + program + ' ' + str(numpoints) + ' ' + mydata + ' ' + mytimes)

    # if error in calculating metrics, return -1
    if ('nan' in TUmetrics or 'inf' in TUmetrics) :
        return -1

    # return metrics
    return TUmetrics

