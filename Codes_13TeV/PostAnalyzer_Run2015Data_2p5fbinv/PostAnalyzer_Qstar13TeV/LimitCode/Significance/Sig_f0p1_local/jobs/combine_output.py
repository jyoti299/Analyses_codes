#!/usr/bin/env python

import os, sys, re, glob
from ROOT import TMath
from array import array


#######################################################
def getQuantiles(limits, median, onesigma, twosigma):

    nit = len(limits)
    if nit==0:
        return

    # sort the list with limits
    limits = sorted(limits)
    limits_array = array('d',limits)

    # median for the expected limit
    median[0] = TMath.Median(nit, limits_array)

    # quantiles for the expected limit bands
    prob = array('d',[]) # array with quantile boundaries
    prob.append(0.021)
    prob.append(0.159)
    prob.append(0.841)
    prob.append(0.979)

    # array for the results

    quantiles = array('d',[0., 0., 0., 0.])
    TMath.Quantiles(nit, 4, limits_array, quantiles, prob) # evaluate quantiles

    onesigma[0] = quantiles[1]
    onesigma[1] = quantiles[2]
    twosigma[0] = quantiles[0]
    twosigma[1] = quantiles[3]

    return
#######################################################


def main():

    res_type = sys.argv[1]
    masses = range(int(sys.argv[2]),int(sys.argv[3])+int(sys.argv[4]),int(sys.argv[4]))

    for mass in masses:
        # hadd .root files
        os.system('hadd ' + 'stats_' + str(mass) + '_' + res_type + '.root stats_' + str(mass) + '_' + res_type + '_*.root')

        # combine log files
        new_log = open('stats_' + str(mass) + '_' + res_type + '.log','w')
        observed_bound = ""
        expected_lower_bounds = []
        expected_upper_bounds = []
        pwd = os.getcwd()
        logs = glob.glob(os.path.join(pwd,"*.log"))

        job_numbers = []
        for log in logs:
            if not re.search('stats_' + str(mass) + '_' + res_type + '_\d+.log$', log): continue
            job_numbers.append(int(log.replace('.log','').split('_')[-1]))

        job_numbers = sorted(job_numbers)

        for num in job_numbers:
            log_file = open('stats_' + str(mass) + '_' + res_type + '_' + str(num) + '.log','r')
            outputlines = log_file.readlines()
            log_file.close()

            store = True
            for line in outputlines:
                if re.search("observed bound =", line):
                    store = False
                    if observed_bound=="":
                        observed_bound = line
                if re.search("expected bound", line):
                    expected_lower_bounds.append(float(line.split()[4]))
                    expected_upper_bounds.append(float(line.split()[6]))
                if store:
                    new_log.write(line)

        #print observed_bound
        #print expected_lower_bounds
        #print expected_upper_bounds

        median = [0.]
        onesigma = [0., 0.]
        twosigma = [0., 0.]

        getQuantiles(expected_upper_bounds, median, onesigma, twosigma)

        new_log.write('\n')
        new_log.write(observed_bound)
        new_log.write('\n')
        new_log.write('***** expected upper bounds *****\n')
        new_log.write('median: ' + str(median[0]) + '\n')
        new_log.write('+/-1 sigma band: [ ' + str(onesigma[0]) + ' , ' + str(onesigma[1]) + ' ]\n')
        new_log.write('+/-2 sigma band: [ ' + str(twosigma[0]) + ' , ' + str(twosigma[1]) + ' ]\n')
        new_log.close()

        # delete .root and .log file
        os.system('rm stats_' + str(mass) + '_' + res_type + '_*.root')
        os.system('rm stats_' + str(mass) + '_' + res_type + '_*.log')


if __name__ == '__main__':
    main()
