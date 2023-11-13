"""
Functions to plot the temperature distributions obtained for the purpose of analysis
"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import gc

# gaussian distribution
def gauss(data, min, max):
    mean, std = stats.norm.fit(data)
    T = np.linspace(min, max, 100)
    return stats.norm.pdf(T, mean, std), mean, std

# load the monthly avg temp data
hist_monthly_avg = np.load('hist_monthly_avg_lf.npy')
nat_monthly_avg = np.load('nat_monthly_avg_lf.npy')

# load the daily data
hist_daily_june_avg = np.load('hist_daily_june_avg.npy')
hist_daily_june_max = np.load('hist_daily_june_max.npy')
hist_daily_june_min = np.load('hist_daily_june_min.npy')
nat_daily_june_avg = np.load('nat_daily_june_avg.npy')
nat_daily_june_max = np.load('nat_daily_june_max.npy')
nat_daily_june_min = np.load('nat_daily_june_min.npy')

def plot_dist(thresh, min_temp, max_temp, hist, nat):
    T = np.linspace(min_temp, max_temp, 100)

    # print summary statistics
    print('Hist. mean: '+str(gauss(hist, min_temp, max_temp)[1]))
    print('Hist. std: '+str(gauss(hist, min_temp, max_temp)[2]))
    print('Hist. max: '+str(max(hist)))

    print('Nat. mean: '+str(gauss(nat, min_temp, max_temp)[1]))
    print('Nat. std: '+str(gauss(nat, min_temp, max_temp)[2]))
    print('Nat. max: '+str(max(nat)))

    plt.figure(dpi=(150), tight_layout=True)
    plt.yticks([])

    # plot hist data
    plt.plot(T, gauss(hist, min_temp, max_temp)[0], label = 'hist', color = 'r')
    plt.hist(hist, bins=20, alpha = 0.1, color = 'r', density = True)

    # plot nat data
    plt.plot(T, gauss(nat, min_temp, max_temp)[0], label = 'nat', color = 'b')
    plt.hist(nat, bins=20, alpha = 0.1, color = 'b', density = True)

    plt.axvline(thresh, color = 'k', linestyle = '--')
    plt.xlabel(r'Temperature $(^{\circ}C)$', fontsize = 20)
    plt.xticks(fontsize = 15)
    plt.legend(fontsize = 15)
    #plt.savefig('monthly_june_lf.png')
    plt.close()
    gc.collect()

# plot temp across each day in june
def plot_daily_temps():
    t = np.arange(30)
    plt.figure(dpi = (150), tight_layout=True)
    plt.xlabel('Days since 01/06/2023', fontsize = 20)
    plt.ylabel(r'Temperature $(^{\circ}C)$', fontsize = 20)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize=15)
    plt.xlim(0, 29)
    plt.plot(t, hist_daily_june_avg, color = 'r', label = 'hist')
    plt.fill_between(t, hist_daily_june_max, hist_daily_june_min, color = 'r', alpha = 0.1)
    plt.plot(t, nat_daily_june_avg, color = 'b', label = 'nat')
    plt.fill_between(t, nat_daily_june_max, nat_daily_june_min, color = 'b', alpha = 0.1)
    plt.legend(fontsize = 15)
    plt.savefig('daily_temp.png')

plot_daily_temps()