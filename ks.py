from tkinter import font
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

#read real file and random file to dataframe
def prepareKS(real_file, random_file):
    #read tab delimited file to dataframe
    real_data = pd.read_csv(real_file, sep='\t')
    random_data = pd.read_csv(random_file, sep='\t')
    #get unique codon in real data to list
    real_codon = ["GCA","GCC","GCG","GCT","TGC","TGT","GAC","GAT","GAA","GAG","TTC","TTT","GGA","GGC","GGG","GGT","CAC","CAT","ATA","ATC","ATT","AAA","AAG","CTA","CTC","CTG","CTT","TTA","TTG","ATG","AAC","AAT","CCA","CCC","CCG","CCT","CAA","CAG","AGA","AGG","CGA","CGC","CGG","CGT","AGC","AGT","TCA","TCC","TCG","TCT","TAA","TAG","TGA","ACA","ACC","ACG","ACT","GTA","GTC","GTG","GTT","TGG","TAC","TAT"]
    # for list of unique codon, extract data from real data and random data and merge them
    fig,axs = plt.subplots(8, 8,sharey=True)
    for i in range(8):
        for j in range(8):
            #get codon from list
            codon = real_codon[i*8+j]
            #get data from real data
            real_data_codon = real_data[real_data['codon'] == codon]
            #get data from random data
            random_data_codon = random_data[random_data['codon'] == codon]
            #merge data
            data_codon = pd.merge(real_data_codon, random_data_codon, on='pos')
            data_codon = data_codon.rename(columns={'count_x': 'count_real', 'count_y': 'count_random'})
            data_codon_count = data_codon[['count_real', 'count_random']]
            result = stats.ks_2samp(data_codon_count['count_real'], data_codon_count['count_random'])
            pvalue = result[1]
            #change pvalue digit to scientific notation
            pvalue = "{:.2e}".format(pvalue)
            values, base = np.histogram(data_codon_count['count_real'])
            values_random, base_random = np.histogram(data_codon_count['count_random'])
            cumulative = np.cumsum(values)
            cumulative_random = np.cumsum(values_random)
            #plot data
            axs[i,j].set_title("codon: "+codon + "\np-value: " + str(pvalue), fontsize=6,x=0.7,y=0.1)
            # axs[i,j].text(0.3, 0.8,codon,bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 0.1},fontsize=4)
            axs[i,j].plot(base[:-1], cumulative, c='#04bfc4')
            axs[i,j].plot(base_random[:-1], cumulative_random, c='#f9766d')
            #set title and small font
            # axs[i,j].text(0.1, 0.1, codon , fontsize=8)
            # axs[i,j].text(0.5, 0.1, 'p-value: ' + str(result[1]), horizontalalignment='center', verticalalignment='center', fontsize=8)
    
    # for ax in axs.flat:
    #     ax.set(xlabel='Count', ylabel='ECDF')

    for ax in fig.get_axes():
        ax.label_outer() 
    
    #save plot to pdf
    # plt.savefig('ks.pdf')
    plt.show()

    # # for codon in real_codon:
    # codon = "ATA"
    # #extract data from real data
    # real_data_codon = real_data[real_data['codon'] == codon]
    # #extract data from random data
    # # print(real_data_codon)
    # random_data_codon = random_data[random_data['codon'] == codon]
    # # print(random_data_codon)
    # #merge data
    # data_codon = pd.merge(real_data_codon, random_data_codon, on='pos')
    # #extract count data from merged data
    # #rename count_x to count_real and count_y to count_random
    # data_codon = data_codon.rename(columns={'count_x': 'count_real', 'count_y': 'count_random'})
    # data_codon_count = data_codon[['count_real', 'count_random']]
    # result = stats.ks_2samp(data_codon_count['count_real'], data_codon_count['count_random'])
    # pvalue = result[1]
    # #change pvalue digit to scientific notation
    # pvalue = "{:.2e}".format(pvalue)

    # values, base = np.histogram(data_codon_count['count_real'])
    # values_random, base_random = np.histogram(data_codon_count['count_random'])
    # cumulative = np.cumsum(values)
    # cumulative_random = np.cumsum(values_random)
    # plt.plot(base[:-1], cumulative, c='blue')
    # plt.plot(base_random[:-1], cumulative_random, c='green')
    # #plot X-axis label
    # plt.xlabel("Count")
    # #plot Y-axis label
    # plt.ylabel("ECDF")
    # #plot legend to box on the right
    # # plt.legend(['real', 'random'], loc='upper left')
    # #plot p-value on bottom right corner
    # plt.text(0.7, 0.2, 'codon: ' + codon, horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes)
    # plt.text(0.7, 0.1, 'p-value: ' + str(pvalue), horizontalalignment='left', verticalalignment='bottom', transform=plt.gca().transAxes)
    # #export plot to png file
    # plt.savefig(codon + '.pdf')
    # # plt.show()
    # plt.close()

prepareKS('real.tab', 'random.tab')