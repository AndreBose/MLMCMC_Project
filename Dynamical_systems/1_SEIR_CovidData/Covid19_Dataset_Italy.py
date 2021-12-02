#########################
#                       #
# Covid19 Italy Dataset #
#                       #
#########################


# Relevant libraries

import os
import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt


# Settings

T0 =   0
Tf = 170

new_dataset_filename = 'covid_dataset_IR.csv'


# Import dataset

covid_data_0 = pd.read_csv('https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv')


# Select relevant columns

covid_data_1 = covid_data_0[['totale_positivi', 'dimessi_guariti', 'deceduti']]
covid_data_1 = covid_data_1.sort_index()


# Create time, I, R columns

covid_data_2 = covid_data_1.copy(deep=True)
covid_data_2['Time'] = np.arange(1,len(covid_data_1)+1,1)

covid_data_2['I'] = covid_data_2['totale_positivi']
covid_data_2['R'] = covid_data_2['dimessi_guariti'] + covid_data_2['deceduti']

covid_data_2 = covid_data_2[['Time', 'I', 'R']]


# Select time interval

covid_data_3 = covid_data_2[T0:Tf]

# Save dataframe in the current folder

covid_data_3.to_csv(os.path.join(os.getcwd(), new_dataset_filename), index=False)


# Plot data

plt.figure(dpi=300)

plt.plot(covid_data_3['Time'], covid_data_3['I'], label = 'I')
plt.plot(covid_data_3['Time'], covid_data_3['R'], label = 'R')

plt.title('Infected and Recovered')
plt.legend()
plt.show


# DataFrame.to_csv(path_or_buf=None, sep=',', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, mode='w', encoding=None, compression='infer', quoting=None, quotechar='"', line_terminator=None, chunksize=None, date_format=None, doublequote=True, escapechar=None, decimal='.', errors='strict', storage_options=None)