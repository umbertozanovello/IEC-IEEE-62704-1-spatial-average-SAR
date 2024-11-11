import numpy as np

reportDirectory = "./CodeTesting/"
reportFilename = "iecReport_unsorted.txt"

data = np.loadtxt(reportDirectory+reportFilename)

arg_sort = np.argsort(data[:,0])

data = data[arg_sort]

np.savetxt(reportDirectory+"iecReport_sorted_10g_temp.txt", data[:,1:], fmt = ["%d", "%d", "%d", "%d", "%.6e", "%.6e", "%d", "%.6e", "%.6e"])