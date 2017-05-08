import FlowCal
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import os
import pandas as pd
import glob
import numpy as np

def noZeros(file):
    '''Remove zero values and log10 transform'''
    os.chdir(path)

    fcs = FlowCal.io.FCSData(file)
    fcs = FlowCal.gate.high_low(fcs)
    for channel in ('FSC-H', 'SSC-H', 'SSC-A'):
        mask = fcs[:, channel] > 0
        fcs = fcs[mask, :]

    log_fcs = np.log10(fcs[:, channels])
    return(log_fcs)


def gating(file):
    '''Return density gated fcs file and the plot'''
    fcs_gate, mask, contour = FlowCal.gate.density2d(file,
                                                     channels = ['FSC-H', 'SSC-H'],
                                                     gate_fraction = 0.70,
                                                     full_output= True)

    return(fcs_gate, contour)

def plot(fcs, fcs_gate, contour):
    fcs_plot = FlowCal.plot.density_and_hist(fcs,
                                             gated_data = fcs_gate,
                                             gate_contour= contour,
                                             density_channels= ['FSC-H', 'SSC-H'],
                                             density_params= {'mode':'scatter'},
                                             hist_channels= ['BL1-H']
                                              )
    plt.tight_layout()

    return fcs_plot



def processData(path, atc = False, iptg = False):
    '''Extract the data from every single FCS file and create a txt file'''
    os.chdir(path)

    if atc:
        filenamesATC = glob.glob('*atc*.fcs')
        dataATC = pd.DataFrame(columns = ('mCherry', 'GFP', 'Conc', 'Replicate'))
        for filename in filenamesATC:
            print(filename)
            singleFCS = FlowCal.io.FCSData(filename)
            num = len(singleFCS[:, 1])
            splitname = os.path.splitext(filename)
            splitname = splitname[0].split('_')
            Conc = [float(splitname[2])] * num
            Replicate = [int(splitname[3])] * num
            data = pd.DataFrame({'mCherry': singleFCS[:, 'YL2-H'],
                                 'GFP': singleFCS[:, 'BL1-H'],
                                 'Conc': Conc,
                                 'Replicate': Replicate})
            dataATC = dataATC.append(data)
        dataATC.to_csv('dataATC.txt', "\t", header = True, columns = ['Conc', 'Replicate', 'mCherry', 'GFP'])

    if iptg:
        filenamesIPTG = glob.glob('*iptg*.fcs')
        dataIPTG = pd.DataFrame(columns = ('mCherry', 'GFP', 'Conc', 'Replicate'))
        for filename in filenamesIPTG:
            print(filename)
            singleFCS = FlowCal.io.FCSData(filename)
            num = len(singleFCS[:, 1])
            splitname = os.path.splitext(filename)
            splitname = splitname[0].split('_')
            Conc = [float(splitname[2])] * num
            Replicate = [int(splitname[3])] * num

            data = pd.DataFrame({'mCherry': singleFCS[:, 'YL2-H'],
                                 'GFP': singleFCS[:, 'BL1-H'],
                                 'Conc': Conc,
                                 'Replicate': Replicate})
            dataIPTG = dataIPTG.append(data)

        dataIPTG.to_csv('dataIPTG.txt', "\t", header = True, columns = ['Conc', 'Replicate', 'mCherry', 'GFP'])

    print dataATC.Conc.unique()
    print dataIPTG.Conc.unique()



