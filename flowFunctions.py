import FlowCal
import matplotlib.pyplot as plt
import os
import pandas as pd
import glob
import numpy as np

def noZeros(file):
    '''Remove zero values and log10 transform'''


    fcs = FlowCal.io.FCSData(file)
    fcs = FlowCal.gate.high_low(fcs)
    for channel in ('FSC-H', 'SSC-H', 'SSC-A'):
        mask = fcs[:, channel] > 0
        fcs = fcs[mask, :]

    channels = list(fcs.channels)
    #log_fcs = np.log10(fcs[:, channels])
    return(fcs)


def gating(file):
    '''Return density gated fcs file and the plot'''
    fcs_gate, mask, contour = FlowCal.gate.density2d(file,
                                                     channels = ['FSC-H', 'SSC-H'],
                                                     gate_fraction = 0.60,
                                                     full_output= True)

    return(fcs_gate, contour)

def plot(fcs, fcs_gate, contour, imagePath):
    fcs_plot = FlowCal.plot.density_and_hist(fcs,
                                             gated_data = fcs_gate,
                                             gate_contour= contour,
                                             density_channels= ['FSC-H', 'SSC-H'],
                                             density_params= {'mode':'scatter',
                                                              # 'xlim':[1e0, 1e1],
                                                              # 'ylim': [1e0, 1e1],
                                                              # 'xscale':'log',
                                                              # 'yscale':'log'
                                                               },
                                             hist_channels = ['BL1-H', 'YL2-H'],
                                             savefig = imagePath
                                             )
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



def beadsCalibration(fcs, beadsFCS, mef_values = 0):
    fcs_beads = FlowCal.io.FCSData(beadsFCS)

    fcs_beads = FlowCal.transform.to_rfi(fcs)

    fcs_beads_g, __, contour = FlowCal.gate.density2d(b,
                                                   channels= ['FSC', 'SSC'],
                                                   gate_fraction= 0.3,
                                                   full_output= True)

    beadsPlot = FlowCal.plot.density_and_hist(b,
                                  gated_data=b_g,
                                  gate_contour= contour,
                                  density_channels=['FSC', 'SSC'],
                                  density_params= {'mode':'scatter',
                                                   'xlim': [1e2, 1e3],
                                                   'ylim':[1e2, 1e3],
                                                   'sigma': 5},
                                  hist_channels=['FL1'])

    mef_fun = FlowCal.mef.get_transform_fxn(fcs_beads_g,
                                            mef_values= mef_values,
                                            mef_channels = ['BL1-H', 'YL2-H'],
                                            plot = True)


    fcs = FlowCal.io.FCSData(fcs)

    #Transform to au and MEF
    fcs = FlowCal.transform.to_rfi(fcs)
    fcs = FlowCal.gate.density2d(fcs, channels = ['BL1-H', 'YL2-H'])

    FlowCal.plot.hist1d(fcs, channel = ['BL1-H', 'YL2-H'])


    return(fcs)





os.chdir('/Users/lucarosa/Downloads/FlowCal-master/examples')
b, b_g = beadsCalibration('FCFiles/Beads006.fcs')


