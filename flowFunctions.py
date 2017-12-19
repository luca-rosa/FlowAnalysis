import FlowCal
import matplotlib.pyplot as plt
import os
import pandas as pd
import glob
import numpy as np
import warnings
import pprint
import time

def noZeros(fcs):
    '''Remove zero values'''

    for channel in ('FSC-H', 'SSC-H', 'SSC-A'):
        mask = fcs[:, channel] > 0
        fcs = fcs[mask, :]
    mask = fcs [:, "SSC-H"] > 200
    fcs = fcs[mask, :]


    return(fcs)


def gating(file):
    '''Return density gated fcs file and the plot'''
    fcs_gate, mask, contour = FlowCal.gate.density2d(file,
                                                     channels = ['FSC-H', 'SSC-H'],
                                                     gate_fraction = 0.70,
                                                     full_output= True)

    return(fcs_gate, contour)

def plot(fcs, fcs_gate, contour, imageName, imagePath):
    os.chdir(imagePath)
    fcs_plot = FlowCal.plot.density_and_hist(fcs,
                                             gated_data = fcs_gate,
                                             gate_contour= contour,
                                             density_channels= ['FSC-H', 'SSC-H'],
                                             density_params= {'mode':'scatter',
                                                              'xlim':[1e2, 1e5],
                                                              'xscale':'log',
                                                              'yscale':'log'
                                                               },
                                             hist_channels = ['BL1-H'],
                                             # hist_channels=['BL1-H', 'YL2-H'],
                                             savefig = imageName
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

            data = data.byteswap().newbyteorder()
            data = pd.Series(data)
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


def processDataDirect(fcs, filename):
    '''Extract the data from every single FCS file and create a txt file'''
    #os.chdir(path)

    # if filename.find("atc") > -1:
    #     inducer = "atc"
    # elif filename.find("iptg") > -1:
    #     inducer = "iptg"
    # elif filename.find("AHL") > -1:
    #     inducer = "AHL"
    # elif filename.find("lpatoL") > -1:
    #     inducer = "lpatoL"

    inducer = "lpatoL"

    if inducer == "AHL":

        if os.path.isfile("data_ahl.txt"):
            dataahl = pd.read_csv("data_ahl.txt")
        else:
            dataahl= pd.DataFrame(columns = ('Concentration', 'Replicate', 'mCherry', 'GFP'))

        num = len(fcs[:, 1])
        splitname = os.path.splitext(filename)

        splitname = splitname[0].split('_')
        Conc = [str(splitname[2])] * num

        Replicate = [int(splitname[3])] * num

        data = pd.DataFrame({'Concentration': Conc,
                             'Replicate': Replicate,
                            'mCherry': fcs[:, 'YL2-H'],
                             'GFP': fcs[:, 'BL1-H']
                             })



        #dataahl = pd.concat([dataahl, data])
        dataahl = dataahl.append(data)
        dataahl.to_csv('data_ahl.txt', ",", header = True, columns = ['Concentration', 'Replicate', 'mCherry', 'GFP'], index = False)

    if inducer == "lpatoL":

        if os.path.isfile("data_lpatoL.txt"):
            datalpatoL = pd.read_csv("data_lpatoL.txt")
        else:
            datalpatoL= pd.DataFrame(columns = ('Concentration', 'Replicate', 'GFP'))

        num = len(fcs[:, 1])
        splitname = os.path.splitext(filename)

        splitname = splitname[0].split('_')
        Conc = [str(splitname[2])] * num

        Replicate = [int(splitname[3])] * num

        data = pd.DataFrame({'Concentration': Conc,
                             'Replicate': Replicate,
                             # 'mCherry': fcs[:, 'YL2-H'],
                             'GFP': fcs[:, 'BL1-H']
                             })



        #datalpatoL = pd.concat([datalpatoL, data])
        datalpatoL = datalpatoL.append(data)
        datalpatoL.to_csv('data_lpatoL.txt', ",", header = True, columns = ['Concentration', 'Replicate', 'GFP'], index = False)




    if inducer == "atc":

        if os.path.isfile("data_atc.txt"):
            dataATC = pd.read_csv("data_atc.txt")
        else:
            dataATC= pd.DataFrame(columns = ('Concentration', 'Replicate', 'mCherry', 'GFP'))

        num = len(fcs[:, 1])
        splitname = os.path.splitext(filename)

        splitname = splitname[0].split('_')
        Conc = [float(splitname[2])] * num

        Replicate = [int(splitname[3])] * num

        data = pd.DataFrame({'Concentration': Conc,
                             'Replicate': Replicate,
                            'mCherry': fcs[:, 'YL2-H'],
                             'GFP': fcs[:, 'BL1-H']
                             })



        #dataATC = pd.concat([dataATC, data])
        dataATC = dataATC.append(data)
        dataATC.to_csv('data_atc.txt', ",", header = True, columns = ['Concentration', 'Replicate', 'mCherry', 'GFP'], index = False)


    if inducer == "iptg":

        if os.path.isfile("data_iptg.txt"):
            dataIPTG = pd.read_csv("data_iptg.txt")
        else:
            dataIPTG = pd.DataFrame(columns = ('Conc', 'Replicate', 'mCherry', 'GFP'))

        num = len(fcs[:, 1])
        splitname = os.path.splitext(filename)

        splitname = splitname[0].split('_')
        Conc = [float(splitname[2])] * num
        Replicate = [int(splitname[3])] * num

        data = pd.DataFrame({'Concentration': Conc,
                             'Replicate': Replicate,
                            'mCherry': fcs[:, 'YL2-H'],
                             'GFP': fcs[:, 'BL1-H']
                             })

        dataIPTG = dataIPTG.append(data)
        dataIPTG.to_csv('data_iptg.txt', ",", header = True, columns = ['Concentration', 'Replicate', 'mCherry', 'GFP'], index = False)


def beadsCalibration(imagePath):
    '''Create au-mef function'''
    warnings.filterwarnings('ignore')

    #Load beads data
    beads_path = raw_input("Insert beads FCS path:")
    os.chdir(beads_path)
    fcs_beads = glob.glob("*.fcs")
    fcs_beads = FlowCal.io.FCSData(fcs_beads[0])
    fcs_beads = FlowCal.transform.to_rfi(fcs_beads)
    os.chdir(imagePath)
    #Gating of the beads
    fcs_beads_g, __, contour = FlowCal.gate.density2d(fcs_beads,
                                                   channels= ['FSC-H', 'SSC-H'],
                                                   gate_fraction= 0.4,
                                                   full_output= True)

    #Plotting of the gated beads
    beadsPlot = FlowCal.plot.density_and_hist(fcs_beads,
                                  gated_data=fcs_beads_g,
                                  gate_contour= contour,
                                  density_channels=['FSC-H', 'SSC-H'],
                                  density_params= {'mode':'scatter',
                                                   'sigma': 5},
                                  hist_channels=['BL1-H', 'YL2-H'],
                                  savefig='gatedBeads.png')

    # mef_values [(gfp/fitc/mefl/BL1),
    #             (mcherry/pe/mepe/YL2)]
    mef_values = np.array([[0, 806, 2159, 5640, 19900, 52630, 172155, 345870],
                           [0, 409, 1250, 3428, 12229, 34294, 113118, 256134]] )

    #Create transformation function
    mef_fun = FlowCal.mef.get_transform_fxn(fcs_beads_g,
                                            mef_values= mef_values,
                                            mef_channels = ['BL1-H', 'YL2-H'],
                                            plot = True,
                                            plot_dir = imagePath
                                            )
    print("MEF transformation function generated")
    return(mef_fun)







