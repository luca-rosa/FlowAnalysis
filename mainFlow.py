import glob
import fcswrite
import os
import os.path

from flowFunctions import *

def main():
    path = '/Users/lucarosa/Documents/PhD/Data/Raw/Flow/Concentration_Assay/17_05_17'
    # path = raw_input("Insert FCS path: ")

    os.chdir(path)
    filenames = glob.glob("*.fcs")
    folderName = os.getcwd() + '_gated'

    os.mkdir(folderName)

    fun = False;
    for file in filenames:
        os.chdir(path)
        print('Processing ' + file)

        fcs = FlowCal.io.FCSData(file)

        # Density gate and return contour for plotting
        fcs_gate, contour = gating(fcs)
        imageName = os.path.splitext(file)[0] + '.png'
        imageGatingPath = folderName + '/gating_plots'
        imageMEFPath = folderName + '/MEF_conversion'
        if not os.path.isdir(imageGatingPath):
            os.mkdir(imageGatingPath)
            os.mkdir(imageMEFPath)
        plot(fcs, fcs_gate, contour, imageName, imageGatingPath)

        # Remove zero values and (not) log everything
        #fcs_noZero = noZeros(fcs)

        # Create transformation function in the first iteration
        if fun == False:
            os.chdir(path)
            mef_fun = beadsCalibration(imageMEFPath)
            fun = True

        # Convert to MEF
        fcs_MEF = mef_fun(fcs_gate, channels=['BL1-H', 'YL2-H'])
        os.chdir(folderName)
        processDataDirect(fcs, file)
        #fcsname = os.path.splitext(file)[0] + '_gated.fcs'
        #channels = list(fcs_gate.channels)
        #gated_path = folderName + "/" + fcsname

        #fcswrite.write_fcs(filename = gated_path, chn_names = channels, data = fcs_gate)







main()


#first gate
#tweak the beads reading on the attune
#analyse with R vs Alex gating






















