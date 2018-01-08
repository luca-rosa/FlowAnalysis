
import os
import os.path
import sys
import shutil

from flowFunctions import *

def main():
    #path = '/Users/lucarosa/Documents/PhD/Data/Raw/Flow/Concentration_Assay/17_06_29/Old'
    path = raw_input  ("Insert FCS path: ")

    os.chdir(path)
    filenames = glob.glob("*.fcs")
    print(sys.argv[0])
    folderName = os.getcwd() + '_gated'
    print(sys.argv[0])
    if os.path.isdir(folderName):
        answer = raw_input("The folder with the gated files already exist, do yo want to overwrite it? (Y/N)")
        if (answer == "Y") or (answer == "yes") or (answer == "Yes") or (answer == 'y') or (answer == "yes"):
            shutil.rmtree(folderName)
            os.mkdir(folderName)
        else:
            sys.exit("Ok, byeeee")
    else:
        os.mkdir(folderName)



    imageGatingPath = folderName + '/gating_plots'
    imageMEFPath = folderName + '/MEF_conversion'
    if not os.path.isdir(imageGatingPath):
        os.mkdir(imageGatingPath)
        os.mkdir(imageMEFPath)

    mef_fun = beadsCalibration(imageMEFPath)


    for file in filenames:
        os.chdir(path)
        print('Processing ' + file)

        fcs = FlowCal.io.FCSData(file)
        # Remove zero values and (not) log everything
        fcs_noZero = noZeros(fcs)

        # Density gate and return contour for plotting
        fcs_gate, contour = gating(fcs_noZero)
        imageName = os.path.splitext(file)[0] + '.png'

        #Save the the gate plot
        plot(fcs, fcs_gate, contour, imageName, imageGatingPath)

        # Convert to MEF
        fcs_MEF = mef_fun(fcs_gate, channels=['BL1-H', 'YL2-H'])
        os.chdir(folderName)

        fcs_MEF = np.log10(fcs_MEF)
        processDataDirect(fcs_MEF, file)

        # fcs_gate = np.log10(fcs_gate)
        # processDataDirect(fcs_gate, file)



main()


























