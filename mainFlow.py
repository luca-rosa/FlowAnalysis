
import os
import os.path
import sys
import shutil

from flowFunctions import *


def main():

    if len(sys.argv) > 1:
        if float(sys.argv[1]) < 20 or float(sys.argv[1]) > 90:
            sys.exit("Error: Insert  between 20 and 90")

        print("Gating at " + str(sys.argv[1]) + "%")
        gating_per = float(sys.argv[1]) / 100

    else:
        gating_per = 0.50
        print("Gating at default 50%")

    path = raw_input("Insert FCS path: ")

    os.chdir(path)
    filenames = glob.glob("*.fcs")
    folderName = os.getcwd() + '_gated'

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
        # Remove zero values
        fcs_noZero = noZeros(fcs)

        # Density gate and return contour for plotting
        fcs_gate, contour = gating(fcs_noZero, gating_per)
        imageName = os.path.splitext(file)[0] + '.png'

        # Save the the gate plot
        plot(fcs, fcs_gate, contour, imageName, imageGatingPath)

        # Convert to MEF
        fcs_MEF = mef_fun(fcs_gate, channels=['BL1-H', 'YL2-H'])
        os.chdir(folderName)

        fcs_MEF = np.log10(fcs_MEF)

        # Outputs csv file
        processData(fcs_MEF, file)


main()
