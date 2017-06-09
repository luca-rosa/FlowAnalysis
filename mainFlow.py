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
    folderName = os.getcwd() + '_gatedd'

    # if not os.listdir(folderName):
    #     os.rmdir(folderName)

    os.mkdir(folderName)
    #os.chdir(folderName)


    fun = False;
    for file in filenames:
        os.chdir(path)
        print('Processing ' + file)
        fcs = FlowCal.io.FCSData(file)

        # Remove zero values and (not) log everything
        fcs_noZero = noZeros(fcs)

        imagePath = folderName + '/images'
        if not os.path.isdir(imagePath):
            os.mkdir(imagePath)

        # Create transformation function in the first iteration
        if fun == False:
            mef_fun = beadsCalibration(imagePath)
            fun = True


        # Convert to MEF
        fcs_MEF = mef_fun(fcs_noZero, channels=['BL1-H', 'YL2-H'])

        # Density gate and return contour for plotting
        fcs_gate, contour = gating(fcs_MEF)

        fcsname = os.path.splitext(file)[0] + '_gated.fcs'

        channels = list(fcs_gate.channels)
        #print(channels)
        gated_path = folderName + "/" + fcsname

        fcswrite.write_fcs(filename = gated_path, chn_names = channels, data = fcs_gate)

        imageName =  os.path.splitext(file)[0] + '.png'
        plot(fcs_noZero, fcs_gate, contour, imageName, imagePath)
        print('Image saved')

return(folderName)


folderName = main()


#first gate
#tweak the beads reading on the attune
#analyse with R vs Alex gating






















