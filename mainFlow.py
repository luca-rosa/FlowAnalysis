import glob
import fcswrite

from flowFunctions import *

def main():

    path = raw_input("Insert FCS path: ")

    os.chdir(path)
    filenames = glob.glob("*.fcs")

    folderName = os.getcwd() + '_gated'
    os.mkdir(folderName)
    # os.chdir(folderName)



    for fcs in filenames:

        print('Trimming' + fcs)

        # Remove zero values and log everything
        fcs_noZero = noZeros(fcs)

        # Density gate and return contour for plotting
        fcs_gate, contour = gating(fcs_noZero)

        fcsname = (os.path.splitext(fcs)[0]) + '_gated.fcs'
        channels = list(fcs_gate.channels)

        print(channels)
        print(fcs_gate.shape)

        gated_path = folderName + "/" + fcsname
        print(gated_path)

        fcswrite.write_fcs(filename = gated_path, chn_names = channels, data = fcs_gate)

        imageName =  os.path.splitext(fcs)[0] + '.png'
        imagePath = folderName + '/' + imageName

        plot(fcs_noZero, fcs_gate, contour, imagePath)
        print('Image saved')




main()





















