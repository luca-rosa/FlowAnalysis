
import os
import os.path
import sys
import shutil
import unicodedata
import ast

from configparser import ConfigParser
from flowFunctions import *


def main():

    parser = ConfigParser()
    parser.read_file(open("config.ini"))

    # path = raw_input("Insert FCS path: ")

    # Read config file
    parser = ConfigParser()
    parser.read_file(open("config.ini"))

    path = unicodedata.normalize('NFKD', parser.get("path", "files_path")).encode('ascii', 'ignore')
    beads_path = unicodedata.normalize('NFKD', parser.get("path", "beads_path")).encode('ascii', 'ignore')

    # need to check values here
    ssch_ts = int(parser.get("threshold", "SSC-H"))
    fsch_ts = int(parser.get("threshold", "FSC-H"))

    channels_list = ast.literal_eval(parser.get("channels", "colour_channels"))
    gating_per = float(parser.get("channels", "gating_percentage")) / 100  # need to check this value

    # if len(sys.argv) > 1:
    #     if float(sys.argv[1]) < 20 or float(sys.argv[1]) > 90:
    #         sys.exit("Error: Insert  between 20 and 90")

    #     print("Gating at " + str(sys.argv[1]) + "%")
    #     gating_per = float(sys.argv[1]) / 100

    # else:
    #     gating_per = 0.50
    #     print("Gating at default 50%")

    def cls():
        os.system('cls' if os.name=='nt' else 'clear')

    cls()

    print(' ')
    print(" ")
    print("CONFIG FILE:")
    print(" ")
    print("THRESHOLDS:")
    print(" " + "FSC-H" + str(fsch_ts))
    print(" " + "SSC-H" + str(ssch_ts))
    print(" ")
    print("PATH:")
    print(" " + "Files path:" + " " + path)
    print(" " + "Beads path:" + " " + beads_path)
    print(" ")
    print("CHANNELS")
    print("Colour channles:" + " " + str(channels_list))
    print("Gating percentage:" + " " + str(gating_per * 100))
    print(" ")

    os.chdir(path)
    filenames = glob.glob("*.fcs")
    folderName = os.getcwd() + '_gated'

    if os.path.isdir(folderName):
        if sys.version_info[0] == 2:
            answer = raw_input("The folder with the gated files already exist, do yo want to overwrite it? (Y / N)")
        elif sys.version_info[0] >= 3:  # hopefully python 4 doesnt change it again
            answer = input("The folder with the gated files already exist, do yo want to overwrite it? (Y / N)")
        else:
            sys.exit("don't use python 1")
        if answer.lower() in ['y', 'yes']:
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

    mef_fun = beadsCalibration(imageMEFPath, beads_path)

    for file in filenames:
        os.chdir(path)
        if len(file.split("_")) > 4:
            # couldn't find a nicer way to do this
            new_file_name = file.split("_")[7] + "_" + file.split("_")[8] + "_" + file.split("_")[9] + "_" + file.split("_")[10]
            os.rename(file, new_file_name)
            file = new_file_name

        print('Processing ' + file)

        fcs = FlowCal.io.FCSData(file)
        # Remove zero values
        fcs_noZero = noZeros(fcs, ssch_ts, fsch_ts, channels_list)

        # Density gate and return contour for plotting
        fcs_gate, contour = gating(fcs_noZero, gating_per)
        imageName = os.path.splitext(file)[0] + '.png'

        # Save the the gate plot
        plot(fcs, fcs_gate, contour, imageName, imageGatingPath)

        # Convert to MEF
        fcs_MEF = mef_fun(fcs_gate, channels=channels_list)
        os.chdir(folderName)

        fcs_MEF = np.log10(fcs_MEF)

        # Outputs csv file
        processData(fcs_MEF, file)


main()
