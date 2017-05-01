"""
This is a proof-of-concept.

I wanna be able to take some data and fit


"""

import os
import subprocess

root_directory_of_data = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/")
data_location = os.path.join(root_directory_of_data, "Partially_Reduced_Spectra")

def fit_line_and_return_raw_output(
    filename = "3b-averaged.hifi", freq=863071, 
    line_params = "0 0 0 3.5 0 7"):

    original_directory = os.getcwd()

    try:
        os.chdir(data_location)


        class_script_text = (
            "file in {0}\n" # 3b-averaged.hifi
            "find\n"
            "get first\n"
            "set mode x auto\n"
            "set unit v v\n"
        #   "plot"
            "mod freq {1}\n" # 863071
        #   "plot"
            "set mode x -40 +40\n"
        #   "plot"
            "set win -5 15\n"
              # WINDOW #1 (low) : -.7
              # WINDOW #1 (up ) : 12.2
            "base\n"
        #   "plot"
            "method gauss\n"
            "lines 1 \"{2}\"\n" # 0 0 0 3.5 0 7
             # Using the cursor, type / to keep last values
             #      any other key to set a line boundary
             #      Setting line 1
            "mini\n"
            "sic\\examine R%HEAD%GAU%RESULT\n"
            "sic\\examine R%HEAD%GAU%ERROR\n"

             # Observation 1 RMS of Residuals :  Base =  7.81E-02  Line =  7.55E-02

             #          Fit results

             # Line     Area                Position            Width           Tpeak
             # 1    0.88238     (  0.161)   4.430 (  0.824)     9.495 (  2.243)  8.72988E-02
            # "visu\n"
            "exit\n".format(filename, freq, line_params))

        path_to_script_with_filename = os.path.join(os.curdir, "first_try_linefit.class")
        target = open(path_to_script_with_filename, 'w')
        target.write(class_script_text)
        target.close()


        p = subprocess.Popen(['class', '-nw', "@first_try_linefit"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()

    except Exception as e:
        os.chdir(original_directory)
        raise e
    finally:
        os.chdir(original_directory)

    return output, error