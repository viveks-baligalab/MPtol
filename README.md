# MPtol
@Author: Vivek Srinivas - Baliga Lab, ISB

This is a program to calculate drug tolerance from multiwell plate assay

Program requires two files - 'Labels.csv' and 'Readings.txt' for example

(Refer to examples files for format)

    Step 1: Plot growth
            plot_growth_in_wells(results_text_file, labels, instrument = "BioTek"/"BioAnalyzer")
    Step 2: Identify OD threshold for 'Start of growth point determination'
    Step 3: Calculate drug tolerance with command -
            calculate_plot_tolerance(results_text_file, labels, ODT = __,Drug = __, instrument = "BioTek"/"BioAnalyzer" )
