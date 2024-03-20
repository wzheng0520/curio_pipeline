import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from helper import split_bin, compare_neighbour, calc_metrics, peak_location
from peak_detector import peak_seeker
import sys
from tqdm import tqdm
import datetime


list_of_arguments = sys.argv
tablename = sys.argv[1]
file_directory = sys.argv[2]

# plots made based on needs
if len(sys.argv) > 3:
    make_plot = sys.argv[3]
else:
    make_plot = False

if __name__ == "__main__":
    filetable = pd.read_csv(tablename)
    # whether to test labeled data or make predictions on unlabeled data
    prediction = "peaks" not in filetable.columns
    filenames = filetable.filename.tolist()
    sample_name = filetable.SampleName.tolist()
    binwidth = 0.03
    location = []
    pred_count = []
    if not prediction:
        peaks_verdict = filetable.peaks.tolist()
    # loop through each data file based on the table
    for i in tqdm(range(len(filenames))):
        df = pd.read_csv(file_directory + "/" + filenames[i], sep="\t")
        data = df.log10_UMI.tolist()
        # call the function that identify peaks
        peaks = peak_seeker(data, 36, 12, 5)
        peaks = peak_location(peaks, min(data), binwidth)
        location.append(peaks)
        pred_count.append(len(peaks))

        # make plot when needed, save to the folder named plots
        if make_plot:
            plt.figure()
            plt.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth))
            if not prediction:
                plt.title(
                    "sample #{}: {}\npeak label: {}, peak predict: {}".format(
                        i, sample_name[i], peaks_verdict[i], len(peaks)
                    )
                )
            else:
                plt.title(
                    "sample #{}: {}\npeak predict: {}".format(
                        i, sample_name[i], len(peaks)
                    )
                )
            plt.xlabel("log10 (Number of UMI)")
            plt.ylabel("Number of Bead Barcodes")

            [plt.axvline(p + 0.015, c="red", linewidth=0.9) for p in peaks]
            plt.savefig("plots/reads_" + sample_name[i] + ".png", dpi=200)

    # print out related metrics
    if not prediction:
        acc, TP, TN, FP, FN = calc_metrics(peaks_verdict, pred_count)
    # store the output table
    if prediction:
        output = pd.DataFrame(
            {
                "SampleName": filetable.SampleName,
                "PeakLocation": location,
                "PeakModel": pred_count,
            }
        )
    # if testing the accuracy of the model
    else:
        output = pd.DataFrame(
            {
                "SampleName": filetable.SampleName,
                "PeakLocation": location,
                "PeakModel": pred_count,
                "PeakGroundTruth": peaks_verdict,
            }
        )
    output.to_csv(
        "output_{}_{}.csv".format(len(pred_count), datetime.datetime.now()), index=False
    )
