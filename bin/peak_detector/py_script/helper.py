import math
import numpy as np


def split_bin(num):
    """
    This function helps split the log10(Number of UMI) data into bins with desinated size.
    output counts, a list that contains the number of barcodes falling in each bin.
    """
    min_num = min(num)
    max_num = max(num)
    binwidth = 0.03
    # initialize each bin
    counts = [0] * math.ceil((max_num - min_num) / binwidth)
    # loop through each UMI value, put them into bins
    for i in num:
        index = math.floor((i - min_num) / binwidth)
        counts[index] = counts[index] + 1
    return counts


def compare_neighbour(histdata, peakdata):
    """
    compare the height of each peak with its neighbors. This can prevent
    a peak only being a local spike.
    histdata: a list of numbers of barcodes falling in each bin.
    peakdata: positions of current peaks.
    return a list that contains positions of qualified peaks.
    """
    # only qualified peaks will be put into this list
    peaks = []
    # loop through each peak
    for p in peakdata:
        # obtain the height of the peak, using p-1 because of modifications made in find_peaks function
        top = histdata[p - 1]
        ifpeak = True
        zero_count = 0
        for i in range(p - 4, p + 3):
            # all the nearby neighbors cannot be higher than the peak
            if i >= 0 and i < len(histdata) and histdata[i] > top:
                ifpeak = False
            # there should not be many 0s nearby
            if i >= 0 and i < len(histdata) and histdata[i] == 0:
                zero_count += 1
        # too many 0
        if zero_count >= 3:
            ifpeak = False
        # this peak is at the edge of the plot, with at least one neighbor has 0 height
        if p == len(histdata) and zero_count >= 1:
            ifpeak = False
        # for qualified peak, add it to the list
        if ifpeak:
            peaks.append(p)
    return np.array(peaks)


def calc_metrics(verdict, predict):
    """
    calculate the accuracy
    verdict: ground truth, in string format, values can be '0', '1', '2', '1or2'
    predict: model prediction, in integer format, values can be 0, 1, 2
    output the accuracy value
    return accuracy, true positive rate, true negative rate,
    false positive rate, false negative rate
    """
    point = 0
    tp_count = 0
    tn_count = 0
    fp_count = 0
    fn_count = 0

    print("-------------------------------")

    for i in range(len(verdict)):
        # if the ground truth is '1or2'
        if len(verdict[i]) > 1:
            if int(predict[i] == 1) or int(predict[i] == 2):
                point += 1
                tp_count += 1
                tn_count += 1
        else:
            if verdict[i] == "1" and predict[i] == 1:
                point += 1
                tn_count += 1
            elif verdict[i] == "1" and predict[i] != 1:
                fp_count += 1
                print(
                    "mismatch: sample #{}, ground truth: {} peak(s), prediction: {} peak(s)".format(
                        i, verdict[i], predict[i]
                    )
                )
            elif verdict[i] == "2" and predict[i] == 2:
                point += 1
                tp_count += 1
            elif verdict[i] == "2" and predict[i] != 2:
                fn_count += 1
                print(
                    "mismatch: sample #{}, ground truth: {} peak(s), prediction: {} peak(s)".format(
                        i, verdict[i], predict[i]
                    )
                )
            # in case where there is 0 peak or more than 2 peaks
            else:
                if int(verdict[i]) == predict[i]:
                    point += 1
                else:
                    print(
                        "mismatch: sample #{}, ground truth: {} peak(s), prediction: {} peak(s)".format(
                            i, verdict[i], predict[i]
                        )
                    )

    # calculate the metrics
    accuracy = point / len(predict)
    TP = tp_count / (tp_count + fn_count)
    TN = tn_count / (tn_count + fp_count)
    FP = fp_count / (fp_count + tn_count)
    FN = fn_count / (fn_count + tp_count)

    print("-------------------------------")
    print(
        "accuracy: {}\ntrue positive rate: {}\ntrue negativerate: {}\nfalse positive rate: {}\nfalse negative rate: {}".format(
            accuracy, TP, TN, FP, FN
        )
    )
    print("-------------------------------")

    return accuracy, TP, TN, FP, FN


def peak_location(peaks, start, binwidth):
    """
    convert the location of peaks back to original scale
    """
    locations = []
    for p in peaks:
        locations.append(start + (p - 1) * binwidth)
    return locations
