import numpy as np
from scipy.signal import find_peaks, peak_prominences
from helper import split_bin, compare_neighbour, calc_metrics, peak_location


def peak_seeker(data, dist, pro_sens, pro_min):
    """
    input data from original datafile, calculate how many peaks may exist in the binned data
    dist: minimal distance requirement between peaks when there are multiple peaks existing
    pro_sens: the lowest height of the lower peak compared to the higher peak (used ratio), controls the
    sensitivity of peak recognition
    pro_min: minimum height of peak as a fixed number
    return the locations of peaks in an array
    """
    # check if there is only one UMI read
    if len(data) == 1:
        return np.array([])

    # split the data into bins, with the binwidth to be 0.03
    counts = split_bin(data)

    distance = int(len(counts) / 2)
    # first round of peak check
    peaks, _ = find_peaks(
        np.concatenate(([min(counts)], counts, [min(counts)])),
        distance=distance,
        prominence=max(int(max(counts) / pro_sens), pro_min),
    )

    bin_range_cutoff = 75
    if len(peaks) == 1 and len(counts) >= bin_range_cutoff:
        c1 = counts[: int(len(counts) / 2)]
        c2 = counts[int(len(counts) / 2) :]

        peaks1, _ = find_peaks(
            np.concatenate(([min(c1)], c1, [min(c1)])),
            distance=max(len(c1) - 2, 1),
            prominence=pro_min,
        )

        peaks2, _ = find_peaks(
            np.concatenate(([min(c2)], c2, [min(c2)])),
            distance=max(len(c2) - 2, 1),
            prominence=pro_min,
        )

        if len(peaks1) + len(peaks2) > 1:
            peaks = np.concatenate((peaks1, peaks2 + len(c1)))

    # helps prevent peaks being too close
    if len(peaks) == 2:
        if peaks[1] - peaks[0] < dist:
            peak_heights = peak_prominences(
                np.concatenate(([min(counts)], counts, [min(counts)])), peaks
            )[0]
            if peak_heights[0] > peak_heights[1]:
                peaks = np.array([peaks[0]])
            else:
                peaks = np.array([peaks[1]])

    # helps prevent peaks being too close
    if len(peaks) == 3:
        if peaks[2] - peaks[1] < dist:
            peak_heights = peak_prominences(
                np.concatenate(([min(counts)], counts, [min(counts)])), peaks
            )[0]
            if peak_heights[1] > peak_heights[2]:
                peaks = np.array([peaks[0], peaks[1]])
            else:
                peaks = np.array([peaks[0], peaks[2]])

    # make sure the peak is higher than its 6 (3 on left 3 on right) closest neighbors
    peaks = compare_neighbour(counts, peaks)

    return peaks
