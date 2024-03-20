# This functions takes three required parameters py_dir, data, data_dir
# py_dir, the directory of r_transition.py
# data vector that contains log10(UMI) values
# data_dir the place where data vector is stored temperarily
# and three optional parameters for the peak_seeker function, if not given, using the given fine-tuned values
# return the count of peaks in the log10 UMI file
call_py_script <- function(py_dir, data, data_dir, dist = 36, pro_sens = 12, pro_min = 5) {
    write.csv(data, data_dir)
    content = paste('python', py_dir, data_dir, dist, pro_sens, pro_min)
    peak_count = system(content, intern = TRUE)
    system(paste('rm', data_dir))
    return (peak_count)
}

