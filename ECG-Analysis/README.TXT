##############################################################################################################
INTRODUCTION
##############################################################################################################

This code was submitted for assessment as part of the subject COMP90072 at The University of Melbourne. The goal was to efficiently and accurately compute the respiratory rate of an individual using electrocardiography. 

##############################################################################################################
PROGRAM STRUCTURE
##############################################################################################################

This program requires Python3, Pyplot, and SciPy (for FFTs).

This program is used for the analysis of ECG data. The program is structured as follows. We have one superclass DataSegment that is used for simple manipulation of time series data. The data can be obtained in two ways: by reading .dat files or by passing an array. For example, if we have an array of data called my_array, we can create an instance of the class DataSegment by defining my_segment = DataSegment("", 0, 0, my_array). The first three arguments are not important and can be set to anything since we already have an array of data. If we wish to use data obtained from a .dat file, we create an instance of the class my defining my_segment = DataSegment(filename, offset, events). Here, filename is the name of the .dat file, offset is the starting point in the .dat file (in terms of number of events), and events is the number of events to be read. 

Given an instance of DataSegment called my_segment, we can obtain the offset by using self.offset, and we can obtain the data by using self.data. The class DataSegment has four methods: mymax, mymin, mean, and stddev. The methods mymax and mymin return the maximum and minimum values in self.data respectively. The methods mean and stddev return the mean and standard deviation of self.data respectively. 

The class DataSegment has two subclasses: ECG and Respiratory. The subclass ECG is used for studying data obtained from an ECG signal. This subclass has a method find_peaks, which locates the R-peaks of ECG data. The algorithm used in find_peaks is discussed in the final report. The method find_peaks takes an argument called tolerance, which defines the threshold for a peak in terms of the standard deviation of the data. For example, if we choose tolerance = 2 (which is the default value), the algorithm requires the start of a peak to be at least 2 standard deviations above the mean of the ECG data. The method find_peaks returns an array of indices giving the location of the peaks found with the algorithm. Note that the indices returned are not adjusted by the offset of the data. For example, let us assume that we define an instance of the class ECG with offset 1000 and that the algorithm locates a peak at index 100 relative to this offset. Then the index returned will be 1100. Note also that the index returned is the location of the maximum amplitude value of the R-peak.

The method RR_intervals returns an array consisting of the R-R intervals in self.data (where the interval is defined in terms of number of events, not time in seconds). The R-R interval at index i is the difference between the indices of the (i+1)th and ith peaks. 

The method resp_peaks_amplitude is used to estimate the respiratory rate from the ECG signal by studying the amplitude modulation of the ECG R-peaks. The algorithm used is described in the final report. It has an argument tolerance, which is used to set a statistical threshold for the definition of a breath. A tolerance of 2 means that we require the height of a "peak" (not an R-peak but a peak in the R-peak data) to be greater than 2 times the standard deviation of the amplitude data in the given data segment. This method outputs an array of indices corresponding to the location of the breath. This location corresponds to the end of the peak found.

The method resp_troughs_frequency is used to estimate the respiratory rate from the ECG signal by studying the frequency modulation of the ECG R-peaks. The algorithm used is described in the final report. Similar to resp_peaks_amplitude, the argument tolerance is used to set the threshold for a breath. A tolerance of 2 means that we require the height of a "trough" (a trough in the R-R interval data) to be greater than 2 times the standard deviation of the R-R interval data in the given data segment. This method outputs an array of indices corresponding to the location of the breath. This location corresponds to the end of the trough found.

The subclass Respiratory of DataSegment is used to study oronasal airflow data. Similar to the class ECG, it requires a filename, offset, and number of events. The offset and number of events here is the number of oronasal airflow events (not the number of events that includes all 4 respiratory channels). When we initialise self.offset and self.events, we multiply by a factor of 4 to account for the 4 channels. Then, self.data is initialised to the data in only the oronasal airflow channel of the respiratory data. The class Respiratory contains the method find_airflow_troughs, which is used to locate breaths in the oronasal airflow data. The algorithm used for this is described in the final report. The variable tolerance is used to check if we have reached the start of a trough in the airflow data (corresponding to a breath). We assume that we have reached the start of a trough if the airflow data falls below mean - (tolerance * stddev) where mean is the mean of the airflow data in the given data segment and stddev is the standard deviation of the airflow data in a given segment. 

Outside of the class DataSegments, we have some additional functions. The function low_pass_filter(data, scale, attenuation) filters the array data using a Gaussian kernel and returns the filtered data. The scale and attenuation are the contents appearing in the kernel. 

The function extract_ecg_data(filename, respfilename, outfilename1, outfilename2, outfilename3, offset, events, resp_rate_interval = 120) is used to test our AM and FM algorithms on large samples of data and extract relevant data including our computed respiratory rate. The function creates 3 output files with the following data:

Output file 1 format = CSV:
R-PEAK INDEX, CORRESPONDING TIME, CORRESPONDING AMPLITUDE
Output file 2 format = CSV:
R-PEAK AMPLITUDE PEAK INDEX, CORRESPONDING TIME, CORRESPONDING AMPLITUDE, COMPUTED RESPIRATORY RATE (BREATHS PER MINUTE), BENCHMARK RESPIRATORY RATE (BREATHS PER MINUTE) 
Output file 3 format = CSV:
R-R INTERVAL TROUGHS INDEX, CORRESPONDING TIME, CORRESPONDING AMPLITUDE, COMPUTED RESPIRATORY RATE (BREATHS PER MINUTE), BENCHMARK RESPIRATORY RATE (BREATHS PER MINUTE)

In other words, the first CSV file contains index, time and amplitude information of the ECG R-peaks. The next file contains index, time, and amplitude information of the peaks of the ECG R-peaks data (for amplitude modulation method). It also contains the computed respiratory rate (in breaths per minute) for each interval (given in seconds by resp_rate_interval) along with the benchmark respiratory rate. The final file contains index, time, and amplitude information of the troughs of the ECG R-R intervals (for frequency modulation method). It also contains the computed and benchmark respiratory rates. The respiratory rate is computed every resp_rate_interval seconds. This has been set to 120 by default but can be changed.

The function extract_ecg_data also prints data for testing. This includes the root mean squared (RMS) of the computed airflow respiratory rate, the RMS of the difference between the airflow respiratory rate and the AM respiratory rate, and the RMS of the difference between the airflow respiratory rate and the FM respiratory rate. The function also prints the number of breaths computed using our AM and FM algorithms and the number of breaths computed using our benchmark (airflow) algorithm.

There are additional simple functions with brief descriptions that are commented in the code. 

##############################################################################################################
INPUT/OUTPUT EXAMPLES
##############################################################################################################

To use this program, modify the main() function. The below examples show the main function used to produce the output files in the relevant directories. We begin with some plotting examples.

______________________________________________________________________________________________________________
Example1

This example produces three plots. The first two are the amplitude and frequency data for the first 4000 events in a01.dat with no peaks or troughs shown. The third plot is the airflow data for the first 4000 events in a01.dat with troughs (i.e. breaths) shown. The ECG data has not been overlaid onto the respiratory data in this example.
______________________________________________________________________________________________________________

def main():
	myecg = ECG("Data_Files/a01.dat", 0, 4000)
	myresp = Respiratory("Data_Files/a01r.dat", 0, 4000)
	plot_ECG_amplitude_data(myecg, "Output_Examples/Plots/Example1/Amplitude_No_Peaks.png", False, False)
	plot_ECG_frequency_data(myecg, "Output_Examples/Plots/Example1/Frequency_No_Troughs.png", False)
	plot_airflow_data(myresp, "Output_Examples/Plots/Example1/Airflow_Troughs_No_Overlay.png")
	return 0
______________________________________________________________________________________________________________


______________________________________________________________________________________________________________
Example2

This example produces three plots. Two of these plots are the amplitude and frequency data from events 300000 to 300000+4000 in a02.dat. The R-peaks are shown in red. The breaths computed using our amplitude/frequency modulation algorithms are in blue. The third plot is the oronasal airflow data corresponding to events 300000 to 300000+4000 in a02.dat. The read dots are the breaths found using our benchmark algorithm. We have overlaid the plot of the amplitude data to compare the amplitude modulation algorithm to the benchmark algorithm. Notice the order the functions have been called in. The function plot_ECG_amplitude_data comes directly before plot_airflow_data since we are overlaying the amplitude data plot onto the airflow data plot.
______________________________________________________________________________________________________________

def main():
	myecg = ECG("Data_Files/a02.dat", 300000, 4000)
	myresp = Respiratory("Data_Files/a02r.dat", 300000, 4000)
	plot_ECG_amplitude_data(myecg, "Output_Examples/Plots/Example2/Amplitude_Peaks.png")
	plot_airflow_data(myresp, "Output_Examples/Plots/Example2/Airflow_Troughs_Overlay.png", True, True)
	plot_ECG_frequency_data(myecg, "Output_Examples/Plots/Example2/Frequency_Troughs.png")
	return 0
______________________________________________________________________________________________________________


______________________________________________________________________________________________________________
Example3

This example is the same as Example2 except that we have instead overlaid the frequency plot onto the airflow plot (rather than the amplitude plot onto the airflow plot). Notice the order in which the functions are called. 
______________________________________________________________________________________________________________

def main():
	myecg = ECG("Data_Files/a02.dat", 300000, 4000)
	myresp = Respiratory("Data_Files/a02r.dat", 300000, 4000)
	plot_ECG_frequency_data(myecg, "Output_Examples/Plots/Example3/Frequency_Troughs.png")
	plot_airflow_data(myresp, "Output_Examples/Plots/Example3/Airflow_Troughs_Overlay.png", True, True)
	plot_ECG_amplitude_data(myecg, "Output_Examples/Plots/Example3/Amplitude_Peaks.png")
	return 0
______________________________________________________________________________________________________________



We will now show some examples of using the function extract_ecg_data.
______________________________________________________________________________________________________________
Example1

This example produces three CSV files (see extract_ecg_data function description for more details) corresponding to the R-peaks data, AM respiratory data, and FM respiratory data. The data used was from a01.dat and a01r.dat, and we have read events 0 to 2957000 (i.e the entire file).
______________________________________________________________________________________________________________

def main():
	outfile1 = "Output_Examples/CSV_Files/Example1/R-Peaks_Data.csv"
	outfile2 = "Output_Examples/CSV_Files/Example1/AM_Resp_Data.csv"
	outfile3 = "Output_Examples/CSV_Files/Example1/FM_Resp_Data.csv"
	extract_ecg_data("Data_Files/a01.dat", "Data_Files/a01r.dat", outfile1, outfile2, outfile3, 0, 2957000)
	return 0
______________________________________________________________________________________________________________
Console output:

RMS difference (AM): 2.3370602579631212
RMS difference (FM): 2.192491984910345
RMS of airflow: 4.390169665463789
Total number of breaths (AM): 2731
Total number of breaths (FM): 1597
Total number of breaths (airflow): 2059
______________________________________________________________________________________________________________



______________________________________________________________________________________________________________
Example2

The CSV files produced here use the data a02.dat and a02r.dat. We read events 100,000 to 100,000 + 1,000,000.
______________________________________________________________________________________________________________

def main():
	outfile1 = "Output_Examples/CSV_Files/Example2/R-Peaks_Data.csv"
	outfile2 = "Output_Examples/CSV_Files/Example2/AM_Resp_Data.csv"
	outfile3 = "Output_Examples/CSV_Files/Example2/FM_Resp_Data.csv"
	extract_ecg_data("Data_Files/a02.dat", "Data_Files/a02r.dat", outfile1, outfile2, outfile3, 100000, 1100000)
	return 0
______________________________________________________________________________________________________________
Console output:

RMS difference (AM): 6.399834406825118
RMS difference (FM): 6.96018373383172
RMS of airflow: 12.056287823661002
Total number of breaths (AM): 1216
Total number of breaths (FM): 1291
Total number of breaths (airflow): 2039
______________________________________________________________________________________________________________



Finally, we provide a couple of miscellaneous examples. In practice, when studying small segments of data it is easier to not use the function extract_ecg_data. In this case, it is easier to use the peak/trough finding methods directly. For example, to compute the number of R-peaks between events 0 and 10000 in a02.dat, we can define:

def main():
	myecg = ECG("Data_Files/a02.dat", 0, 10000)
	print(len(myecg.find_peaks()))
	return 0

This prints 149 and so we have found 149 R-peaks here. 

If we wish to print the number of breaths in this data (using AM, FM, and benchmark), we simply define:

def main():
	myecg = ECG("Data_Files/a02.dat", 0, 10000)
	myresp = Respiratory("Data_Files/a02r.dat", 0, 10000)
	print(len(myecg.resp_peaks_amplitude()))
	print(len(myecg.resp_troughs_frequency()))
	print(len(myresp.find_airflow_troughs()))
	return 0

This prints the following:

17
15
17

Hence, we have found 17 breaths using our AM algorithm, 15 peaks using our FM algorithm, and 17 peaks using our benchmark (oronasal airflow) algorithm. We can easily modify the tolerance of the corresponding algorithms. For example, to use a tolerance of 2.5 standard deviations for the AM algorithm (the default is 2), we can use the line print(len(myecg.resp_peaks_amplitude(2.5))). A tolerance of 2.5 for AM turns out to work better in the data with apnoea.

The easiest way to get familiar with the full functionality of the program is to experiment by calling various functions in main(). Have fun!
