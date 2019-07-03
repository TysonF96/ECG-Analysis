import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.fftpack import ifft
import math

"""The main function can be modified to extract and plot data, find peaks, apply filters and more."""
def main():
	#Enter code here (see README.TXT for examples)
	return 0

"""Define a data segment class. Each instance corresponds to a given segment of data (which could be ECG or otherwise)
The methods of this class are useful as they allow us to set offsets, find maxima/minima, compute means and more
By default, this class is used for data extracted directly from some file (for this, we use data = [])
To use the methods of this class on other lists (for example, when computing R-R intervals), we can
change the default value for the entry data by passing the array we would like to work with."""
class DataSegment:
	def __init__(self, filename, offset, events, data = []):
		if data == []:
			self.filename = filename
			self.offset = offset
			self.events = events
			self.data = read_data(self.filename, self.offset, self.events)
		else:
			self.data = data
	
	"""This method computes the maximum element of self.data"""
	def mymax(self):
		current_max = self.data[0]
		for element in self.data:
			if element > current_max:
				current_max = element
		return current_max

	"""This method computes the minimum element of self.data"""
	def mymin(self):
		current_min = self.data[0]
		for element in self.data:
			if element < current_min:
				current_min = element
		return current_min

	"""This method computes the mean of self.data"""
	def mean(self):
		mysum = 0
		for element in self.data:
			mysum += element
		mean = float(mysum)/len(self.data)
		return mean

	"""This method computes the standard deviation of self.data"""
	def stddev(self):
		mean = self.mean()
		sum_squares = 0
		for element in self.data:
			sum_squares += (mean - element) ** 2
		stddev = (sum_squares/len(self.data)) ** 0.5
		return stddev

"""Define a subclass of DataSegment that has additional methods used specifically for ECG data"""
class ECG(DataSegment):
	def __init__(self, filename, offset, events):
		self.filename = filename
		self.offset = offset
		self.events = events
		self.data = read_data(self.filename, self.offset, self.events)
		if len(self.data) == 0:
			raise Exception("Offset should not exceed file length and events should be greater than 0.")

	"""This is method for finding peaks in linear time. We set a threshold (defined by the variable tolerance) for the start of a peak depending on the
	mean and standard deviation of the data. The function returns an array where each element is the index of a peak (note that the incides have
	been shifted by the offset, which depends on how the class is initialised).
	UPDATE: We have since added some code to account for potential false T wave peak detections that occur in a02.dat. The code we have added includes
	the edge case at the start of the for loop, the final conditional in the for loop, and the use of found_end_peak."""
	def find_peaks(self, tolerance = 2):
		found_peak = False
		found_end_peak = False
		mean = self.mean()
		stddev = self.stddev()
		peaks = []
		for i in range(0, len(self.data)):
			#Edge case. This will delete a false T wave peak detection at the end of the data segment. It may also delete
			#a true peak if the following trough occurs after the end of the data segment. However, this is not as an issue
			#as this true peak will be picked up when using a data buffer to analyse larger segments of data.
			if (i == len(self.data) - 1) and (found_end_peak):
				del peaks[-1]
			#If we have found a peak, we continue until the data falls back below peak threshold
			if found_peak and self.data[i] <= mean + (tolerance * stddev):
				#Reached end of peak
				found_peak = False
				found_end_peak = True
				end_peak = i
				peaks.append(start_peak + index_of_max(self.data[start_peak:end_peak+1]) + self.offset)
			#If we have not yet found a peak (or we are between peaks), we keep searching until we find the next peak
			elif (not found_peak) and self.data[i] > mean + (tolerance * stddev):
				found_peak = True
				#Here, we delete the previously detected peak if found_end_peak is True (which means the previously detected peak is not followed by a significant
				#trough and so we expect it to be a falsely detected T wave peak)
				if found_end_peak:
					del peaks[-1]
				start_peak = i
			#If the below statement is true, the previous peak we found was an R wave peak (and not a T wave peak)
			#and so we do not need to remove any false peaks from the array peaks
			if found_end_peak and self.data[i] <= mean - stddev:
				found_end_peak = False
		return peaks

	def RR_intervals(self):
		RR_data = []
		peaks = self.find_peaks()
		for i in range(0, len(peaks) - 1):
			RR_data.append(peaks[i + 1] - peaks[i])
		return RR_data

	"""This method finds the peaks of the result of the find_peaks method. This allows us to compute respiratory rate from
	amplitude modulation of the ECG peaks. The function returns a list containing the locations of the found peaks
	as an integer representing its index in the ECG data. The location corresponds to the end of the found peak.
	Recommended range of tolerances is discussed in the final report. Higher tolerance works better for apnoea data and 
	lower tolerance works better for non-apnoea data."""
	def resp_peaks_amplitude(self, tolerance = 2):
		peaks = self.find_peaks()
		amplitude_data = []
		location_data = []
		for peak in peaks:
			amplitude_data.append(self.data[peak - self.offset])
		#Obtain the mean and standard deviation of amplitude data
		amplitude = DataSegment("", 0, 0, amplitude_data)
		amplitude_data_stddev = amplitude.stddev()
		i = 0
		while i < (len(amplitude_data) - 3):
			#If data at index i+1 greater than data at index i, and either data at index i+2 greater than data at index i+1 or index i+3 greater than data
			#at index i+1, we assume we have reached the start of a peak
			if (amplitude_data[i + 1] > amplitude_data[i]) and ((amplitude_data[i + 2] > amplitude_data[i + 1]) or (amplitude_data[i + 3] > amplitude_data[i + 1])):
				#Search for end of peak. We assume that we are still at a peak if the data continues to increase (i.e. data[i+1]>data[i] or data[i+2]>data[i])
				count = 0
				while (i < (len(amplitude_data) - 4)) and ((amplitude_data[i + 3] > amplitude_data[i + 2]) or (amplitude_data[i + 4] > amplitude_data[i + 2])):
					count += 1
					i += 1
				#Check if peak is sufficiently large
				if abs(amplitude_data[i + 2] - amplitude_data[i - count]) > (tolerance * amplitude_data_stddev):
					location_data.append(peaks[i + 2])
				i += 2
			i += 1
		return location_data

	"""This method finds the troughs of the R-R intervals. This allows us to compute respiratory rate from
	frequency modulation of the ECG peaks. The function returns a list containing the locations of the found troughs
	as an integer representing its index in the ECG data. The location corresponds to the end of the found trough. Recommended
	values for tolerance are discussed in the final report."""
	def resp_troughs_frequency(self, tolerance = 1):
		RR_data = self.RR_intervals()
		peaks = self.find_peaks()
		location_data = []
		RR = DataSegment("", 0, 0, RR_data)
		RR_data_stddev = RR.stddev()
		i = 0
		while i < (len(RR_data) - 1):
			#If data at index i+1 smaller than data at index i, we assume we have reached the start of a trough
			if RR_data[i + 1] < RR_data[i]:
				#Search for end of trough. We assume that we are still at a trough if the data continues to decrease (i.e. data[i+1]<data[i] or data[i+2]<data[i])
				count = 0
				while (i < len(RR_data) - 3) and ((RR_data[i + 2] < RR_data[i + 1]) or (RR_data[i + 3] < RR_data[i + 1])):
					count += 1
					i += 1
				#Check if trough is sufficiently large
				if abs(-RR_data[i + 1] + RR_data[i - count]) > (tolerance * RR_data_stddev):
					location_data.append(peaks[i + 1])
				i += 2
			i += 1
		return location_data

"""Define a subclass of DataSegment that has additional methods used specifically for Respiratory data"""
class Respiratory(DataSegment):
	def __init__(self, filename, offset, events):
		self.filename = filename
		self.offset = 4 * offset
		self.events = 4 * events
		self.respdata = read_data(self.filename, self.offset, self.events)
		self.data = []
		for i in range(0, len(self.respdata)):
			if (i + 2) % 4 == 0:
				self.data.append(self.respdata[i])
		if len(self.data) == 0:
			raise Exception("Offset should not exceed file length and events should be greater than 0." + str(self.offset) + " " +str(self.events))

	"""This function is used to find troughs in the oronasal airflow data. These troughs should correspond to inhalation. 
	We use this function to compute the respiratory rate using the oronasal airflow data as a benchmark for our
	algorithm that computes the respiratory rate from the ECG signal. Note that the algorithm used here is essentially
	the same as the algorithm used for the find_peaks function except that we look for data that falls below a certain
	threshold rather than above it.
	The variable tolerance refers to how many standard deviations we would like the data to reach below the mean
	before counting as a trough. Recommended tolerance = 2 for apnoea and tolerance = 1 for no apnoea in a01.dat. We set
	the tolerance to be 1.5 so that it works well for both apnoea and no apnoea data."""
	def find_airflow_troughs(self, tolerance = 1.5):
			found_trough = False
			mean = self.mean()
			stddev = self.stddev()
			troughs = []
			for i in range(0, len(self.data)):
				#If we have found a trough, we continue until the data increases back above trough threshold
				if found_trough and self.data[i] >= mean - (tolerance * stddev):
					#Reached end of trough
					found_trough = False
					end_trough = i
					troughs.append(start_trough + index_of_min(self.data[start_trough:end_trough+1]) + self.offset)
				#If we have not yet found a trough (or we are between troughs), we keep searching until we find the next trough
				elif (not found_trough) and self.data[i] < mean - (tolerance * stddev):
					found_trough = True
					start_trough = i
			return troughs	

"""This function filters data using a Gaussian kernel (low pass). It is useful for removing noise from signals.
The variable scale is a constant factor that scales the attenuated data.
The variable attenuation determines how much the high frequencies are attenuated.
Higher attenuation corresponds to attentuating the high frequencies more."""
def low_pass_filter(data, scale, attenuation):
	#Take the discrete Fourier transform of the data
	transformed_data = fft(data)
	attenuated_data = []
	#Attenuate unwanted frequencies
	for i in range(0, len(data)):
		attenuated_data.append(transformed_data[i] * scale * math.exp(-((attenuation * i) ** 2)))
	#Take the inverse Fourier transform of the data
	filtered_data = ifft(attenuated_data)
	return filtered_data

"""This function can be used to plot ECG data and the associated peaks. The variable ecg represents
an instance of the class ECG where the corresponding data has a sample rate of 100Hz. If the sample rate 
changes, the only variable that would need to change is xlabel, which can be done manually. The boolean peaks
and amplitude_peaks are set to True by default. They correspond to adding a subplot of the R-peaks and the
breaths computed from amplitude modulation respectively. To remove either of these subplots, one can set the
corresponding variable to False."""
def plot_ECG_amplitude_data(ecg, outfilename, peaks = True, amplitude_peaks = True, xlabel = "Time (10ms)", ylabel = "Voltage (A/D)"):
	plt.clf()
	plt.plot(ecg.data, color = "black")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.tight_layout()
	plt.subplot(1, 1, 1)
	if peaks:
		for element in ecg.find_peaks():
			plt.plot(element - ecg.offset, ecg.data[element - ecg.offset], color = "red", ls = "", marker = "o", label = "points")
	if amplitude_peaks:
		for element in ecg.resp_peaks_amplitude():
			plt.plot(element - ecg.offset, ecg.data[element - ecg.offset], color = "blue", ls = "", marker = "o", label = "points")
	plt.savefig(outfilename)
	return 0

"""This function can be used to plot frequency domain data such as R-R intervals. Plotting R-R intervals
will be the primary purpose of this function. The variable ecg represents an instance of the class ECG. The boolean
peaks is set to True by default and corresponds to adding a subplot containing the breaths computed from frequency
modulation. To remove this subplot, set the variable troughs to False. The R-R interval plotted at some time t corresponds
to the time interval between the R-peak at time t and the next R-peak."""
def plot_ECG_frequency_data(ecg, outfilename, troughs = True, xlabel = "Time (10ms)", ylabel = "R-R interval (10ms)"):
	plt.clf()
	RR_data = {}
	peaks = ecg.find_peaks()
	for i in range(0, len(peaks) - 1):
		RR_data[peaks[i]] = peaks[i + 1] - peaks[i]
	RR_intervals = RR_data.values()
	frequency_troughs = ecg.resp_troughs_frequency()
	plt.plot([element - ecg.offset for element in peaks[0:-1]], RR_intervals, color = "black")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.tight_layout()
	plt.subplot(1, 1, 1)
	if troughs:
		for element in frequency_troughs:
			plt.plot(element - ecg.offset, RR_data[element], color = "blue", ls = "", marker = "o", label = "points")
	plt.savefig(outfilename)
	return 0

"""This function can be used to plot oronasal airflow data and the associated troughs. To plot data without troughs,
set the variable troughs to False. We have the option to overlay ECG data onto this plot by setting overlay to
True. To overlay the ECG data, call the corresponding ECG plot function before calling this function and leave 
overlay as its default value False."""
def plot_airflow_data(ecg, outfilename, troughs = True, overlay = False, xlabel = "Time (10ms)", ylabel = "Pressure"):
	if not overlay:
		plt.clf()
	plt.plot(ecg.data, color = "black")
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.tight_layout()
	plt.subplot(1, 1, 1)
	if troughs:
		for element in ecg.find_airflow_troughs():
			plt.plot(element - ecg.offset, ecg.data[element - ecg.offset], color = "red", ls = "", marker = "o", label = "points")
	plt.savefig(outfilename)
	return 0

"""This function returns the index of the maximum value in an array of data"""
def index_of_max(array):
	current_max = array[0]
	max_index = 0
	for i in range(0, len(array)):
		if array[i] > current_max:
			current_max = array[i]
			max_index = i
	return max_index

"""This function returns the index of the minimum value in an array of data"""
def index_of_min(array):
	current_min = array[0]
	min_index = 0
	for i in range(0, len(array)):
		if array[i] < current_min:
			current_min = array[i]
			min_index = i
	return min_index

"""This function reads data (ECG or otherwise). (Function obtained from Roger with slight modification)"""
def read_data(filename, offset, events):   # extract data
    file = open(filename, 'rb')                # returns data in A/D units
    file.seek(offset * 2)
    data = []
    for i in range(0, events):
        D = file.read(2)
        #If we have reached end of file, then break
        if not D:
        	break
        N = int.from_bytes(D, byteorder = 'little', signed = True)
        data.append(N)
    file.close()
    return data

"""This function is used to extract data from an ECG file including heatbeats and amplitude as a function of time.
Here, we use a data buffer to read the data in segments. The variables filename and respfilename are the file names of the 
ECG data and oronasal airflow data respectively. The variables outfilename1, outfilename2, and outfilename3
correspond to files containing R-peaks data, AM respiratory data + benchmark, and FM respiratory data + benchmark respectively. 
This function also prints data. Descriptions of the data printed and output files are given in readme.txt."""
def extract_ecg_data(filename, respfilename, outfilename1, outfilename2, outfilename3, offset, events, resp_rate_interval = 120):
	"""resp_rate_interval is the time interval (in seconds) for which we compute the respiratory rate
	segment_size is the number of events per segment of data (which is defined by the resp_rate_interval)
	overlap is the number of overlapping events between two adjacent segments. This has been chosen sufficiently large so
	that we don't miss detecting peaks that occur between segments."""
	amp_count = 0
	freq_count = 0
	airflow_count = 0
	segment_size = resp_rate_interval * 100
	overlap = segment_size // 5
	number_of_segments = ((events - segment_size) // (segment_size - overlap)) + 2
	current_event = offset
	current_segment = ECG(filename, current_event, segment_size)
	current_segment_airflow = Respiratory(respfilename, current_event, segment_size)
	file1 = open(outfilename1, "w+")
	file1.write("R-PEAK INDEX, TIME, AMPLITUDE\n")
	file2 = open(outfilename2, "w+")
	file2.write("R-PEAK AMPLITUDE PEAK INDEX, TIME, AMPLITUDE, COMPUTED RESPIRATORY RATE (BREATHS PER MINUTE), BENCHMARK RESPIRATORY RATE (BREATHS PER MINUTE)\n")
	file3 = open(outfilename3, "w+")
	file3.write("R-R INTERVAL TROUGHS INDEX, TIME, AMPLITUDE, COMPUTED RESPIRATORY RATE (BREATHS PER MINUTE), BENCHMARK RESPIRATORY RATE (BREATHS PER MINUTE)\n")
	rms_difference_amplitude = 0
	rms_difference_frequency = 0
	rms_airflow = 0
	for i in range(0, number_of_segments - 1):
		prev_segment = current_segment
		prev_segment_airflow = current_segment_airflow
		current_event += (segment_size - overlap)
		current_segment = ECG(filename, current_event, segment_size)
		current_segment_airflow = Respiratory(respfilename, current_event, segment_size)
		#Compare first peak in current_segment to last peak in previous segment to ensure that we aren't double counting peaks due to the events in the overlap.
		prev_segment_peaks = prev_segment.find_peaks()
		prev_segment_amplitude_peaks = prev_segment.resp_peaks_amplitude()
		prev_segment_frequency_troughs = prev_segment.resp_troughs_frequency()
		prev_segment_airflow_troughs = prev_segment_airflow.find_airflow_troughs()
		current_segment_peaks = current_segment.find_peaks()
		current_segment_amplitude_peaks = current_segment.resp_peaks_amplitude()
		current_segment_frequency_troughs = current_segment.resp_troughs_frequency()
		current_segment_airflow_troughs = current_segment_airflow.find_airflow_troughs()
		#Look at the first peak in current_segment and check if it is in prev_segment. If it is, truncate the prev_segment data so that there are no repeating peaks
		j = 1
		k = 1
		l = 1
		m = 1
		#Truncate prev_segment data
		while ((current_segment_peaks != []) and (len(prev_segment_peaks) >= j)) and (current_segment_peaks[0] <= prev_segment_peaks[-j]):
			j += 1
		prev_segment_peaks = prev_segment_peaks[0:len(prev_segment_peaks) - (j - 1)]
		while ((current_segment_amplitude_peaks != []) and (len(prev_segment_amplitude_peaks) >= k)) and (current_segment_amplitude_peaks[0] <= prev_segment_amplitude_peaks[-k]):
			k += 1
		prev_segment_amplitude_peaks = prev_segment_amplitude_peaks[0:len(prev_segment_amplitude_peaks) - (k - 1)]
		while ((current_segment_frequency_troughs != []) and (len(prev_segment_frequency_troughs) >= l)) and (current_segment_frequency_troughs[0] <= prev_segment_frequency_troughs[-l]):
			l += 1
		prev_segment_frequency_troughs = prev_segment_frequency_troughs[0:len(prev_segment_frequency_troughs) - (l - 1)]
		while ((current_segment_airflow_troughs != []) and (len(prev_segment_airflow_troughs) >= m)) and (current_segment_airflow_troughs[0] <= prev_segment_airflow_troughs[-m]):
			m += 1
		prev_segment_airflow_troughs = prev_segment_airflow_troughs[0:len(prev_segment_airflow_troughs) - (m - 1)]
		#Write prev_segment data to file
		for element in prev_segment_peaks:
			#Time in seconds
			time = float(element) / 100
			#Amplitude in A/D units
			amplitude = prev_segment.data[element - (current_event - (segment_size - overlap))]
			file1.write(str(element) + "," + str(time) + "," + str(amplitude) + "\n")
		for element in prev_segment_amplitude_peaks:
			time = float(element) / 100
			amplitude = prev_segment.data[element - (current_event - (segment_size - overlap))]
			file2.write(str(element) + "," + str(time) + "," + str(amplitude) + ", , " + "\n")
		for element in prev_segment_frequency_troughs:
			time = float(element) / 100
			amplitude = prev_segment.data[element - (current_event - (segment_size - overlap))]
			file3.write(str(element) + "," + str(time) + "," + str(amplitude) + ", , " + "\n")
		#Write respiratory rate data to files
		resp_rate_amplitude = (100 * 60 * float(len(prev_segment_amplitude_peaks))) / (segment_size - overlap)
		resp_rate_frequency = (100 * 60 * float(len(prev_segment_frequency_troughs))) / (segment_size - overlap)
		resp_rate_airflow = (100 * 60 * float(len(prev_segment_airflow_troughs))) / (segment_size - overlap)
		file2.write(" , , , " + str(resp_rate_amplitude) + ", " + str(resp_rate_airflow) + "\n")
		file3.write(" , , , " + str(resp_rate_frequency) + ", " + str(resp_rate_airflow) + "\n")
		#Compute root mean squared difference between algorithm and its benchmark for testing
		rms_difference_amplitude += (resp_rate_amplitude - resp_rate_airflow) ** 2
		rms_difference_frequency += (resp_rate_frequency - resp_rate_airflow) ** 2
		rms_airflow += resp_rate_airflow ** 2
		amp_count += len(prev_segment_amplitude_peaks)
		freq_count += len(prev_segment_frequency_troughs)
		airflow_count += len(prev_segment_airflow_troughs)
	#Write current_segment data to file. This is an edge case.
	for element in current_segment_peaks:
		time = float(element) / 100
		amplitude = current_segment.data[element - current_event]
		file1.write(str(element) + "," + str(time) + "," + str(amplitude) + "\n")
	for element in current_segment_amplitude_peaks:
		time = float(element) / 100
		amplitude = current_segment.data[element - current_event]
		file2.write(str(element) + "," + str(time) + "," + str(amplitude) + "\n")
	for element in current_segment_frequency_troughs:
		time = float(element) / 100
		amplitude = current_segment.data[element - current_event]
		file3.write(str(element) + "," + str(time) + "," + str(amplitude) + "\n")
	resp_rate_amplitude = (100 * 60 * float(len(current_segment_amplitude_peaks))) / (segment_size)
	resp_rate_frequency = (100 * 60 * float(len(current_segment_frequency_troughs))) / (segment_size)
	resp_rate_airflow = (100 * 60 * float(len(current_segment_airflow_troughs))) / (segment_size)
	file2.write(" , , , " + str(resp_rate_amplitude) + ", " + str(resp_rate_airflow) + "\n")
	file3.write(" , , , " + str(resp_rate_frequency) + ", " + str(resp_rate_airflow) + "\n")
	rms_difference_amplitude += (resp_rate_amplitude - resp_rate_airflow) ** 2
	rms_difference_frequency += (resp_rate_frequency - resp_rate_airflow) ** 2
	rms_airflow += resp_rate_airflow ** 2
	rms_difference_amplitude = math.sqrt(rms_difference_amplitude / number_of_segments)
	rms_difference_frequency = math.sqrt(rms_difference_frequency / number_of_segments)
	rms_airflow = math.sqrt(rms_airflow / number_of_segments)
	amp_count += len(current_segment_amplitude_peaks)
	freq_count += len(current_segment_frequency_troughs)
	airflow_count += len(current_segment_airflow_troughs)
	print("RMS difference (AM): " + str(rms_difference_amplitude))
	print("RMS difference (FM): " + str(rms_difference_frequency))
	print("RMS of airflow: " + str(rms_airflow))
	print("Total number of breaths (AM): " + str(amp_count))
	print("Total number of breaths (FM): " + str(freq_count))
	print("Total number of breaths (airflow): " + str(airflow_count))
	file1.close()
	file2.close()
	file3.close()
	return 0

"""This function is used to test the function find_peaks by searching for R-R intervals that differ significantly from the mean
in a given segment of data.
It returns array of tuples (peak number, R-R interval) where peak number identifies a peak corresponding to an R-R interval that differs 
significantly from the mean R-R interval in the data segment"""
def test_find_peaks(peaks, number_of_peaks):
		anomalies = []
		intervals = []
		for i in range(0, len(peaks) - 1):
			intervals.append(peaks[i + 1] - peaks[i])
		sums = 0
		for interval in intervals:
			sums += interval
		interval_mean = float(sums) / len(intervals)
		sum_squares = 0
		for interval in intervals:
			sum_squares += (interval_mean - interval) ** 2
		interval_stddev = (float(sum_squares) / len(intervals)) ** 0.5
		for i in range(0, len(intervals)):
			#check if R-R interval is an anomaly
			if (intervals[i] >= interval_mean + 3 * interval_stddev) or (intervals[i] <= interval_mean - 3 * interval_stddev):
				anomalies.append((number_of_peaks - (len(intervals) - i) + 1, intervals[i]))
		return anomalies

"""This function computes the running variance of some data.
The running variance at each point is the variance of the
n data samples ahead of that point. Here, n = segment_size."""
def compute_running_variance(data, segment_size = 1000):
	variance_data = []
	for i in range(0, len(data) - segment_size):
		sums = 0
		sum_squares = 0
		for j in range(0, segment_size):
			sums += data[i + j]
		mean = float(sums) / segment_size
		for j in range(0, segment_size):
			sum_squares += (mean - data[i + j]) ** 2
		variance = float(sum_squares) / segment_size
		variance_data.append(variance)
	return variance_data

"""Call main() function on execution."""
if __name__ == "__main__":
	main()
