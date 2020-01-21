#####################################################
# @Author: Abhilash Sarwade
# @Date:   2019-11-13 09:55:47
# @email: sarwade@isac.gov.in
# @File Name: binary_read.py
# @Project: PySOLEXS_master

# @Last Modified time: 2020-01-21 13:48:18
#####################################################

import os
import numpy as np

class solexs_header():
    def __init__(self,hdr_data_arr):
        n_data_packets = hdr_data_arr.shape[0]

        #checking for packet header F9 A4 2B B1
        self.hdr_check = (hdr_data_arr[:,0]==249) & (hdr_data_arr[:,1]==164) & (hdr_data_arr[:,2]==43) & (hdr_data_arr[:,3]==177)

        self.frame_id = hdr_data_arr[:,4] # fifth byte

        #sixth byte
        self.det_id = np.bitwise_and(hdr_data_arr[:,5],1)
        self.flare_trigger = np.right_shift(np.bitwise_and(hdr_data_arr[:,5],2),1)
        self.shaping_time = np.right_shift(np.bitwise_and(hdr_data_arr[:,5],12),2)
        self.hv_enable = np.right_shift(np.bitwise_and(hdr_data_arr[:,5],16),4)
        self.reset_enable = np.right_shift(np.bitwise_and(hdr_data_arr[:,5],32),5)

        #day counter
        self.day_count = np.right_shift(np.bitwise_and(hdr_data_arr[:,6],120),3)

        #reference counter
        ref_count1 = np.fliplr(hdr_data_arr[:,6:10]).copy()
        ref_count1.dtype = 'uint32'

        self.ref_count = np.bitwise_and(ref_count1[:,0],134217727) 

        #shaped pulse baseline upload
        shaped_pulse_baseline1 = np.fliplr(hdr_data_arr[:,10:12]).copy()
        shaped_pulse_baseline1.dtype = 'uint16'
        self.shaped_pulse_baseline = shaped_pulse_baseline1[:,0]

        #timing_channel_energy_selection_window_threshold_lower
        timing_channel_thresh_lower1 = np.fliplr(hdr_data_arr[:,12:14]).copy()
        timing_channel_thresh_lower1.dtype = 'uint16'
        self.timing_channel_thresh_lower = timing_channel_thresh_lower1[:,0]

        #cooler current
        self.cooler_current = hdr_data_arr[:,14]        

        #back_contact
        self.back_contact = hdr_data_arr[:,15]

        #sdd temperature
        self.sdd_temp = hdr_data_arr[:,16]

        #timing_channel_energy_selection_window_threshold_higher
        timing_channel_thresh_higher1 = np.fliplr(hdr_data_arr[:,18:20]).copy()
        timing_channel_thresh_higher1.dtype = 'uint16'
        self.timing_channel_thresh_higher = timing_channel_thresh_higher1[:,0]

        #fpga data subtraction
        self.input_data_subtraction = np.left_shift(np.bitwise_and(hdr_data_arr[:,5],192),2) + np.left_shift(np.bitwise_and(hdr_data_arr[:,6],128),1) + np.bitwise_and(hdr_data_arr[:,17],128)

        #gain selection
        self.gain = np.bitwise_and(hdr_data_arr[:,17],127)




        # #sixth byte
        # sixth_byte = np.unpackbits(hdr_data_arr[:,5].reshape(n_data_packets,1),axis=1)

        # self.event_count = np.packbits(np.hstack((np.zeros((n_data_packets,5),dtype='uint8'),sixth_byte[:,:3]))) # 3 MSB bits
        # self.det_id = sixth_byte[:,-1]
        # self.reset_pulse = sixth_byte[:,3]
        # self.hv_enable = sixth_byte[:,4]

        # #seventh byte
        # seventh_byte = np.unpackbits(hdr_data_arr[:,6].reshape(n_data_packets,1),axis=1) #seventh 
        # self.day_count = np.packbits(np.hstack((np.zeros((n_data_packets,3),dtype='uint8'),seventh_byte[:,:5]))) # 5 MSB bits

        
        # eighth_byte = np.unpackbits(hdr_data_arr[:,7].reshape(n_data_packets,1),axis=1) #eighth 
        # ninth_byte = np.unpackbits(hdr_data_arr[:,8].reshape(n_data_packets,1),axis=1) #ninth 
        # tenth_byte = np.unpackbits(hdr_data_arr[:,9].reshape(n_data_packets,1),axis=1) #tenth

        # # self.ref_counter = np.packbits(np.hstack((np.zeros((n_data_packets,5),dtype='uint32'),seventh_byte[:,5:],eighth_byte,ninth_byte,tenth_byte)))
        #self.ref_counter = 



class solexs_spectrum():
    n_channels = 340
    def __init__(self,spectral_data_arr,n_channels=n_channels):
        # n_channels = 340
        n_data_packets = spectral_data_arr.shape[0]
        spectral_data = np.zeros((n_channels,n_data_packets))

        #spectral_data_arr = self.data_full[:,hdr_size:hdr_size+spectral_data_size]

        for i in range(n_channels):
            spectral_data[i,:] = spectral_data_arr[:,2*i]*2**8 + spectral_data_arr[:,2*i+1] #
            # spectral_data[i,:] = np.left_shift(spectral_data_arr[:,2*i],8) + spectral_data_arr[:,2*i+1]
        # self.SDD0 = spectral_data[:,det_id==0]        
        # self.SDD1 = spectral_data[:,det_id==1]
        self.spectra = spectral_data

    def add_spectra(self):
        self.full_spectrum = self.spectra.sum(axis=1)

    def make_lightcurve(self,channel_start=0,channel_stop=n_channels):
        self.lightcurve = self.spectra[channel_start:channel_stop,:].sum(axis=0)




class solexs_lightcurve():
    def __init__(self,temporal_data_arr):
        n_data_packets = temporal_data_arr.shape[0]

        temporal_high = np.zeros((10,n_data_packets))
        temporal_med = np.zeros((10,n_data_packets))
        temporal_low = np.zeros((10,n_data_packets))

        for i in range(10):
            temporal_high[i,:] = temporal_data_arr[:,6*i]*2**8 + temporal_data_arr[:,6*i+1]
            temporal_med[i,:] = temporal_data_arr[:,6*i+2]*2**8 + temporal_data_arr[:,6*i+3]
            temporal_low [i,:]= temporal_data_arr[:,6*i+4]*2**8 + temporal_data_arr[:,6*i+5]

        self.high = temporal_high.T.reshape(n_data_packets*10)
        self.med = temporal_med.T.reshape(n_data_packets*10)
        self.low = temporal_low.T.reshape(n_data_packets*10)


class SDD_data_structure():
    def __init__(self,data_sdd): #sdd_id = 0 or 1
        # det_id = np.bitwise_and(data_full[:,5],1)
        # data_sdd = data_full[det_id==sdd_id,:]
        # data_sdd1 = data_full[det_id==1,:]

        hdr_size = 20 #bytes
        spectral_data_size = 680 #bytes
        timing_data_size = 60 #bytes        

        hdr_data_arr_sdd = data_sdd[:,:hdr_size]
        
        spectral_data_arr_sdd = data_sdd[:,hdr_size:hdr_size+spectral_data_size]

        temporal_data_arr_sdd = data_sdd[:,-timing_data_size:]

        
        self.hdr_data = solexs_header(hdr_data_arr_sdd)
        self.spectral_data = solexs_spectrum(spectral_data_arr_sdd)
        self.temporal_data = solexs_lightcurve(temporal_data_arr_sdd)
    


class read_solexs_binary_data():
    """
    Read SoLEXS binary data and segregate header, spectral and timing data.
    Output: SDD data structure
    """

    def __init__(self,input_filename):
        self.input_filename = input_filename
        self.data_type = 'Raw'

        hdr_size = 20 #bytes
        spectral_data_size = 680 #bytes
        timing_data_size = 60 #bytes

        self.packet_size = hdr_size + spectral_data_size + timing_data_size

        # if not os.path.isfile(filename):

        file_size = os.path.getsize(input_filename)
        self.n_data_packets = int(np.floor(os.path.getsize(input_filename)/self.packet_size))

        if np.mod(file_size,self.packet_size)==0:
            self.left_over_data_flag = 0
        else:
            self.left_over_data_flag = 1

        data_full = self.read_file()

        det_id = np.bitwise_and(data_full[:,5],1)
        data_sdd1 = data_full[det_id==0,:]
        data_sdd2 = data_full[det_id==1,:]


        self.SDD1 = SDD_data_structure(data_sdd1)
        self.SDD2 = SDD_data_structure(data_sdd2)


    def read_file(self):
        # fid = open(self.filename,'rb')
        # bin_data_full = fid.read()
        # bin_data_full_arr = bytearray(bin_data_full)

        data_full = np.fromfile(self.input_filename,dtype='uint8')

        if self.left_over_data_flag:
            print('Left over bytes to be ignored.')
            extra_bytes = np.mod(self.n_data_packets,self.packet_size)
            # bin_data_full_arr = bin_data_full_arr[:-extra_bytes]
            data_full = data_full[:-extra_bytes]

        # data_full = np.array(bin_data_full_arr)

        data_full = data_full.reshape(self.n_data_packets,self.packet_size)
        return data_full

