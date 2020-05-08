# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 23:23:49 2018

@author: Vecsei

Spatial calibration of the for W7-X alkali BES diagnostic module for flap
"""
#import matplotlib
#matplotlib.use('Agg')

import os.path
import h5py
import numpy as np
from scipy import misc
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit
from functools import partial
from matplotlib import pyplot as plt
from matplotlib import gridspec
from collections import OrderedDict
import flap
import flap_w7x_abes

class ShotSpatCal(flap.DataObject):
    def __init__(self, shotID, options=None):
        self.shotID = shotID
    
    def full_calib(self, options={}):
        CalcCalibration(self.shotID[:4], options=options)
        
    def read(self, options=None):
        # read the spatcal hdf5 file
        if options is not None:
            options_list = options.keys()
        else:
            options_list = dict()
        if 'Shot spatcal save dir' in options_list:
            spatcal_dir = options['Shot spatcal save dir']
        else:
            spatcal_dir ='./'
        filename = self.shotID + '_spat.cal'
        fn = os.path.join(spatcal_dir,filename)
        with h5py.File(fn, "r", libver='earliest') as h5File: 
            channel_names = np.array(h5File['/Channels'].value)
            fibres = h5File["/Fibres"].value
            self.calibration = h5File["/Calibration"].value[0].decode("utf-8") 
            self.beam_start_im = h5File["/Beam_start_im"].value
            self.beam_end_im = h5File["/Beam_end_im"].value
            fibre_coords_im = h5File["/Fibre_coords_im"].value
            self.beam_start_beam = h5File["/Beam_start_beam"].value
            self.beam_end_beam = h5File["/Beam_end_beam"].value
            fibre_coords_beam = h5File["/Fibre_coords_beam"].value
            self.beam_start_xyz = h5File["/Beam_start_xyz"].value
            self.beam_end_xyz = h5File["/Beam_end_xyz"].value
            fibre_coords_xyz = h5File["/Fibre_coords_xyz"].value
            self.beam_plane_vector = h5File["/Beam_plane_vector"].value
        
        if 'Channels' in options.keys():
            channels = [bytes(channel, "utf-8") for channel in options["Channels"]]
        else:
            channels=[bytes('ABES-'+str(index), 'utf-8') for index in range(1,41)]

        relative_loc = np.zeros(len(channel_names))
        index=0
        for ch in channels:
            relative_loc[np.where(channel_names==ch)[0][0]] = index+1
            index=index+1
        relative_loc = [x-1 if x > 0 else len(channel_names) for x in relative_loc]
        relative_loc = np.argsort(relative_loc)

        channel_names = channel_names[relative_loc]
        self.fibres = fibres[relative_loc]
        self.fibre_coords_im = fibre_coords_im[:, relative_loc]
        fibre_coords_beam = fibre_coords_beam[:,relative_loc]
        fibre_coords_xyz = fibre_coords_xyz[:,relative_loc]
        
        self.data={
                "Channel name": np.array(channel_names),
                "Device x": fibre_coords_xyz[0,:],
                "Device y": fibre_coords_xyz[1,:],
                "Device z": fibre_coords_xyz[2,:],
                "Device R": np.sqrt(fibre_coords_xyz[0,:]**2+fibre_coords_xyz[1,:]**2),
                "Device Z": fibre_coords_xyz[2,:],
                "Beam axis": fibre_coords_beam[0,:]
                }
        return

    def create_coordinate_object(self, dimension_list, coord_name, channel_names=None):
        if not (channel_names is None):
            data = np.zeros(len(channel_names))
            index = 0
            for channel in channel_names:
                data[index] = (self.data[coord_name])[np.where(self.data["Channel name"] == bytes(str(channel), 'utf-8'))[0][0]]
                index = index+1
        else:
            data = self.data[coord_name]

        coord_object = flap.Coordinate(name=coord_name, unit='m', values=data, shape=np.shape(data),
                          mode=flap.CoordinateMode(equidistant=False), dimension_list=dimension_list)

        return coord_object
   
    def generate_shotdata(self, options = {}):
        options_default = {
            'year': self.shotID[:4],
            'Spatial calib. source': os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        'spatcal'),
            'Shot spatcal save dir': os.path.join(os.path.dirname(os.path.abspath(__file__))),
            'Overwrite': False,
            'Datapath': flap.config.get("Module W7X_ABES","Datapath"),
            'Plot': True}
        options = {**options_default, **options}
        
        shotspatcal_filename = os.path.join(options['Shot spatcal save dir'],self.shotID + '_spat.cal')
        print(shotspatcal_filename)
        
        # opening the relevant calibration files 
        filename = os.path.join(options['Spatial calib. source'], options['year'], 'Geometry', 'apdcam_to_cmos.hdf5')
        apdcam_to_cmos = flap.load(filename)

        cmos_apdcam_trans_mat = apdcam_to_cmos['Transformation matrix']
        chan_cent = apdcam_to_cmos['APDCAM channel centers']
        filename = os.path.join(options['Spatial calib. source'], options['year'], 'Geometry', 'cmos_to_real.hdf5')
        projections = flap.load(filename)
        cmos_to_real = projections['CMOS to machine projection']
        real_to_cmos = projections['Machine to CMOS projection']
        opt_axis_norm = projections['Optical axis vector']
        obs_point = projections['Observation point']
        image_x_vector = projections['Direction of image x axis']
        image_y_vector = projections['Direction of image y axis']
        midplane_crosspoint = projections['Midplane crosspoint']
        
        # Get micrometer settings
        datapath_base = options['Datapath']
        datapath = os.path.join(datapath_base, self.shotID)
        xmlfile = os.path.join(datapath, self.shotID + '_config.xml')
        xml = flap.FlapXml()
        xml.read_file(xmlfile)
        try:
            if (xml.head.tag != 'ShotSettings'):
                raise ValueError("No ShotSettings entry found in XML file " + xmlfile + ".")
            if (xml.head.attrib['Experiment'] != "W7-X A-BES"):
                raise ValueError(xmlfile + " is not a W7-X ABES config file.")
        except Exception:
            raise ValueError("File format error in " + xmlfile)
        config = flap_w7x_abes.abes_get_config(xml)
        
        # Reading the beam coordinates
        filename = os.path.join(options['Spatial calib. source'], options['year'], 'Geometry', 'beam.dat')
        with open(filename) as beam_data:
            lines = beam_data.readlines()
            beam_start = np.asarray([float(data) for data in lines[2].rstrip().split(' ') if data != ''])
            beam_end = np.asarray([float(data) for data in lines[3].rstrip().split(' ') if data != ''])
            beam_norm = beam_end-beam_start
            beam_norm = beam_norm/np.linalg.norm(beam_norm)
            beam_plane_vector = np.cross(beam_norm, np.asarray([0,0,1]))
        
        #Calculating the position of the channels on the CMOS image
        shot_channel_cent = OrderedDict()
        for channel in config['signal_list']:
            new_channel_name = channel
            if channel[:4] == 'ABES':
                new_channel_name = channel[5:]
            if new_channel_name in chan_cent.keys():
                lab_config = chan_cent[new_channel_name]
                if type(lab_config[0]) is list:
                    lab_config = lab_config[0]
                h_del = config['H-Micrometer'] - lab_config[0]
                v_del = config['V-Micrometer'] - lab_config[1]
                chan_cent_move = np.dot(cmos_apdcam_trans_mat, np.array([h_del, v_del]))
                curr_chan_cent = np.asarray([lab_config[2]+chan_cent_move[0],
                                             lab_config[3]+chan_cent_move[1]])
                shot_channel_cent.update({channel: curr_chan_cent})
        
        #Calculating the position of the channels in the machine, assuming that they are in the z=0 plane
        for key in shot_channel_cent.keys():
            point_XY = get_points_projection(np.asarray([shot_channel_cent[key]]),
                                                         cmos_to_real)
            point_machine_coord = image_x_vector*point_XY[0][0]+\
                                  image_y_vector*point_XY[0][1]+\
                                  midplane_crosspoint
            # Getting the cross point of observation point - point_machine_coord line on the z=0 plane
            connecting_vector = obs_point- point_machine_coord
            length_along_vector = -obs_point[2]/connecting_vector[2]
            shot_channel_cent[key] = {'Fibre coords im': shot_channel_cent[key],
                                      'Device xyz': obs_point + connecting_vector * length_along_vector}
            # Adding the beam axis coord
            coord_from_beam_start = shot_channel_cent[key]['Device xyz']-beam_start
            shot_channel_cent[key].update({'Beam coord': np.asarray([np.dot(coord_from_beam_start,beam_norm),
                                                                     np.linalg.norm(coord_from_beam_start-beam_norm*np.dot(coord_from_beam_start,beam_norm))])})
        
        # Everything of interest is done. The next part is just saving everything the same manner as it was 
        # for the idl code
        calibconf = MachineCalibConfig()
        beam_points = calibconf.get_image_XY_coord([beam_start, beam_end], obs_point=obs_point, opt_axis_norm=opt_axis_norm,
                                           image_x_vector=image_x_vector,
                                           image_y_vector=image_y_vector)
        beam_im = get_points_projection(np.asarray(beam_points),real_to_cmos)

        object_dict = {'Calibration': str(options['year']),
                       'Channels': np.asarray(shot_channel_cent.keys),
                       'Beam_start_im': beam_im[0,:],
                       'Beam_end_im': beam_im[1,:],
                       'Beam_start_xyz': beam_start,
                       'Beam_end_xyz': beam_end,
                       'Beam_start_beam': np.asarray([0, 0]),
                       'Beam_end_beam': beam_end-beam_start,
                       'Beam_plane_vector': beam_plane_vector}
        channels = []
        fibres = []
        fibre_coords_xyz = []
        fibre_coords_beam = []
        fibre_coords_im = []
        for chan in shot_channel_cent.keys():
            if chan[:4] == 'ABES':
                channels = channels + [chan]
                fibre_coords_xyz = fibre_coords_xyz + [list(shot_channel_cent[chan]['Device xyz'])]
                fibre_coords_im = fibre_coords_im + [list(shot_channel_cent[chan]['Fibre coords im'])]
                fibre_coords_beam = fibre_coords_beam + [shot_channel_cent[chan]['Beam coord']]
        fibre_coords_im = np.transpose(np.asarray(fibre_coords_im))
        fibre_coords_xyz = np.transpose(np.asarray(fibre_coords_xyz))
        fibre_coords_beam = np.transpose(np.asarray(fibre_coords_beam))
        plt.figure()
#        points = np.sort(np.sqrt(fibre_coords_xyz[0]**2+fibre_coords_xyz[1]**2))
#        plt.plot(points, label='new')
        plt.scatter(fibre_coords_im[0], fibre_coords_im[1], color='red')
        plt.plot(np.asarray([beam_im[0,0],beam_im[1,0]]),[beam_im[0,1],beam_im[1,1]], color='red')
        plt.scatter(beam_im[0,0], beam_im[0,1], color='red')
#        plt.scatter(fibre_coords_xyz[0], fibre_coords_xyz[1], color='red')
#        plt.scatter([beam_start[0],beam_end[0]],[beam_start[1],beam_end[1]], color='red')
#        plt.scatter(beam_start[0], beam_start[1], color='red')
        plt.show(block=False)
        plt.pause(0.01)
        plt.legend()
        
#                    channel_names = np.array(h5File['/Channels'].value)
#            fibres = h5File["/Fibres"].value
#            self.calibration = h5File["/Calibration"].value[0].decode("utf-8") 
#            self.beam_start_im = h5File["/Beam_start_im"].value
#            self.beam_end_im = h5File["/Beam_end_im"].value
#            fibre_coords_im = h5File["/Fibre_coords_im"].value
#            self.beam_start_beam = h5File["/Beam_start_beam"].value
#            self.beam_end_beam = h5File["/Beam_end_beam"].value
#            fibre_coords_beam = h5File["/Fibre_coords_beam"].value
#            self.beam_start_xyz = h5File["/Beam_start_xyz"].value
#            self.beam_end_xyz = h5File["/Beam_end_xyz"].value
#            fibre_coords_xyz = h5File["/Fibre_coords_xyz"].value
#            self.beam_plane_vector = h5File["/Beam_plane_vector"].value
        
    
class CalcCalibration:
    def __init__(self, year, options={}):
        options_default = {
            'Spatial calib. source': os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        'spatcal'),
            'type': 'Points',
            'Get APDCAM to CMOS': True,
            'Get CMOS to machine': True,
            'Overwrite': False}
        options = {**options_default, **options}

        self.year = year
        self.fibre_calib_list = []

        if options['Get APDCAM to CMOS'] is True:
            self.process_apdcam_to_cmos_calib(options)
            # self.chan_cent is a directory with keys of the channel names and contains h0, v0, cx0, cy0 
            data_to_save = {'Transformation matrix': self.cmos_apdcam_trans_mat,
                            'APDCAM channel centers': self.chan_cent}
            filename = os.path.join(options['Spatial calib. source'], self.year, 'Geometry', 'apdcam_to_cmos.hdf5')
            print('saving '+filename)
            flap.save(data_to_save, filename)
            print("done")

        if options['Get CMOS to machine'] is True:
            self.process_cmos_to_machine_calib(options=options)
            data_to_save = {'CMOS to machine projection': self.cmos_to_real,
                            'Machine to CMOS projection': self.real_to_cmos,
                            'Optical axis vector': self.calibconf.opt_axis_norm,
                            'Observation point': self.calibconf.obs_point,
                            'Direction of image x axis': self.calibconf.image_x_vector,
                            'Direction of image y axis': self.calibconf.image_y_vector,
                            'Midplane crosspoint': self.calibconf.midplane_crosspoint}

            filename = os.path.join(options['Spatial calib. source'], self.year, 'Geometry', 'cmos_to_real.hdf5')
            print('saving '+filename)
            flap.save(data_to_save, filename)
            print("done")

    def read_fibre_calib_list(self, options):
        options_default={}
        options = {**options_default, **options}

        filename = os.path.join(options['Spatial calib. source'], self.year, 'fibre_calib_images.dat')
        with open(filename) as fibre_list:
            lines = fibre_list.readlines()
            for line_index in range(1,len(lines)):
                line_data = [data for data in lines[line_index].split("\t") if data is not ""]
                self.fibre_calib_list += [[line_data[0], float(line_data[1]), float(line_data[2]), line_data[3][:-1]]]
        return self.fibre_calib_list

    def process_apdcam_to_cmos_calib(self, options):
        options_default={'limit': 100,
                         'calc angle': False,
                         'plot': True}
        options = {**options_default, **options}

        if self.fibre_calib_list == []:
            self.read_fibre_calib_list(options)

        for image_params in self.fibre_calib_list:
            if image_params[-1][-3:] == 'bmp':
                split_filename = image_params[-1].split('.')
                split_filename[-1] = 'png'
                image_params[-1] = '.'.join(split_filename)            
            image = plt.imread(os.path.join(options['Spatial calib. source'], self.year , image_params[-1]))
            image = np.asarray(image)
            # removing the background and the dead pixels
            image = image/np.var(image)
            image = image-np.median(image)
            image[image<options['limit']]=0
            image =  median_filter(image, size=3)

            # adding the location of the center of the channel to the self.fibre_calib_list variable
            x_weight = np.sum(image, axis=0)
            y_weight = np.sum(image, axis=1)

            image_params += [np.average(np.arange(len(x_weight)), weights=x_weight),
                             np.average(np.arange(len(y_weight)), weights=y_weight)]
            if options['plot'] is True:
                plt.contourf(image)
                plt.scatter(np.average(np.arange(len(x_weight)), weights=x_weight), np.average(np.arange(len(y_weight)),
                                       weights=y_weight), color='red')
                plt.show(block=False)
                plt.pause(0.01)

        # Get the transformation matrix with the v and h micrometer
        # Finding the channels with multiple measurements and storing the corresponding data
        channel_meas = np.asarray(self.fibre_calib_list)[:,0]
        multi_chan_meas_id, chan_count = np.unique(channel_meas, return_counts=True)
        multi_chan_meas_id = multi_chan_meas_id[chan_count>1]
        # Collecting the data for the channels with multiple measurements

        chan_cent = {}
        for chan in multi_chan_meas_id:
            chan_cent[chan] = list()
        for image_params in self.fibre_calib_list:
            if image_params[0] in multi_chan_meas_id:
                chan_cent[image_params[0]] += [[image_params[1],image_params[2],image_params[4],image_params[5]]]
            else:
                chan_cent[image_params[0]] = [[image_params[1],image_params[2],image_params[4],image_params[5]]]
        
        # Creating an equation system for obtaining the elements of the transformation matrix
        for chan in multi_chan_meas_id:
            index = 0
            h0, v0, cx0, cy0 = chan_cent[chan][0]
            # The notation is the following A*[dh, dv] = [dcx, dcy] the goal is to find A
            for meas_id in range(1,len(chan_cent[chan])):
                curr_eq = [[chan_cent[chan][meas_id][0]-h0, chan_cent[chan][meas_id][1]-v0, 0, 0],
                           [0, 0, chan_cent[chan][meas_id][0]-h0, chan_cent[chan][meas_id][1]-v0]]
                curr_cent_move = [chan_cent[chan][meas_id][2]-cx0, chan_cent[chan][meas_id][3]-cy0]
                if meas_id == 1 and index == 0:
                    eqsys = np.array(curr_eq)
                    cent_move = np.array(curr_cent_move)
                else:
                    eqsys = np.concatenate((eqsys, curr_eq), axis=0)
                    cent_move = np.concatenate((cent_move,curr_cent_move))
        trans_mat_sol = np.linalg.solve(np.dot(np.transpose(eqsys),eqsys),
                                        np.dot(np.transpose(eqsys),cent_move))
        trans_mat = np.array([trans_mat_sol[:2],trans_mat_sol[2:]])
        self.cmos_apdcam_trans_mat = trans_mat
        
        if options['calc angle'] is True:
            # Ideally the rows of self,trans_mat should be orthogonal. The angle can be calculated as
            cosa = (trans_mat[0,0]*self.trans_mat[1,0]+trans_mat[0,1]*trans_mat[1,1])/\
                   (np.sqrt(trans_mat[0,0]**2+self.trans_mat[0,1] **2)*np.sqrt(trans_mat[1,0]**2+trans_mat[1,1]**2))
            angle = np.arccos(cosa)*180/np.pi # should be around 90
            print("The angle between the horizontal and vertical micrometer: "+str(angle)+" degrees")
        
        # Saving the center of the coordinates - it the APDCAM side is what we need, so the fibre-head configuraiton is
        # needed
        with open(os.path.join(options['Spatial calib. source'], self.year, 'head_fibre_config.dat'), 'r') as head_config:
            line=head_config.readline()
            head_config_map = {}
            while line:
                line = [name for name in head_config.readline().rstrip().split(' ') if name != '']
                if len(line) == 1:
                    line = [name for name in line[0].split('\t') if name != '']
                if len(line) != 0:
                    head_config_map.update({line[1]:line[0]})
        self.chan_cent = {}
        for chan in chan_cent.keys():
            if chan in multi_chan_meas_id:
                curr_chan_data = {head_config_map[chan]: chan_cent[chan][0]}
            else:
                curr_chan_data = {head_config_map[chan]: chan_cent[chan]}
            self.chan_cent.update(curr_chan_data)
        return self.cmos_apdcam_trans_mat

    def process_cmos_to_machine_calib(self, options={}):
        # TURN off Spyder graphics activation!
        # The goal is to find the x,y coordinate of the pixels of the CMOS, assuming that they are on the z=0 plane
        options_default={'Spatial calib. source': os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        'spatcal'),
                         'type': 'Points'}
        options = {**options_default, **options}
        
        # Read the point coordinates in real space
        filename = os.path.join(options['Spatial calib. source'], self.year, 'Geometry', 'calib_image.dat')
        with open(filename, 'r') as pointdata:
            line=pointdata.readline()
            points_rzt = []
            while line:
                line = pointdata.readline().rstrip()
                if len(line) != 0:
                    points_rzt += [[float(data) for data in line.split(' ') if data != '']]
        points_rzt = np.asarray(points_rzt)
        
        # Click on an image to signal the point coordinates in a photo

        self.points_cmos = np.array([])
        plt.ion()
        f = plt.figure()
        gs = gridspec.GridSpec(1,3, width_ratios=[1, 3, 2]) 
        # Adding the 'buttons'
        ax = plt.subplot(gs[0])
        buttons = OrderedDict()
        buttons['Start\nCalibration'] = [True, partial(self.start_calibration, points_rzt, options)]
        buttons['Delete Last\nPoint'] = [True, self.delete_last_selected_point]
        buttons['Add Point'] = [False, None]
        self.clickbuttons = ClickButtonList(ax, buttons)
        # adding the calibration image with the event handling
        filename = os.path.join(options['Spatial calib. source'], self.year, 'Geometry', 'calib_image.png')
        self.calib_image = np.asarray(plt.imread(filename))
        plt.subplot(gs[1])
        self.calib_interface = plt.imshow(self.calib_image)
        self.calib_interface.figure.canvas.mpl_connect('button_press_event', self.add_to_selected_points)
        # adding the reference image
        filename = os.path.join(options['Spatial calib. source'], self.year, 'Geometry', 'reference_image.png')
        ref_image = np.asarray(plt.imread(filename))
        plt.subplot(gs[2])
        plt.imshow(ref_image)
        f.tight_layout()
        plt.pause(0.01)
        plt.show(block=True)

        # Do some kind of fitting the output should be a function with CMOS coordinate -> real space
        # This part is done when clicking on Start calibration which calls the start_calibration() routine

    def add_to_selected_points(self, event):
        axes = self.calib_interface.axes
        if (event.inaxes != axes) or (self.clickbuttons.buttons['Add Point'].active is False): return
        self.points_cmos = np.array(list(self.points_cmos)+[[event.xdata, event.ydata]])
        scatterplot = axes.scatter(self.points_cmos[:,0], self.points_cmos[:,1], color='#FF0000', marker='x')
        axes.draw_artist(scatterplot)
        canvas = self.calib_interface.figure.canvas
        canvas.blit(axes.bbox)
    
    def add_vector_of_points(self, cmos_XY, color = "#00FF00"):
        axes = self.calib_interface.axes
        scatterplot = axes.scatter(cmos_XY[:,0], cmos_XY[:,1], color=color, marker='x')
        axes.draw_artist(scatterplot)
        canvas = self.calib_interface.figure.canvas
        canvas.blit(axes.bbox)
    
    def delete_last_selected_point(self):
        if len(self.points_cmos) == 0: return
        self.points_cmos = self.points_cmos[:-1]

        axes = self.calib_interface.axes
        axes.cla()
        canvas = self.calib_interface.figure.canvas
        self.calib_interface = axes.imshow(self.calib_image)
        scatterplot = axes.scatter(self.points_cmos[:,0], self.points_cmos[:,1], color='#FF0000', marker='x')
        axes.draw_artist(self.calib_interface)
        axes.draw_artist(scatterplot)
        canvas.blit(axes.bbox)

    def start_calibration(self, points_rzt, options={}):
        # Reading the location parameters of the observation point
        if options['type'] == 'Points':
            self.cmos_to_real, self.real_to_cmos = self.calib_points(points_rzt, options=options)
        else:
            raise NotImplementedError('Spatial calibration with fiducial curves not implemented yet.')
    
    def calib_points(self, points_rzt, options = {}):
        options_default = {'Elliptical symmetry': False}
        options = {**options_default, **options}

        if points_rzt.shape[0] > self.points_cmos.shape[0]:
            raise ValueError('Not enough points selected on image')
        if points_rzt.shape[0] < self.points_cmos.shape[0]:
            too_many_points_clicked = True
            while too_many_points_clicked is True:
                self.delete_last_selected_point()
                too_many_points_clicked = points_rzt.shape[0] < self.points_cmos.shape[0]
        self.calibconf = MachineCalibConfig(year=self.year, options=options)
        self.calibconf.get_optical_axis_midplane_crosspoint()
        # Calculating the part of the points perpendicular to the optical axis:
        point_XY = np.zeros([points_rzt.shape[0], 2])
        for index in range(points_rzt.shape[0]):
            point_xyz = np.asarray([points_rzt[index,0]*np.cos(points_rzt[index,2]*np.pi/180),
                                    points_rzt[index,0]*np.sin(points_rzt[index,2]*np.pi/180),
                                    points_rzt[index,1]])
            point_onplane = self.calibconf.get_proj_to_image_plane(point_xyz)

            #Calculating the XY coordinates
            perp_to_optax = point_onplane-self.calibconf.midplane_crosspoint
            point_XY[index,0] = np.dot(perp_to_optax, self.calibconf.image_x_vector)
            point_XY[index,1] = np.dot(perp_to_optax, self.calibconf.image_y_vector)

        real_to_cmos = solve_warp_equation(point_XY, self.points_cmos, options=options)

        if options['Elliptical symmetry'] is False and options['Spherical symmetry'] is False:
            Kx00, Kx10, Kx01, Kx11, Ky00, Ky10, Ky01, Ky11 = real_to_cmos
            cmos_to_real = solve_warp_equation(self.points_cmos, point_XY, options=options)
        else:
            # It is possible that the image is flipped
            cmos_check = get_points_projection(point_XY, real_to_cmos)
#            if np.sum(np.abs(1082-cmos_check[:,1]-self.points_cmos[:,1])) < np.sum(np.abs(cmos_check[:,1]-self.points_cmos[:,1])):
#                real_to_cmos[4] = 1082-real_to_cmos[4]
#                real_to_cmos[5] = -real_to_cmos[5]
#                real_to_cmos[6] = -real_to_cmos[6]
                
#                pass
#            real_to_cmos[0] = 1312-real_to_cmos[0]
#            real_to_cmos[1] = -real_to_cmos[1]
#            real_to_cmos[2] = -real_to_cmos[2]
            Kx00, Kx10, Kx01, Kx11, Ky00, Ky10, Ky01, Ky11 = real_to_cmos
            cmos_to_real = [(-Kx00+Ky00*Kx01/Ky01)/(Kx10-Ky10*Kx01/Ky01),1/(Kx10-Ky10*Kx01/Ky01),-Kx01/Ky01/(Kx10-Ky10*Kx01/Ky01), 0,
                            (-Kx00+Ky00*Kx10/Ky10)/(Kx01-Ky01*Kx10/Ky10),1/(Kx01-Ky01*Kx10/Ky10),-Kx10/Ky10/(Kx01-Ky01*Kx10/Ky10), 0]

        real_check = get_points_projection(self.points_cmos, cmos_to_real)
        print(point_XY)
        print(real_check)
        cmos_check = get_points_projection(point_XY, real_to_cmos)
        print(self.points_cmos)
        print(cmos_check)
        self.add_vector_of_points(cmos_check)

        return cmos_to_real, real_to_cmos

class MachineCalibConfig:
    def __init__(self, year=None, options={}):
        options_default = {}
        options = {**options_default, **options}
    
        if year is not None:
            self.read(year, options)
    
    def read(self, year, options):
        self.year=year
        # Reading the location parameters of the observation point
        filename = os.path.join(options['Spatial calib. source'], self.year, 'Geometry', 'observation.dat')
        with open(filename, 'r') as obsdata:
            # skipping the comments
            for i in range(4):
                obsdata.readline()
            #Getting the observation point
            self.obs_point = np.asarray([float(coordval) for coordval in obsdata.readline().rstrip().split(' ')
                                    if coordval !=''])
            # Getting the vector of the optical axis
            point_on_opt_axis = np.asarray([float(coordval) for coordval in obsdata.readline().rstrip().split(' ')
                                            if coordval !=''])
            self.opt_axis_norm = point_on_opt_axis-self.obs_point
            self.opt_axis_norm = self.opt_axis_norm/np.linalg.norm(self.opt_axis_norm)
            # Getting vector along image x
            self.image_x_vector = np.asarray([float(coordval) for coordval in obsdata.readline().rstrip().split(' ')
                                         if coordval !=''])
            self.image_x_vector = self.image_x_vector - self.opt_axis_norm * np.dot(self.image_x_vector, self.opt_axis_norm)
            self.image_x_vector = self.image_x_vector/np.linalg.norm(self.image_x_vector)
            self.image_y_vector = -np.cross(self.image_x_vector, self.opt_axis_norm)
        
    def get_optical_axis_midplane_crosspoint(self, obs_point=None, opt_axis_norm=None):
        # Calculating where the optical axis crosses the z=0 plane, we will create a plane at this point perpendicular
        # to the optical axis
        if obs_point is not None:
            self.obs_point = obs_point
        if opt_axis_norm is not None:
            self.opt_axis_norm = opt_axis_norm

        length_along_optax = -self.obs_point[2]/self.opt_axis_norm[2]
        self.midplane_crosspoint = self.obs_point + self.opt_axis_norm * length_along_optax
        # The equation of the plane at the self.midplane_crosspoint crossing this point is
        # 0 = (x-self.midplane_crosspoint[0])*self.opt_axis_norm[0] +\
        #     (y-self.midplane_crosspoint[1])*self.opt_axis_norm[1] +\
        #     (z-self.midplane_crosspoint[2])*self.opt_axis_norm[2]

    def get_proj_to_image_plane(self, point_xyz, obs_point=None, opt_axis_norm=None):
        # image plane is a plane in perpendicular to the optical axis, including the point where the optical axis crosses
        # the z=0 plane
        if obs_point is not None:
            self.obs_point = obs_point
        if opt_axis_norm is not None:
            self.opt_axis_norm = opt_axis_norm

        if hasattr(self, 'midplane_crosspoint') is False:
            self.get_optical_axis_midplane_crosspoint()
        # The line connecting point_xyz and self.obs_point crosses the plane defined previously at
        x_onplane = (self.midplane_crosspoint[0]*self.opt_axis_norm[0]+\
                         (point_xyz[0]/(self.obs_point[0]-point_xyz[0])*(self.obs_point[1]-point_xyz[1])+self.midplane_crosspoint[1]-point_xyz[1])*self.opt_axis_norm[1]+\
                         (point_xyz[0]/(self.obs_point[0]-point_xyz[0])*(self.obs_point[2]-point_xyz[2])+self.midplane_crosspoint[2]-point_xyz[2])*self.opt_axis_norm[2])/\
                         (self.opt_axis_norm[0]+self.opt_axis_norm[1]*(self.obs_point[1]-point_xyz[1])/(self.obs_point[0]-point_xyz[0])+self.opt_axis_norm[2]*(self.obs_point[2]-point_xyz[2])/(self.obs_point[0]-point_xyz[0]))
        y_onplane = (self.obs_point[1]-point_xyz[1])/(self.obs_point[0]-point_xyz[0])*(x_onplane-point_xyz[0])+point_xyz[1]
        z_onplane = (self.obs_point[2]-point_xyz[2])/(self.obs_point[0]-point_xyz[0])*(x_onplane-point_xyz[0])+point_xyz[2]
        point_onplane = np.array([x_onplane, y_onplane, z_onplane])
        #Checking, all of the following should be 0
#        obs_to_point = point_onplane-self.obs_point
#        print((x_onplane-self.midplane_crosspoint[0])*self.opt_axis_norm[0]+
#        (y_onplane-self.midplane_crosspoint[1])*self.opt_axis_norm[1]+
#        (z_onplane-self.midplane_crosspoint[2])*self.opt_axis_norm[2])
#        print((x_onplane-point_xyz[0])/(self.obs_point[0]-point_xyz[0])-(y_onplane-point_xyz[1])/(self.obs_point[1]-point_xyz[1]))
#        print((x_onplane-point_xyz[0])/(self.obs_point[0]-point_xyz[0])-(z_onplane-point_xyz[2])/(self.obs_point[2]-point_xyz[2]))
#        print(np.dot(point_onplane-self.midplane_crosspoint, self.opt_axis_norm))
        return point_onplane
    
    def get_image_XY_coord(self, points_xyz, obs_point=None, opt_axis_norm=None, image_x_vector=None, image_y_vector=None):
        if obs_point is not None:
            self.obs_point = obs_point
        if opt_axis_norm is not None:
            self.opt_axis_norm = opt_axis_norm
        if image_x_vector is not None:
            self.image_x_vector = image_x_vector
        if image_y_vector is not None:
            self.image_y_vector = image_y_vector

        if hasattr(self, 'midplane_crosspoint') is False:
            self.get_optical_axis_midplane_crosspoint()
        
        if type(points_xyz) is not list:
            points_xyz = [points_xyz]
        
        points_XY = []
        for point in points_xyz:
            proj = self.get_proj_to_image_plane(point)
            #Calculating the XY coordinates
            perp_to_optax = proj-self.midplane_crosspoint
            point_XY = np.asarray([np.dot(perp_to_optax, self.image_x_vector), np.dot(perp_to_optax, self.image_y_vector)])
            points_XY = points_XY + [point_XY]

        return points_XY

def solve_warp_equation(source_XY, proj_XY, options = {}):
    options_default = {'Elliptical symmetry': False}
    options = {**options_default, **options}

    source_XY = np.asarray(source_XY)
    proj_XY = np.asarray(proj_XY)
    
    if options['Elliptical symmetry'] is True:
        K_mat=np.zeros([np.shape(source_XY)[0], 3])
        K_mat[:,0] = 1
        K_mat[:,1] = source_XY[:,0]
        K_mat[:,2] = source_XY[:,1]
        Kx_vec=proj_XY[:,0]
        Ky_vec=proj_XY[:,1]
    
        Kx00, Kx10, Kx01 = np.linalg.solve(np.dot(np.transpose(K_mat), K_mat), np.dot(np.transpose(K_mat), Kx_vec))
        Ky00, Ky10, Ky01 = np.linalg.solve(np.dot(np.transpose(K_mat), K_mat), np.dot(np.transpose(K_mat), Ky_vec))
        res=np.array([Kx00, Kx10, Kx01, 0, Ky00, Ky10, Ky01, 0])
    elif options['Spherical symmetry'] is True:
        K_mat=np.zeros([np.shape(source_XY)[0]*2, 4])
        K_mat[:,0] = np.concatenate([source_XY[:,0], -source_XY[:,1]])
        K_mat[:,1] = np.concatenate([source_XY[:,1], source_XY[:,0]])
        K_mat[:np.shape(source_XY)[0],2] = 1
        K_mat[np.shape(source_XY)[0]:,3] = 1
        K_vec=np.concatenate([proj_XY[:,0], proj_XY[:,1]])

        Kx10, Kx01, Kx00, Ky00 = np.linalg.solve(np.dot(np.transpose(K_mat), K_mat), np.dot(np.transpose(K_mat), K_vec))
        res=np.array([Kx00, Kx10, Kx01, 0, Ky00, Kx01, -Kx10, 0])
    else:
        K_mat=np.zeros([np.shape(source_XY)[0], 4])
        K_mat[:,0] = 1
        K_mat[:,1] = source_XY[:,0]
        K_mat[:,2] = source_XY[:,1]
        K_mat[:,3] = source_XY[:,0]*source_XY[:,1]
        Kx_vec=proj_XY[:,0]
        Ky_vec=proj_XY[:,1]
    
        Kx00, Kx10, Kx01, Kx11 = np.linalg.solve(np.dot(np.transpose(K_mat), K_mat), np.dot(np.transpose(K_mat), Kx_vec))
        Ky00, Ky10, Ky01, Ky11 = np.linalg.solve(np.dot(np.transpose(K_mat), K_mat), np.dot(np.transpose(K_mat), Ky_vec))
        res=np.array([Kx00, Kx10, Kx01, Kx11, Ky00, Ky10, Ky01, Ky11])
    return res

def get_points_projection(source_XY, proj_vect):
    Kx00, Kx10, Kx01, Kx11, Ky00, Ky10, Ky01, Ky11 = proj_vect
    proj_XY =  np.asarray([Kx00+source_XY[:,0]*Kx10+source_XY[:,1]*Kx01+source_XY[:,0]*source_XY[:,1]*Kx11,
                           Ky00+source_XY[:,0]*Ky10+source_XY[:,1]*Ky01+source_XY[:,0]*source_XY[:,1]*Ky11])
    proj_XY = np.transpose(proj_XY)
    return proj_XY

#-----------------------------------------Graphics----------------------------------------------------------------------
# The following two classes are needed for selecting points if the cmos to machine calibration is called
class ClickButtonList:
    def __init__(self, axes, button_dict):
        num_of_buttons = len(button_dict.keys())
        buttons = axes.barh(np.arange(num_of_buttons), np.ones(num_of_buttons), height=0.5)
        self.buttons = []
        index = 0
        self.buttons = OrderedDict()
        for key in button_dict:
            self.buttons[key] = Button(buttons[index], key,
                                       single_event=button_dict[key][0],
                                       connected_function=button_dict[key][1])
            index += 1
        self.active_button_id = None
        

class Button:
    def __init__(self, button, button_name, single_event=True, connected_function = None):
        self.button = button
        self.button_name = button_name
        self.connected_function = connected_function
        self.active = False
        self.single_event = single_event
        self.connect()
        plt.pause(0.1)
        self.draw_button(state='off')

    def connect(self):
        self.cidpress = self.button.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        if self.single_event is True:
            self.cidpress = self.button.figure.canvas.mpl_connect(
                'button_release_event', self.on_release)

    def on_press(self, event):
        if event.inaxes != self.button.axes:
            return
        contains, attrd = self.button.contains(event)
        if contains is False:
            self.active = False
            self.draw_button(state='off')
        else:
            if self.active == False:
                self.active = True
                self.draw_button(state='on')
            else:
                self.active = False
                self.draw_button(state='off')
                
    def on_release(self, event):
        if self.active is False: return
        self.connected_function()
        self.active = False
        self.draw_button(state='off')
    
    def draw_button(self, state='off'):
        if state is 'off':
            color = '#1f77b4'
        else:
            color = '#800000'
        self.button.set_color(color)
        canvas = self.button.figure.canvas
        axes = self.button.axes
        bbox = axes.get_window_extent().transformed(self.button.figure.dpi_scale_trans.inverted())
        width = int(bbox.width * self.button.figure.dpi/10)
        text = self.button.axes.text(self.button.get_x()+self.button.get_width()/2,
                                    self.button.get_y()+self.button.get_height()/2,
                                    self.button_name, ha='center', va='center', fontsize=width, color='white')
        axes.draw_artist(self.button)
        axes.draw_artist(text)
        axes.axis('off')
        canvas.blit(axes.bbox)