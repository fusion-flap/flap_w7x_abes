#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 18:01:52 2023

@author: refydani
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
from PIL import Image
from tqdm import tqdm
import xmltodict
import xml.etree.ElementTree as ET
import json
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import multivariate_normal
import scipy

class W7X_ABES_beam_alignment():
    def __init__(self, shot=None):
        self.shot=shot
        self.search_dir='/data/W7-X/APDCAM/'
        self.read_cmos_camera_data()
        # self.int_calib_cmos_camera_data() #divides camera data with exposure time
        self.calc_plasma_time()
        self.locate_beam() #locates camera frames when the beam was on and also its position in the picture
        # self.subtract_background()
        # self.plot_average_frame()
        # self.plot_frames()
        #self.plot_cross()        
        
    def read_xml(self):
        dl=glob.glob(self.search_dir+'*'+self.shot+'*')
        if len(dl) > 0:
            for search_dir in dl:
                fl=glob.glob(search_dir+'/'+self.shot+'*.xml')
                if len(fl) == 0: 
                    raise ValueError("Camera data not found in '{:s}'".format(search_dir))
        else:
            raise ValueError("Shot directory not found in '{:s}'".format(search_dir))
        im=[]
        for file in fl:
            tree = ET.parse(file)
            xml_data = tree.getroot()
            xmlstr = ET.tostring(xml_data, encoding='utf-8', method='xml')
            data_dict = dict(xmltodict.parse(xmlstr))
            # xml = xmltodict.parse(file, encoding='utf-8')
        self.xml=data_dict        
    def get_camerainfo(self):
        self.read_xml()
        self.frametime=float(self.xml['ShotSettings']['CMOS']['Frametime']['@Value'])
        self.exptime=float(self.xml['ShotSettings']['CMOS']['Exptime']['@Value'])
        self.framenumber=int(self.xml['ShotSettings']['CMOS']['FrameNumber']['@Value'])
        
    def read_cmos_camera_data(self):
        dl=glob.glob(self.search_dir+'*'+self.shot+'*')
        if len(dl) > 0:
            for search_dir in dl:
                fl=glob.glob(search_dir+'/*.bmp')
                if len(fl) == 0: 
                    raise ValueError("Camera data not found in '{:s}'".format(search_dir))
        else:
            raise ValueError("Shot directory not found in '{:s}'".format(search_dir))
        im=[]
        print('Reading images')
        for file in tqdm(np.sort(fl)):
            img=Image.open(file).convert('L').rotate(-7)
            im.append(np.array(img)[320:820,370:900])
        self.get_camerainfo()
        self.im=np.array(im)
        self.time=np.arange(self.im.shape[0])*self.frametime*1e-3
    
    def int_calib_cmos_camera_data(self):
        self.im=self.im/self.exptime
        
    def plot_average_frame(self):
        plt.figure()    
        plt.imshow(np.average(self.im,axis=0), cmap='gray', vmin=0, vmax=100)  
        plt.show()

    def plot_average_frame_intensity(self):
        plt.figure()    
        plt.plot(self.time,np.average(np.average(self.im,axis=1),axis=1))  
        plt.title('CMOS Camera avergae frame intensity')
        plt.ylabel('int [digit]')
        plt.xlabel('time [s]')
        plt.show()
        
    def calc_plasma_time(self):
        frameint=np.average(np.average(self.im,axis=1),axis=1)
        ind=np.where(frameint>np.average(frameint))
        self.plasma_on_frames=range(np.amin(ind),np.amax(ind))
        self.maxframe=np.where(frameint==np.amax(frameint))[0][0]
        
    def plot_frame(self,frame):
        plt.figure()    
        plt.imshow(self.im[frame,:,:], cmap='gray', vmin=0, vmax=100)  
        plt.show()
        
    def plot_frames(self, vmax=100):
        plt.figure(1)  
        for i in self.plasma_on_frames:
            plt.imshow(self.im[i,:,:], cmap='gray', vmin=0, vmax=vmax)
            plt.annotate('Frame:'+str(i)+' Time: {:3.1f} ms'.format(self.time[i]*1e3),[0,0])
            plt.pause(0.2)
            plt.clf()
        plt.close(1)
        
    def plot_beam_overview(self):
        fig, main_ax = plt.subplots(figsize=(5, 5))
        fig.suptitle("CMOS CAMERA shot#"+self.shot)
        divider = make_axes_locatable(main_ax)
        top_ax = divider.append_axes("top", 1.05, pad=0.1,sharex=main_ax)
        right_ax = divider.append_axes("right", 1.05,pad=0.1,sharey=main_ax)
        
        # make some labels invisible
        top_ax.xaxis.set_tick_params(labelbottom=False)
        right_ax.yaxis.set_tick_params(labelleft=False)
        
        main_ax.set_xlabel('dim 1')
        main_ax.set_ylabel('dim 2')
        top_ax.set_ylabel('Z profile')
        right_ax.set_xlabel('Z profile')
        
        # x, y = np.mgrid[-1:1:.01, -1:1:.01]
        # pos = np.empty(x.shape + (2,))
        # pos[:, :, 0] = x; pos[:, :, 1] = y
        # rv = multivariate_normal([-0.2, 0.2], [[1, 1.5], [0.25, 0.25]])
        # z = rv.pdf(pos)
        x, y = np.mgrid[0:self.im.shape[1]:1, 0:self.im.shape[2]:1]
        # for frame in self.plasma_on_frames:
        for frame,i in zip(self.beam_on_frames,range(len(self.beam_on_frames))):
            z=self.im[frame,:,:]
            z_max = z.max()
            
            # cur_x = 260
            # cur_y = 180
            cur_x=self.maxarr[i,0]
            cur_y=self.maxarr[i,1]
            
            main_ax.imshow(z, origin='lower', cmap='gray')
            main_ax.annotate('Frame:'+str(frame)+' Time: {:3.1f} ms'.format(self.time[frame]*1e3),xy=(0.05, 0.95), xycoords='axes fraction',color='white')
            main_ax.autoscale(enable=False)
            right_ax.autoscale(enable=False)
            top_ax.autoscale(enable=False)
            right_ax.set_xlim(right=z_max)
            top_ax.set_ylim(top=z_max)
            v_line = main_ax.axvline(cur_x, color='r',linewidth=0.75,linestyle='--')
            h_line = main_ax.axhline(cur_y, color='g',linewidth=0.75,linestyle='--')
            v_prof, = right_ax.plot(z[:,int(cur_x)],np.arange(x.shape[0]), 'r-')
            h_prof, = top_ax.plot(np.arange(x.shape[1]),z[int(cur_y),:], 'g-')
            plt.show()
            plt.pause(0.2)
            main_ax.cla()
            top_ax.cla()
            right_ax.cla()
        
        plt.close(fig)
        
    def calc_rot(self):
        dl=glob.glob(self.search_dir+'*'+self.shot+'*')
        if len(dl) > 0:
            for search_dir in dl:
                fl=glob.glob(search_dir+'/*.bmp')
                if len(fl) == 0: 
                    raise ValueError("Camera data not found in '{:s}'".format(search_dir))
        else:
            raise ValueError("Shot directory not found in '{:s}'".format(search_dir))
        plt.figure(2)
        rota=[]
        maxa=[]
        for rot in range(-20,5):
            img=Image.open(np.sort(fl)[self.maxframe]).convert('L').rotate(rot)
            im=(np.array(img)[320:820,370:900])
            rota.append(rot)
            # plt.plot(np.average(im[200:300,:],axis=0))
            # plt.show()
            # plt.pause(0.1)
            maxa.append(np.amax(np.average(im[200:300,:],axis=0)))
            plt.scatter(rot,np.amax(np.average(im[200:300,:],axis=0)))
        rota=np.array(rota)
        maxa=np.array(maxa)
        p=np.polyfit(rota,maxa,2)
        self.rot=-p[1]/(2*p[0])
        print('Camera rotation: {:3.1f} degrees'.format(self.rot))
        # plt.plot(rota,p[0]*rota**2+p[1]*rota+p[2])
            # plt.scatter(rot,np.amax(im[200:300,250]))
    
    # def subtract_background(self):
    #     self.imbg1=self.im[1::2]-self.im[0::2]
    #     self.imbg2=self.im[0::2]-self.im[1::2]
    #     self.timebg=self.time[0::2]
        
    def smooth(self,y, box_pts):
        # box = np.ones(box_pts)/box_pts
        # y_smooth = np.convolve(y, box, mode='same')
        y_smooth=[]
        for i in range(0, len(y)):
            li = int(i - (box_pts-1)/2)
            if li < 0:
                li = 0
            ui = li + box_pts
            if ui > len(y)-1:
                ui = len(y)-1
                li = ui - box_pts
            y_smooth.append(np.mean(y[li:ui]))
        return np.array(y_smooth)
    
    def locate_beam(self):
        # plt.figure(1)
        onframe=[]
        maxarr=[]
        for frame in self.plasma_on_frames:
            im=self.im[frame,:,:]
            bg=self.im[frame,200,:100]
            bgav=np.mean(bg)
            bgdev=np.std(bg)
            xcut=self.im[frame,200,:]
            if xcut.max() > bgav+bgdev*10:
                onframe.append(frame)
                maxx=self.smooth(xcut,10).argmax()
                ycut=self.im[frame,:,maxx]
                maxy=self.smooth(ycut,10).argmax()
                maxarr.append([maxx,maxy])
                # plt.imshow(im)
                # plt.scatter(maxx,maxy)
                # plt.show()
                # plt.pause(0.5)
            self.beam_on_frames=onframe
            self.maxarr=np.array(maxarr)

    def plot_beam_position(self):
        for frame,i in zip(self.beam_on_frames,range(len(self.beam_on_frames))):
            plt.figure(3)
            plt.scatter(self.time[frame],self.maxarr[i,0])
        plt.show()
    
    def calc_bes_signal(self):
        lightprofile_time=[]
        lightprofile=[]
        for frame,i in zip(self.beam_on_frames,range(len(self.beam_on_frames))):
            lightprofile_time.append(self.time[frame])
            prof=np.mean(self.im[frame,:,self.maxarr[i,0]-30:self.maxarr[i,0]+30],axis=1)
            bgprof=np.mean(np.mean(self.im[[frame-1,frame+1],:,self.maxarr[i,0]-30:self.maxarr[i,0]+30],axis=0),axis=1)
            # bgprof=np.mean(self.im[frame-1,:,self.maxarr[i,0]-30:self.maxarr[i,0]+30],axis=1)
            lightprofile.append(prof-bgprof)
        self.lightprofile_time=np.array(lightprofile_time)
        self.lightprofile=np.array(lightprofile)
        
    def plot_bes_signal(self):
        self.calc_bes_signal()
        plt.figure(3)
        for frame,i in zip(self.beam_on_frames,range(len(self.beam_on_frames))):
            plt.plot(self.lightprofile[i,:])
            plt.annotate('Frame:'+str(frame)+' Time: {:3.1f} ms'.format(self.lightprofile_time[i]*1e3),xy=(0.05, 0.95), xycoords='axes fraction',color='black')
            plt.title("ABES lightprofile, shot#"+self.shot)
            plt.show()
            plt.ylim(0,np.amax(self.lightprofile))
            plt.pause(0.2)
            plt.cla()
        plt.close()
if __name__ == '__main__':  
        # shot='20181018.041'
        shot='20230216.066'
        # plt.close('all')
        a=W7X_ABES_beam_alignment(shot=shot)
        # a.plot_average_frame_intensity()
        # a.plot_beam_position()
        a.plot_beam_overview()
        # a.plot_frames(vmax=200)
        # a.calc_rot()
        # a.plot_bes_signal()
        