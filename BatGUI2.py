#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 10:19:39 2020

@author: devan
"""


import numpy as np
from scipy.signal import spectrogram, butter, lfilter
import time
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg)
from tkinter import filedialog

def amplitude(filename):
    data = []
    ldata = []
    rdata = []
    file = open(filename, newline='\n')
    for _ in file:
        i = int(_.strip())
        j = (((i - 0) * (1.65 - (-1.65))) / (4095-0))
        data.append(j)
    lavg = sum(data[0:N])/len(data)//2
    ravg = sum(data[N:])/len(data)//2
    for _ in range(len(data)//2):
        ldata.append(data[_] - lavg)
    for _ in range(len(data)//2, len(data)):
        rdata.append(data[_] - ravg)
    return ldata, rdata

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def build_og(Plot, runtime, amp, side='left'):
    Plot.cla()
    amp = amp[20:]
    runtime = runtime[20:]
    Plot.plot(runtime, amp)
    Plot.set_xlabel('Time')
    if side=='left':
        Plot.set_ylabel('Amplitude')
    return canvas


def ft(rate, setSize, amp):
    xf = np.linspace(0.0, 1.0/(2.0*rate), setSize//2)
    freq = np.abs(np.fft.fft(amp))
    yf = freq[N//20:N//4]
    xf = xf[N//20:N//4]
    base_dB = max(yf)
    dB = []
    for v in yf:
        dB.append(10*np.log10(v/base_dB))
    return xf, dB

def build_sm(Plot, xf, yf, side='left'):
    Plot.cla()
    Plot.plot(xf, yf)
    Plot.set_xlabel('Frequency')
    if side == 'left':
        leftPlot.set_ylabel('Amplitude')

def build_sg(leftPlot, rightPlot, leftData, rightData, rate):
    leftPlot.cla()
    rightPlot.cla()
    leftData = np.array(leftData)
    rightData = np.array(rightData)
    lf, lt, lSxx = spectrogram(leftData, rate)
    rf, rt, rSxx = spectrogram(rightData, rate)
    print(len(lSxx))
    print(len(lSxx[-1]))
    leftPlot.pcolormesh(lt, lf[13:65], lSxx[13:65])
    rightPlot.pcolormesh(rt, rf[13:65], rSxx[13:65])
    leftPlot.set_xlabel('Time[Sec]')
    rightPlot.set_xlabel('Time[Sec]')
    leftPlot.set_ylabel('Frequency[Hz]')

def run():
    tTotal = time.time()
    global file
    i.set(i.get() + 1)
    if i.get() >= 4:
        i.set(0)
    if i.get() == 0:
        file = '20181111_121401467.txt'
    elif i.get() == 1:
        file = '20181111_121329107.txt'
    elif i.get() == 2:
        file = '20181111_121335497.txt'
    elif i.get() == 3:
        file = '20181111_121335587.txt'
    tAmp = time.time()
    lamp, ramp = amplitude(file)
    print("Amplitude: {}".format(time.time() - tAmp))
    if r.get() == 'sg':
        tSg = time.time()
        build_sg(leftPlot, rightPlot, lamp, ramp, fs)
        print("Spectrogram: {}".format(time.time() - tSg))
    elif r.get() == 'sm':
        tSm = time.time()
        lxf, lyf = ft(fs, N, lamp)
        rxf, ryf = ft(fs, N, ramp)
        build_sm(leftPlot, lxf, lyf, 'left')
        build_sm(rightPlot, rxf, ryf, 'right')
        print("Spectrum: {}".format(time.time() - tSm))
    elif r.get() == 'og':
        tOg = time.time()
        flamp = butter_bandpass_filter(lamp, lowcut, highcut, fs)
        framp = butter_bandpass_filter(ramp, lowcut, highcut, fs)
        build_og(leftPlot, timeArray, flamp, 'left')
        build_og(rightPlot, timeArray, framp, 'right')
        print("Oscillogram: {}".format(time.time() - tOg))
    tDraw = time.time()
    canvas.draw()
    print("canvas.draw(): {}".format(time.time() - tDraw))
    print("Total: {}".format(time.time() - tTotal))
        

def init():
    
    file = '20181111_121401467.txt'
    lamp, ramp = amplitude(file)
    fig = plt.Figure(figsize=(10, 4), dpi=80)
    leftPlot = fig.add_subplot(121)
    rightPlot = fig.add_subplot(122)
    build_sg(leftPlot, rightPlot, lamp, ramp, fs)
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().grid(row=3, column=1, columnspan=6)
    return fig, leftPlot, rightPlot, canvas, file


def plotSelection():
    lamp, ramp = amplitude(file)
    if r.get() == 'sg':
        build_sg(leftPlot, rightPlot, lamp, ramp, fs)
    if r.get() == 'sm':
        lxf, lyf = ft(fs, N, lamp)
        rxf, ryf = ft(fs, N, ramp)
        build_sm(leftPlot, lxf, lyf, 'left')
        build_sm(rightPlot, rxf, ryf, 'right')
    if r.get() == 'og':
        build_og(leftPlot, timeArray, lamp, 'left')
        build_og(rightPlot, timeArray, ramp, 'right')
    canvas.draw()
        
def waveformUpload():
    filepath = filedialog.askopenfilename(initialdir='Coding/MuellerLab/BatBotGUI/', title='Select Waveform', filetypes=(("txt files","*.txt"),("all files","*.*")))
    file = filepath
    print(file)
    return file
        
N = 10000 #number of data points
fs = 400000 # sampling frequency
lowcut = 20000 # low end of passed frequency range
highcut = 100000 # high end of passed frequency range
T = 1 / fs; #sampling interval
beginTime = 0 #start time of sampling
endTime = 0.025 #end time of sampling
timeArray = np.arange(beginTime, endTime, T); #array of sample time points
figureWidth = 40 #width of the placeholder text
figureHeight = 25 #height of the placeholder text


root = tk.Tk() #creating the program (window?)
root.title('BatBot GUI') #window title
root.geometry('800x400') #window size (size of jetson screen)
i = tk.IntVar() #iteration of startStop button
i.set(0) #start iteration at 0
r = tk.StringVar()
r.set('sg')

durationLabel = tk.Label(root, text='t: ')
sampleRateLabel = tk.Label(root, text='Fs: ')
iterationLabel = tk.Label(root, text='N: ')

drop = tk.StringVar()
drop.set('1')
iterationDrop = tk.OptionMenu(root, drop, '1','5','10','20','50','100')

fig, leftPlot, rightPlot, canvas, file = init()

startStop = tk.Button(root, text='Start/Stop', command=run)
uploadButton = tk.Button(root, text='Upload Waveform', command=waveformUpload)

sgRadio = tk.Radiobutton(root, text='Spectrogram', variable=r, value='sg', command=plotSelection)
smRadio = tk.Radiobutton(root, text='Spectrum', variable=r, value='sm', command=plotSelection)
ogRadio = tk.Radiobutton(root, text='Oscillogram', variable=r, value='og', command=plotSelection)

iterationInput = tk.Entry(root)

uploadButton.grid(row=7, column=2)
#iterationInput.grid(row=7, column=3)
startStop.grid(row=7, column=3)
sgRadio.grid(row=6, column=4)
smRadio.grid(row=7, column=4)
ogRadio.grid(row=8, column=4)
durationLabel.grid(row=7, column=5)
sampleRateLabel.grid(row=8, column=5)
iterationLabel.grid(row=6, column=5)
iterationDrop.grid(row=6, column=6)

root.mainloop()
