# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 11:44:27 2016

@author: yungkuo
"""

import numpy as np
import os

def listdir(path, startswith, endswith):
    for f in os.listdir(path):
        if f.startswith(startswith):
            if f.endswith(endswith):
                yield f

def findparticle(image, boundary, data, scan, nstd):
    pts = []
    for i in np.arange(scan[0],data[1]-scan[0],1, dtype='int'):
        for j in np.arange(boundary[0],boundary[1],1, dtype='int'):
            if image[i,j] == np.max(image[(i-scan[0]):(i+scan[0]),(j-scan[1]):(j+scan[1])]):
                if np.mean(image[(i-scan[0]):(i+scan[0]),(j-scan[1]):(j+scan[1])]) > np.mean(image[(i-scan[0]):(i+scan[0]),boundary[0]:boundary[1]])+nstd*np.std(image[(i-scan[0]):(i+scan[0]),boundary[0]:boundary[1]], ddof=1):
                    if np.mean(image[(i-scan[0]):(i+scan[0]),(j-scan[1]):(j+scan[1])])> np.mean(image[(i-scan[0]*2):(i+scan[0]*2),(j-scan[1]*2):(j+scan[1]*2)]):
                        pt = [j,i]
                        pts = np.append(pts, pt)
    return np.reshape(pts,[len(pts)/2,2]).astype('int')


def findparticle_round(image, boundary, data, scan, nstd):
    pts = []
    for i in np.arange(scan[0],data[1]-scan[0],1, dtype='int'):
        for j in np.arange(boundary[0],boundary[1],1, dtype='int'):
            if image[i,j] == np.max(image[(i-scan[0]):(i+scan[0]),(j-scan[1]):(j+scan[1])]):
                if np.mean(image[(i-scan[0]):(i+scan[0]),(j-scan[1]):(j+scan[1])]) > np.mean(image[:,boundary[0]:boundary[1]])+nstd*np.std(image[:,boundary[0]:boundary[1]], ddof=1):
                    pt = [j,i]
                    pts = np.append(pts, pt)
    return np.reshape(pts,[len(pts)/2,2])


def get_roi_square(point, pad):
    """Return a square selection of pixels around `point`.
    """
    col, row = point
    mask = (slice(row-pad[0], row+pad[0]+1), slice(col-pad[1], col+pad[1]+1))
    return mask


def get_roi_square_3d(point, pad):
    """Return a square selection of pixels around `point`.
    """
    col, row = point
    mask = (slice(None, None), slice(row-pad[0], row+pad[0]+1), slice(col-pad[1], col+pad[1]+1))
    return mask

def get_timetrace_square(video, point, pad):
    """Returna a timetrace from `video` by averaging a square pixel selection.
    """
    mask = get_roi_square_3d(point, pad)
    timetrace = video[mask].mean(1).mean(1)
    return timetrace


def get_timetrace(movie, pt, scan, return_mean=False):
    timetrace = np.mean(np.mean(movie[:,pt[1]-scan[0]:pt[1]+scan[0],pt[0]-scan[1]:pt[0]+scan[1]], axis=1), axis=1)
    timetrace_mean = timetrace.mean()
    if return_mean:
        return timetrace, timetrace_mean
    else:
        return timetrace

def get_timetrace_sum(movie, pt, scan, return_mean=False):
    timetrace = np.sum(np.sum(movie[:,pt[1]-scan[0]:pt[1]+scan[0],pt[0]-scan[1]:pt[0]+scan[1]], axis=1), axis=1)
    timetrace_mean = timetrace.mean()
    if return_mean:
        return timetrace, timetrace_mean
    else:
        return timetrace