#! /usr/bin/env python
import logging
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import consensusSeq

def arrNumeric(arr, typ):
    if typ == 'nt':
        D = consensusSeq.getNtColours()
    else:
        D = consensusSeq.getAAColours()
    
    keys = list(D.keys())
    colours = list(D.values())
    
    ali_height, ali_width = np.shape(arr)
    
    D = {x: i for i, x in enumerate(keys)}
    
    arr2 = np.empty([ali_height, ali_width])
    for x in range(ali_width):
        for y in range(ali_height):
            arr2[y, x] = D[arr[y, x]]
    cmap = matplotlib.colors.ListedColormap(colours)
    return (arr2, cmap)

def drawMiniAlignment(arr,
                      nams,
                      log,
                      outfile,
                      typ,
                      dpi,
                      title,
                      width,
                      height,
                      markup=False,
                      markupdict=None):

    ali_height, ali_width = np.shape(arr)
    lineweight_h = 10 / ali_height
    lineweight_v = 10 / ali_width
    if  typ == 'nt':
        cD = consensusSeq.getNtColours()
    elif typ == 'aa':
        cD = consensusSeq.getAAColours()
    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot('111')
    a.set_xlim(0, ali_width)
    a.set_ylim(-0.5, ali_height-0.5)
    xpoints = range(0, ali_width)
    arr2, cm = arrNumeric(arr, typ)
    a.imshow(arr2, cmap=cm, aspect='auto')
    a.hlines(np.arange(-0.5, ali_height), -0.5, ali_width, lw=lineweight_h, color='white')
    a.vlines(np.arange(-0.5, ali_width), -0.5, ali_height, lw=lineweight_v, color='white')


    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.spines['left'].set_visible(False)
    a.yaxis.set_ticks_position('none')
    a.yaxis.set_visible(False)

    if title:
        f.suptitle(title)
    for t in a.get_xticklabels():
        t.set_fontsize(8)

    if markup:
        # print(markupdict)
        if "crop_ends" in markupdict:
            colour = "#8470FF"
            for nam, boundary in markupdict['crop_ends'].items():
                i = nams.index(nam) + 1
                if boundary[0].shape[0] > 0:
                    a.hlines((i - 0.5), boundary[0][0], boundary[0][-1] + 1, color=colour, lw=lineweight, zorder=0)
                if boundary[1].shape[0] > 0:
                    a.hlines((i - 0.5), boundary[1][0], boundary[1][-1] + 1, color=colour, lw=lineweight, zorder=0)
        
        if "remove_short" in markupdict:
            colour = "#d7ddf2"
            for nam in markupdict['remove_short']:
                i = nams.index(nam) + 1
                a.hlines((i - 0.5), 0, ali_width, color=colour, lw=lineweight, zorder=0)

        if "remove_insertions" in markupdict:
            colour = "#fece88"
            a.vlines(list(markupdict['remove_insertions']),
                     0, ali_height, color=colour, lw=lineweight, zorder=1)
        if "remove_gaponly" in markupdict:
            colour = "#EEAEEE"
            a.vlines(list(markupdict['remove_gaponly']),
                     0, ali_height, color=colour, lw=0.75, zorder=1)
    f.savefig(outfile, dpi=dpi)
