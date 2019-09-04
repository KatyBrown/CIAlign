#! /usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utilityFunctions
import math


def arrNumeric(arr, typ):
    arr = np.flip(arr, axis=0)
    if typ == 'nt':
        D = utilityFunctions.getNtColours()
    else:
        D = utilityFunctions.getAAColours()
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

def drawMiniAlignment(arr, nams, log, outfile, typ,
                      dpi, title, width, height, markup=False,
                      markupdict=None):
    fontsize = 1500 / dpi
    ali_height, ali_width = np.shape(arr)
    om = math.floor(math.log10(ali_height))
    tickint = 1 if om == 0 else 10 if om == 1 else 100
    lineweight_h = 10 / ali_height
    lineweight_v = 10 / ali_width
    f = plt.figure(figsize=(width, height), dpi=dpi)
    a = f.add_subplot('111')
    a.set_xlim(-0.5, ali_width)
    a.set_ylim(-0.5, ali_height-0.5)
    arr2, cm = arrNumeric(arr, typ)
    a.imshow(arr2, cmap=cm, aspect='auto')
    a.hlines(np.arange(-0.5, ali_height), -0.5,
             ali_width, lw=lineweight_h, color='white',
             zorder=100)
    a.vlines(np.arange(-0.5, ali_width), -0.5,
             ali_height, lw=lineweight_v, color='white',
             zorder=100)
    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.spines['left'].set_visible(False)
   # a.tick_params(left='off', labelleft='on')
    if title:
        f.suptitle(title)
    for t in a.get_xticklabels():
        t.set_fontsize(fontsize)
    a.set_yticks(np.arange(ali_height-1, -1, -tickint))
    if tickint == 1:
        a.set_yticklabels(np.arange(1, ali_height+1, tickint), fontsize=fontsize)
    else:
        a.set_yticklabels(np.arange(0, ali_height, tickint), fontsize=fontsize)

    if markup:
        if "crop_ends" in markupdict:
            colour = "black"
            for nam, boundary in markupdict['crop_ends'].items():
                y = len(nams) - nams.index(nam) - 1
                
                if boundary[0].shape[0] > 0:
                    a.add_patch(matplotlib.patches.Rectangle([boundary[0][0]-0.5, y-0.5],
                                                             (boundary[0][-1] - boundary[0][0]) + 1,
                                                             1, color='black', lw=0, zorder=50))


                if boundary[1].shape[0] > 0:
                    a.add_patch(matplotlib.patches.Rectangle([boundary[1][0]-0.5, y-0.5],
                                                             (boundary[1][-1] - boundary[1][0]) + 1,
                                                             1, color='black', lw=0, zorder=50))
        
        if "remove_badlyaligned" in markupdict:
            colour = '#f434c5'
            for row in markupdict['remove_badlyaligned']:
                y = len(nams) - nams.index(row) - 1.5
                a.add_patch(matplotlib.patches.Rectangle((-0.5, y), ali_width, 1, color=colour, zorder=49, lw=0))
                
        if "remove_insertions" in markupdict:
            colour = "#7bc5ff"
            for col in markupdict['remove_insertions']:
                a.add_patch(matplotlib.patches.Rectangle((col-0.5, -0.5), 1,
                                                         ali_height,
                                                         color=colour,
                                                         zorder=48, lw=0))
   
        if "remove_short" in markupdict:
            colour = '#fff6b3'
            for row in markupdict['remove_short']:
                y = len(nams) - nams.index(row) - 1.5
                a.add_patch(matplotlib.patches.Rectangle((-0.5, y), ali_width, 1, color=colour, zorder=47, lw=0))


        if "remove_gaponly" in markupdict:
            colour = "#f57700"
            for col in markupdict['remove_gaponly']:
                a.add_patch(matplotlib.patches.Rectangle((col-0.5, -0.5), 1,
                                                         ali_height,
                                                         color=colour,
                                                         zorder=46, lw=0))
        legend = plt.figure(figsize=(2, 2), dpi=100)
        l = legend.add_subplot(111)
        colours = ['black', '#f434c5', "#7bc5ff", '#fff6b3', "#f57700"]
        functions = ['Cropped Ends', 'Badly Aligned', 'Insertions', 'Too Short', 'Gap Only']
        for i, c in enumerate(colours):
            l.plot(1, 5-i, marker='.', color=c, markersize=20)
            l.text(2, 5-i, functions[i])
        l.set_xlim(0.5, 3)
        l.set_ylim(-1, 6)
        l.set_axis_off()
        legend.gca().set_axis_off()
        l.margins(0, 0)
        legend.savefig("%s_legend.png" % (outfile.replace(".png", "")),
                       dpi=100, bbox_inches='tight')
        
    f.savefig(outfile, dpi=dpi)
