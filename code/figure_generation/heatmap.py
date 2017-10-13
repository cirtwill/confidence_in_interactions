
# This example demonstrates how to use a color plot.
#
# One important feature of color plots is that drawing objects in
# XMGrace are plotted above the frame of a graph.  As a result, we use
# the Grace.clone_graph method to copy the frame of the graph and make
# things look pretty.
#
# Another important feature of color plots is that specifying a
# particular world coordinate can automatically place drawing objects
# outside of the domain of the graph.
# Graph.remove_extraworld_drawing_objects fixes this problem by
# removing all of the drawing objects that appear in the world
# coordinates and are out of bounds.
#
# It also demonstrates how to use a PyGrace.Style to format figures in
# a customized format

#import math
#from PyGrace.Extensions.panel import Panel,MultiPanelGrace



import random
import sys
import os

from PyGrace.grace import Grace as Plot
from PyGrace.Extensions.colorbar import SolidRectangle, ColorBar
from PyGrace.colors import ColorBrewerScheme
from PyGrace.axis import LINEAR_SCALE, LOGARITHMIC_SCALE
from PyGrace.drawing_objects import DrawBox, DrawText, DrawLine

from PyGrace.Styles.el import ElGraph, ElLogColorBar, ElLinColorBar


#------------------------------------------------------------------------------
# make a nice figure
#------------------------------------------------------------------------------
# get the data (as specified in the command line)
readfile = sys.argv[1]
filename = readfile+'.csv'
troph1 = sys.argv[2]
tr1id = troph1+'_ID'
troph2 = sys.argv[3]
tr2id = troph2+'_ID'
infile = open(filename)
header=infile.readline().strip().split(',')


# read in the data into something we can manipulate easily
data = []
for line in infile:
	sline = line.strip().split(',')
	pt = dict(zip(header,[i for i in sline]))
	data.append(pt)
infile.close()
print data[0:7]

ndata = {}
for d in data:
	if float(d['interact']) == 0:
		if float(d['cooccur']) == 0:
			try:
				ndata[(float(d[tr1id]),float(d[tr2id]))].append(0)
			except KeyError:
				ndata[(float(d[tr1id]),float(d[tr2id]))] = 0
		else:
			try:
				ndata[(float(d[tr1id]),float(d[tr2id]))].append(-float(d['cooccur']))
			except KeyError:
				ndata[(float(d[tr1id]),float(d[tr2id]))] = -float(d['cooccur'])
	else:
		try:
			ndata[(float(d[tr1id]),float(d[tr2id]))].append(float(d['interact']))
		except KeyError:
			ndata[(float(d[tr1id]),float(d[tr2id]))] = float(d['interact'])


print ndata[(63,17)]
print ndata[(63,5)]
print ndata[(63,34)]
print ndata[(63,18)]


# instantiate a sweet figgy fig
colors = ColorBrewerScheme("PRGn",n=253) # this is the maximum number of colors

# you can change the opacity percent of a colorscheme if you want:
# colors.change_opacity(20, exclude_black=False)

grace = Plot(colors=colors)

zmin = min(ndata.values())
zmax = max(ndata.values())
zlim = max(-zmin,zmax)
print zlim


# add a colorbar
colorbar = grace.add_graph(ElLogColorBar,domain=(-zlim,zlim),
                           scale=LINEAR_SCALE,autoscale=False)

# to add some data to graph, just add SolidRectangle datasets
graph = grace.add_graph()
graph.copy_format(ElGraph)


for d in ndata:
	color = colorbar.z2color(ndata[d])

	graph.add_dataset([(d[0]-0.5,d[1]-0.5),(d[0]+0.5,d[1]+0.5)], SolidRectangle, color)


# for (x0,y0,x1,y1,pdf) in data:
    # if pdf > 0.0:
    #     color = colorbar.z2color(pdf)
    #     # you can change the opacity percentage of a single color, as well
    #     # color.change_opacity(60)
    #     graph.add_dataset([(x0,y0), (x1,y1)], SolidRectangle, color)

# move things around
graph.set_view(0.2,0.2,0.9,0.9)
colorbar.set_view(0.95,0.2,1.0,0.9)

# label axes
colorbar.set_label("Interaction/cooccurrence")
graph.xaxis.label.text = troph1
graph.yaxis.label.text = troph2

# all of the autoscale features work, too
graph.autoscalex()
graph.autotickx()

# restrict y-axis scale to the same as the x-axis scale
xmin,ymin,xmax,ymax = graph.get_world()
graph.set_world(xmin,xmin,xmax,xmax)
graph.autoticky()

# get rid of points that are out of bounds
graph.remove_extraworld_drawing_objects()

# print out the grace
outfile=readfile+".jpg"
grace.write_file('plot1.eps')
