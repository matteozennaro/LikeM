from getdist import plots, loadMCSamples
import getdist
import numpy as np
import matplotlib.pyplot as plt
import os

LW = 1.5
FS = 16

plt.rc('axes', linewidth=LW)

plt.rc('xtick', labelsize=FS)
plt.rc('xtick.major', size=6, width=1)
plt.rc('xtick.minor', size=3, width=1)

plt.rc('ytick', labelsize=FS)
plt.rc('ytick.major', size=6, width=1)
plt.rc('ytick.minor', size=3, width=1)

font = {'family' : 'serif',
	'serif'  : 'Computer Modern Roman',
	'weight' : 'normal',
	'size'   : FS}

plt.rc('font', **font)

plt.rc('text', usetex=True)


root1 = 'chains/DEMNUNi_LCDM'
labs=['Tinker $\Lambda$CDM DEMNUni']
col1='#0000ff'

chain='plots/'
outfile = 'DEMNUNi_LCDM'

if os.path.isfile(root1+'.py_mcsamples') == True:
    os.system('rm '+root1+'.py_mcsamples')
    print '  removed file '+root1+'.py_mcsamples'

sample1 = loadMCSamples(root1)

triplot = plots.getSubplotPlotter()
triplot.settings.figure_legend_frame=False

triplot.triangle_plot(
						[sample1],
						filled=True,
						legend_labels=labs,
						legend_loc='upper right',
						contour_colors=[col1],
						contour_lws=[1.5]
						)
triplot.export(chain+outfile+'_triangle.pdf')
os.system('okular '+chain+outfile+'_triangle.pdf &')

marg = plots.getSubplotPlotter(width_inch=8)
marg.settings.figure_legend_frame=False
marg.settings.figure_legend_loc='lower right'
marg.settings.tight_layout=True
marg.settings.norm_prob_label='$\mathcal{P}$'

marg.plots_1d(
				[sample1],
				nx=3,normalized=True,
				colors=[col1],
				ls=['-'],
				legend_labels=labs,
				lws=[1.5]
				)
marg.export(chain+outfile+'_marginalized.pdf')
os.system('okular '+chain+outfile+'_marginalized.pdf &')

os.system('python /home/matteo/CosmoCodes/CosmoMC/python/GetDist.py getdist_hmf.ini')
os.system('create_tex_table tex_table.ini')
os.system('pdflatex table_'+outfile+'.tex')
os.system('mv table_'+outfile+'.pdf tmp.pdf')
os.system('rm table_'+outfile+'.*')
os.system('mv tmp.pdf table_'+outfile+'.pdf')
os.system('pdfcrop table_'+outfile+'.pdf')
os.system('mv table_'+outfile+'-crop.pdf table_'+outfile+'.pdf')
os.system('okular table_'+outfile+'.pdf &')