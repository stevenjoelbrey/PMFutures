fig = plt.figure(figsize=(33,22))
nVar = len(v.keys())
#nVar = 2
frameN = 0

vmin = 10**-9 * 1.
vmax = np.percentile(E, 99)

for i in range(nVar):

	x = v.keys()[i]

	for j in range(nVar):
	
		y = v.keys()[j]

		frameN += 1 

		ax = fig.add_subplot(5,5,frameN, frame_on=False)

		c = ax.scatter(v[x], v[y], c = E,\
	                    marker=".",\
                        s=2,\
	                    edgecolors='none',\
	                    vmin=vmin,vmax=vmax,\
	                    alpha=0.7,\
					    norm=matplotlib.colors.LogNorm()
	                    )
	
		ax.spines['top'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(axis='y', labelsize=15)
		ax.tick_params(axis='x', labelsize=15)


		#cb = plt.colorbar(c, cmap='jet', extend='both')
		#cb.set_label('Emissions [kg/day]')

		plt.xlim([v[x].min(), v[x].max()])
		plt.ylim([v[y].min(), v[y].max()])

		plt.xlabel(u[x], fontsize=20)
		plt.ylabel(u[y], fontsize=20)

fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(c, cax=cbar_ax, cmap='jet', extend='both')
cbar.set_label('Emissions [kg m$^{-2}$ day$^{-1}$]', fontsize=30)

plt.savefig('../Figures/' + AScenario + '_parameter_space.png', dpi=500)

