from __future__ import division
import pandas as pd
import math
import matplotlib as mlab
mlab.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_style("whitegrid")
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import itertools
import networkx as nx
from matplotlib_venn import venn3


from lifelines import KaplanMeierFitter
from lifelines import NelsonAalenFitter



class MY_PLOT:
	##########USING PCTA

	def violin_plt(self, df, xlab='', ylab='', tit='',color_arr=[],annotate=False, filename=''):

		plt.clf()

		#fig, ax = plt.subplots()
		plt.figure(figsize=(5,4))
		plt.title(tit)


		def corrfunc(x,y, **kws):
			#t, pv = stats.ttest_ind(x, y, nan_policy='omit')
			t, pv = stats.ranksums(x, y)
			#print pv
			ax = plt.gca()
			#plt.text(1,0,"rho = {:.2f}".format(r))

			#ax.annotate("t = {:.2f}".format(t), xy=(.8, .95), xycoords=ax.transAxes)
			#ax.annotate("p-value = {:.5f}".format(pv), xy=(.8, .9), xycoords=ax.transAxes)
			plt.annotate("p-value = {:.5f}".format(pv), xy=(0.5, 0), xycoords=('axes fraction', 'figure fraction'),xytext=(0, 10),textcoords='offset points', size=10, ha='center', va='bottom')

		df_columns = df.columns.tolist()
		if len(df_columns)==2 or annotate==True:
			corrfunc(df[df_columns[0]], df[df_columns[1]])

		#fig.set_size_inches(2,3)
		#fig.set_size_inches(4,5)
		ax = sns.violinplot(data=df,color=".8")

		sns.stripplot(data=df, palette=color_arr, jitter=True,alpha=.5, ax=ax)
		#ax = sns.violinplot(data=df,palette="hls")
		#sns.swarmplot(data=df, palette="hls", ax=ax)
		ax.set_xlabel(xlab)
		ax.set_ylabel(ylab)

		fn = ax.get_figure()
		fn.savefig(filename+'.png')
		fn.savefig(filename+'.pdf',format='PDF')
		#fn.savefig(filename+'.pdf',format='pdf')

	def single_histogram(self, df, filename ='',tit='',color='blue', xlab='', ylab=''):

		plt.clf()
		plt.ylim(0, 1.5)

		plt.title(tit)
		axes = sns.kdeplot(df, color=color, shade=False)
		axes.set_xlabel(xlab)
		axes.set_ylabel(ylab)

		fn = axes.get_figure()
		fn.savefig(filename+".png",bbox_inches='tight')
		fn.savefig(filename+".pdf",bbox_inches='tight',format='PDF')
		#fn.savefig(filename+'.pdf',bbox_inches='tight',format='pdf')

	def histogram_group(self, df, hist=False, rug=False , colo=['r','g','b','k','y'], xlab='', ylab='',tit='', cut=0, legend=False, filename=''):

		plt.clf()


		#plt.figure(figsize=(5,4))
		plt.ylim(0, 1.0)
		df_head = df.columns.tolist()

		#if df_head[0].find('PCS')>=-1:
		#	plt.ylim(0, 1.5)

		for a in range(len(df_head)):
			plt.title(tit)
			axes = sns.distplot(df[df_head[a]].dropna(), color=colo[a], hist=False, label=df_head[a])
			if legend==False:
				axes.legend_.remove()

		axes.set_xlabel(xlab)
		axes.set_ylabel(ylab)

		fn = axes.get_figure()
		fn.savefig(filename+".png",bbox_inches='tight')
		fn.savefig(filename+".pdf",bbox_inches='tight',format='PDF')
		#fn.savefig(filename+'.pdf',bbox_inches='tight',format='pdf')

	def rank_plot_sep(self,df_dat, x_dat='',y_dat='' ,tit='',xlab='',ylab='',color=[], color_numb=[], filename='',group=0):

		plt.clf()

		if group==0:
			plt.bar(x_dat, y_dat, align='center', alpha=0.5)
			plt.xlabel(xlab)
			plt.ylabel(ylab)
			plt.title(tit)
			plt.savefig(filename+".png")
			plt.savefig(filename+'.pdf',format='pdf')

		elif group==1:
			legend_arr = []
			plt.figure(figsize=(2*len(y_dat),3))

			df_temp = df_dat[[str(y_dat[0]+"samp"),str(y_dat[0])]].dropna()
			#ax = df_dat.plot(x=y_dat[0]+"samp", y=y_dat[0], color='Black',figsize=(7,3))
			plt.bar(df_temp[y_dat[0]+"samp"], df_temp[y_dat[0]], color='Black',width=1.0, label=str(y_dat[0]))

			for i, item in enumerate(y_dat[1:]):

				df_temp = df_dat[[str(item+"samp"),str(item)]].dropna()
				#df_temp.plot(x=item+"samp", y=item, color=color[i], ax=ax)
				plt.bar(df_temp[item+"samp"], df_temp[item], color=color[i],width=1.0, label=str(item))

			plt.xlabel(xlab)
			plt.ylabel(ylab)

			plt.axhline(0, color='black', linewidth=0.5)
			plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

			plt.title(tit)
			plt.savefig(filename+".png", bbox_inches='tight')
			plt.savefig(filename+".pdf", bbox_inches='tight',format='PDF')
			#plt.savefig(filename+'.pdf',format='pdf')

	def scatter_plot_and_line(self,data, x_dat,y_dat, tit='',color='blue', filename='scatter_plot', annotate = False, annotation_file='', conv_id_loc=0, conv_dat_loc=1, line_sep='\n'):

		plt.clf()

		fig = plt.figure(figsize=(5,4))

		def except_proc(arr):
			try:
				arr = [str(int(float(x))) for x in arr]
			except ValueError:
				arr = [str(x) for x in arr]
			return arr

		def corrfunc(x,y, **kws):
			r, pv = stats.spearmanr(x, y, nan_policy='omit')
			#print pv
			ax = plt.gca()
			#plt.text(1,0,"rho = {:.2f}".format(r))

			ax.annotate("rho = {:.2f}".format(r), xy=(.8, .95), xycoords=ax.transAxes)
			if pv<0.0001:
				ax.annotate("p-value < 0.0001", xy=(.8, .9), xycoords=ax.transAxes)
			else:
				ax.annotate("p-value = {:.5f}".format(pv), xy=(.8, .9), xycoords=ax.transAxes)


		test_data = data[[x_dat,y_dat]].dropna()


		plt.axhline(0, color='black', linewidth=1)
		plt.axvline(0, color='black', linewidth=1)

		ax = sns.regplot(x=x_dat,y=y_dat,line_kws={'color':'red', 'lw':1},scatter_kws={'s':20}, color=color, data=test_data)

		corrfunc(test_data[x_dat], test_data[y_dat])
		r, pv = stats.spearmanr(test_data[x_dat], test_data[y_dat], nan_policy='omit')

		#print 'corr : '+str(r)
		#print 'pval : '+str(pv)

		if annotate==True:
			f = open(annotation_file,"r")
			a_dat = f.read().strip()
			a_dat = [x.split("\t") for x in a_dat.split(line_sep)[1:]]

			conv_ids = [x[conv_id_loc] for x in a_dat]
			conv_dats = [x[conv_dat_loc] for x in a_dat]

			conv_ids = except_proc(conv_ids)
			conv_dats = except_proc(conv_dats)

			dat_index = data.index.tolist()
			dat_conv_index = []
			for x in dat_index:
				temp = [conv_dats[i] for i,item in enumerate(conv_ids) if item==x]
				dat_conv_index.append(temp[0])

			new_data = pd.DataFrame(data.values.tolist(), dat_conv_index, data.columns.tolist())
			new_data = new_data.sort_values([y_dat,x_dat], ascending=[False,True])

			with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
				print(new_data)

			#for a in range(5):
			for a in range(len(new_data.index.tolist())):
				if abs(new_data[x_dat][a]/new_data[y_dat][a]) > 1.05:
					ax.annotate(new_data.index[a], (new_data[x_dat][a], new_data[y_dat][a]),xytext=(15, 50), textcoords='offset points', arrowprops=dict(arrowstyle='-|>', connectionstyle='angle3, angleA=0, angleB=90'))

		ax.set_title(tit)
		fn = ax.get_figure()
		fn.savefig('%s.png'%(filename))
		fn.savefig('%s.pdf'%(filename),format='PDF')


	def scatter_matr(self, df, type=1, color='', filename=''):

		plt.clf()

		if type==1:
			def corrfunc(x,y, **kws):
				r, pv = stats.spearmanr(x, y)
				ax = plt.gca()
				ax.annotate("rho = {:.2f}".format(r), xy=(.7, .95), xycoords=ax.transAxes)
				ax.annotate("p-value = {:.2f}".format(pv), xy=(.7, .9), xycoords=ax.transAxes)

			#titles = ['androgen_response', 'androgen_biosynthesis', 'cholesterol_homeostasis', 'cholesterol_biosynthesis','steroid_biosynthesis']
			#df = pd.DataFrame(arr1, columns=titles)

			axes = sns.PairGrid(df.dropna(), palette=color,size=3.5)
			axes.map_upper(sns.regplot,line_kws={'color': 'red'}, scatter_kws={'s':10})
			axes.map_upper(corrfunc)
			axes.map_diag(sns.distplot, kde="True")
			axes.map_lower(sns.regplot,scatter_kws={'s':10},fit_reg=False)
			axes.savefig(filename+".png",bbox_inches='tight')
			axes.savefig(filename+".pdf",bbox_inches='tight',format='pdf')
		#elif type==2:

	def shared_scatter(self, df_arr, x_dat='',y_dat='',color_arr=[],title=[],filename=''):

		color_arr = ['magenta'] + color_arr
		def corrfunc(x,y,ax, **kws):
			r, pv = stats.spearmanr(x, y, nan_policy='omit')
			#print pv
			#plt.text(1,0,"rho = {:.2f}".format(r))

			#ax.annotate("rho = {:.3f}".format(r), xy=(.75, .95), xycoords=ax.transAxes)
			#ax.annotate("p-value = {:.3f}".format(pv), xy=(.75, .9), xycoords=ax.transAxes)
			#ax.annotate("rho = {:.3f}".format(r), xy=(.75, .95), xycoords=ax.transAxes)
			if pv<0.001:
				ax.annotate("rho = {:.3f}, p-value < 0.001".format(r, pv), xy=(0.5, 0), xytext=(0, -0.5), xycoords=('axes fraction', 'figure fraction'), textcoords='offset points', size=12, ha='center', va='bottom')
			else:
				ax.annotate("rho = {:.3f}, p-value = {:.3f}".format(r, pv), xy=(0.5, 0), xytext=(0, -0.5), xycoords=('axes fraction', 'figure fraction'), textcoords='offset points', size=12, ha='center', va='bottom')


		subplot_numb = len(df_arr)
		fig = plt.figure(figsize=(4.25*subplot_numb,3.5))
		subplot_cont = []

		for a in range(subplot_numb):
			temp = fig.add_subplot(1,subplot_numb,a+1)
			subplot_cont.append(temp)

		for i, item in enumerate(df_arr):
			test_data = item[[x_dat,y_dat]].dropna()
			subplot_cont[i].set_title(title[i])

			corrfunc(test_data[x_dat], test_data[y_dat], ax=subplot_cont[i])
			sns.regplot(x=x_dat,y=y_dat,line_kws={'color':'red', 'lw':1}, color=color_arr[i], scatter_kws={'s':20}, data=item, ax=subplot_cont[i])

		plt.tight_layout()
		fig.savefig('%s.png'%(filename))
		fig.savefig('%s.pdf'%(filename),format='PDF')

	def network_plot(self, df, expr, tit='',filename=''):
		plt.clf()

		def expr_counting(expr):
			dat = expr
			up_dat = dat[dat>0.28].count()
			down_dat = dat[dat<-0.28].count()
			zero_dat = dat[(dat<0.28) & (dat>-0.28)].count()

			return up_dat, down_dat, zero_dat

		def value_min_max(v):
			if v>2.0:
				return 2.0
			elif v<-2.0:
				return -2.0
			else:
				return v

		G = nx.Graph()

		width_list = []
		#corr_list = []
		for subset in itertools.combinations(df.index.tolist(), 2):
			G.add_node(subset[0])
			G.add_node(subset[1])

			s1_up, s1_down, s1_zero = expr_counting(expr.loc[subset[0]])
			s2_up, s2_down, s2_zero = expr_counting(expr.loc[subset[1]])

			chi2, p, dof, expected = stats.chi2_contingency([[s1_up,s1_zero,s1_down], [s2_up, s2_zero, s2_down]])
			if p < 0.01:
				G.add_edge(subset[0],subset[1])
				width_list.append(-np.log(p))

		width_med = np.median(width_list)
		width_list = [m/width_med for m in width_list]

		size_map = df['mapped'].to_dict()
		#val_map = expr.loc[df.index.tolist()].median(axis=1).to_dict()
		val_map = df['fc'].to_dict()
		default_node_size = 50

		#cmap = sns.diverging_palette(130, 10, s=99, l=50,sep=21,center='dark', as_cmap=True)

		values = [val_map.get(node, 0.0) for node in G.nodes()]
		#values = [0 if v<0.28 and v>-0.28 else v for v in values]

		values = [ value_min_max(v) for v in values]
		sizes = [size_map.get(node, default_node_size) for node in G.nodes()]
		sizes = [s*default_node_size for s in sizes]

		plt.title('Samples:'+tit)

		#nb_count = [len(list(G.neighbors(node))) for node in G.nodes()]
		#print nb_count

		pos = nx.spring_layout(G)
		#nx.draw(G,pos, cmap=cmap, node_color=values, node_size=sizes, with_labels=True, width=width_list, vmin=-1, vmax=1)
		nx.draw(G,pos, cmap='bwr', node_color=values, node_size=sizes, with_labels=True, width=width_list, vmin=-0.5, vmax=0.5)

		ax = plt.gca() # to get the current axis
		if values!=[]:
			ax.collections[0].set_edgecolor("black")

		plt.savefig(filename+".png", format="PNG")
		plt.savefig(filename+".pdf", format="PDF")

	def lollipop(self, dat, label_data=[],color_arr=[] ,ylab='',tit='', filename=''):
		plt.clf()
		plt.figure(figsize=(8,3))
		print ylab
		color=color_arr
		counter = 0
		for i, item in enumerate(dat):
			markerline, stemlines, baseline = plt.stem(range(counter,counter+len(item)), item, color[i], markerfmt=' ',label=label_data[i])
			#plt.setp(stemlines, 'color', plt.getp(markerline,'color'))
			#plt.setp(stemlines, 'linestyle', 'dotted')
			plt.setp(baseline, 'color', color[i], 'linewidth', 2)

			counter = counter+len(item)

		plt.legend()
		plt.xlabel('Sorted Samples')
		plt.ylabel(ylab)
		plt.title(tit)
		plt.savefig(filename+".png", bbox_inches='tight')
		plt.savefig(filename+".pdf", bbox_inches='tight',format='PDF')

	def survival_plot_and_cox(self, df_arr, label=[], filename=''):
		plt.clf()
		color = ['red','green','blue','cyan','orange','black']

		kmf = KaplanMeierFitter()
		naf = NelsonAalenFitter()

		for a in range(len(df_arr)):
			df_el = df_arr[a]
			if a==0:
				kmf.fit(df_el['bcrmonth'],df_el['bcrstatus'],label=label[a])
				ax = kmf.plot(show_censors=True, ci_show=False, color=color[a], ylim=(0,1))
			else:
				kmf.fit(df_el['bcrmonth'],df_el['bcrstatus'],label=label[a])
				kmf.plot(ax=ax, show_censors=True, ci_show=False, color=color[a], ylim=(0,1))

		fig = ax.get_figure()
		fig.savefig(filename+'.png')
		fig.savefig(filename+'.pdf',format='PDF')

	##########USING PCTA

	def __init__(self):
		print "PLOT ON!"
