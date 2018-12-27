from celery import shared_task,current_task, task
from PLOTS import MY_PLOT
import numpy as np
import pandas as pd
import gseapy
import glob
from scipy import stats

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from multiprocessing import Process, Queue
from pcta_updated import all_expr_df, pcta_id, mra_set, tcga_expr_df

import json
import os
import ast
import shutil
import os.path
from numpy import random
from time import sleep
from celery import Celery
import itertools

def float_format_change(fv):
	if fv < 0.001:
		new_fv = '<0.001'
	else:
		new_fv = str('{:.3f}'.format(fv))
	return new_fv

def test_and_others_selection(test_name, group_info, group_data, others_name=[]):
	if others_name==[]:
		others_index = [ i for i,item in enumerate(group_info) if item!=test_name]
		others_values = [ group_data[i] for i in others_index ]
		others_values = list(sum(others_values, []))
		test_values = group_data[group_info.index(test_name)]

		return test_values, others_values
	else:
		others_index = [ i for i,item in enumerate(group_info) if item in others_name]
		others_values = [ group_data[i] for i in others_index ]
		others_values = list(sum(others_values, []))
		test_values = group_data[group_info.index(test_name)]

		return test_values, others_values

def gene_set_zscore_single_thr( arr1, gene_set=[] ,sample_status="multiple"):

	def cal_z(temp_df, inter):
		selected = temp_df.loc[inter].tolist()
		selected_adj = [float(x) for x in selected]
		total = temp_df.tolist()
		total_adj = [float(x) for x in total]

		diff_mean = np.nanmean(selected_adj) - np.nanmean(total_adj)
		result = diff_mean*np.sqrt(len(selected_adj))/np.nanstd(total_adj, ddof=1)
		return result

	ft_dat_index = gene_set
	zscore=[]

	arr1_index = arr1.index.tolist()
	inter = list(set(arr1_index).intersection(ft_dat_index))

	if sample_status=="single":
		zscore = [cal_z(arr1, inter)]

	elif sample_status=="multiple":
		diff_mean = arr1.loc[inter].mean(axis=0).subtract(arr1.mean(axis=0))
		len_norm = arr1.std(ddof=1, axis=0).apply(lambda x: np.sqrt(len(inter))/x)
		zscore = diff_mean*len_norm
		zscore = zscore.to_frame()
		zscore.columns = ['Zscore']
		#zscore = [cal_z(arr1[x], inter) for x in arr1.columns.tolist()]
		zscore = zscore['Zscore'].values.tolist()

	return zscore

def table_split(cph):
		cph_d = cph.summary
		cph_d_index = cph_d.index.tolist()

		if len(cph_d_index)>1:
			t_arr = []

			for x in cph_d_index:

				testing = x
				coef = format(cph_d['coef'].loc[x],".2f")
				haz_c = format(cph_d['exp(coef)'].loc[x],".2f")
				haz_c_int = str(format(np.exp(cph_d['lower 0.95'].loc[x]),".2f"))+"-"+str(format(np.exp(cph_d['upper 0.95'].loc[x]),".2f"))
				pval = format(cph_d['p'].loc[x],".3f")
				concor = format(concordance_index(cph.durations, -cph.predict_partial_hazard(cph.data).values.ravel(), cph.event_observed),".2f")

				if str(pval)=='0.000':
					pval = '<0.001'
				else:
					pval = str(pval)
				t_arr.append(["Multivariate",testing, coef, haz_c+"("+haz_c_int+")", pval, concor])

			return t_arr
		else:
			for x in cph_d_index:
				testing = x
				coef = format(cph_d['coef'].loc[x],".2f")
				haz_c = format(cph_d['exp(coef)'].loc[x],".2f")
				haz_c_int = str(format(np.exp(cph_d['lower 0.95'].loc[x]),".2f"))+"-"+str(format(np.exp(cph_d['upper 0.95'].loc[x]),".2f"))
				pval = format(cph_d['p'].loc[x],".3f")
				concor = format(concordance_index(cph.durations, -cph.predict_partial_hazard(cph.data).values.ravel(), cph.event_observed),".2f")

				if str(pval)=='0.000':
					pval = '<0.001'
				else:
					pval = str(pval)
					
			return [["Univariate", testing, coef, haz_c+"("+haz_c_int+")", pval, concor]]


@shared_task
def association_c(input_data, dataset_input_data, group_info, group_samples, plot_color, session_key):
	current_task.update_state(state='PROGRESS', meta={'process_percent': 0})

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+session_key+"/"
	nginx_plot_path= "/home/ubuntu/django_proj/pcta_updated/main/staticimages/"+session_key+"/"

	template_plot_path = "images/"+session_key+"/"
	test_name = "association_test"

	if not os.path.exists(plot_path):
		os.mkdir(plot_path)
		os.mkdir(nginx_plot_path)
	else:
		files = glob.glob(plot_path+"*")
		for f in files:
			os.remove(f)
		files = glob.glob(nginx_plot_path+"*")
		for f in files:
			os.remove(f)

	filecount = 0
	file_list = []

	#####Lollipop
	pt=MY_PLOT()

	if dataset_input_data=='PCTA':
		df = all_expr_df
		df = df.drop('Symbol',axis=1)
	else:
		df = tcga_expr_df

	current_task.update_state(state='PROGRESS', meta={'process_percent': 10})

	if len(input_data)>1:
		df_data = [gene_set_zscore_single_thr(df[samples], input_data, sample_status="multiple") for samples in group_samples]
		ylab = 'Z score'
	elif len(input_data)==1:
		df_data = [df[samples].loc[input_data[0]].values.tolist() for samples in group_samples]
		ylab = 'Expression'

	current_task.update_state(state='PROGRESS', meta={'process_percent': 20})

	pt.lollipop([sorted(d) for d in df_data], label_data=group_info, ylab=ylab, color_arr=plot_color,filename=plot_path+test_name+str(filecount))

	###########Oneway ANNOVA with group data
	if len(df_data)==3:
		fv,pv = stats.f_oneway(df_data[0],df_data[1], df_data[2])
	else:
		fv,pv = stats.f_oneway(df_data[0],df_data[1], df_data[2], df_data[3])

	rw = open(plot_path+'onewayANOVA_result.tsv',"w")
	rw.write("F-value\tP-value\n")
	rw.write(str('{:.3f}'.format(fv))+'\t'+float_format_change(pv)+'\n')
	rw.close()
	###########Oneway ANNOVA with group data

	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1
	#####Lollipop
	current_task.update_state(state='PROGRESS', meta={'process_percent': 30})

	#####Violin or boxplot
	df_data_arr = [ pd.DataFrame(data=d) for d in df_data]
	df_data_arr = pd.concat(df_data_arr, axis=1)
	df_data_arr.columns = group_info

	pt.violin_plt(df_data_arr,ylab=ylab,tit='',color_arr=plot_color, filename=plot_path+test_name+str(filecount))

	###########Ranksum with group data
	rw = open(plot_path+'ranksum_result.tsv',"w")
	rw.write("Test group\tFoldChange\tP-value\n")

	if 'GS=7' in group_info:
		#others_index = [ i for i,item in enumerate(group_info) if item!='mCRPC' ]
		#others_values = [ df_data[i] for i in others_index ]
		#others_values = list(sum(others_values, []))
		#test_values = df_data[group_info.index('mCRPC')]
		if dataset_input_data=='PCTA':
			test_values, others_values = test_and_others_selection('mCRPC', group_info, df_data, others_name=['GS=7','GS<7','GS>7'])
			rv,pv = stats.ranksums(test_values, others_values)
			fc = float(np.mean(test_values)) - float(np.mean(others_values))
			rw.write("mCRPC VS Primary"+'\t'+str('{:.3f}'.format(fc))+'\t'+float_format_change(pv)+'\n')

			test_values, others_values = test_and_others_selection('Benign', group_info, df_data, others_name=['GS=7','GS<7','GS>7'])
			rv,pv = stats.ranksums(test_values, others_values)
			fc = float(np.mean(others_values)) - float(np.mean(test_values))
			rw.write("Primary VS Benign"+'\t'+str('{:.3f}'.format(fc))+'\t'+float_format_change(pv)+'\n')
			rw.close()

		else:
			test_values, others_values = test_and_others_selection('GS>7', group_info, df_data, others_name=['GS=7','GS<7'])
			rv,pv = stats.ranksums(test_values, others_values)
			fc = float(np.mean(test_values)) - float(np.mean(others_values))
			rw.write("GS>7 VS Others"+'\t'+str('{:.3f}'.format(fc))+'\t'+float_format_change(pv)+'\n')

			test_values, others_values = test_and_others_selection('GS<7', group_info, df_data, others_name=['GS=7','GS>7'])
			rv,pv = stats.ranksums(test_values, others_values)
			fc = float(np.mean(others_values)) - float(np.mean(test_values))
			rw.write("GS<7  VS Others"+'\t'+str('{:.3f}'.format(fc))+'\t'+float_format_change(pv)+'\n')
			rw.close()

		#ranksum_test_case = [df_data[group_info.index('mCRPC')], ]
	elif 'PCS1' in group_info:
		test_values, others_values = test_and_others_selection('PCS1', group_info, df_data)
		rv,pv = stats.ranksums(test_values, others_values)
		fc = float(np.mean(test_values)) - float(np.mean(others_values))
		rw.write("PCS1 VS Others"+'\t'+str('{:.3f}'.format(fc))+'\t'+float_format_change(pv)+'\n')
		rw.close()

	elif 'LumB' in group_info:
		test_values, others_values = test_and_others_selection('LumB', group_info, df_data)
		rv,pv = stats.ranksums(test_values, others_values)
		fc = float(np.mean(test_values)) - float(np.mean(others_values))
		rw.write("LumB VS Others"+'\t'+str('{:.3f}'.format(fc))+'\t'+float_format_change(pv)+'\n')
		rw.close()


	current_task.update_state(state='PROGRESS', meta={'process_percent': 40})
	###########Ranksum with group data

	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1
	#####Violin or boxplot

	current_task.update_state(state='PROGRESS', meta={'process_percent': 50})

	#####Histogram
	#pt.histogram_group(df_data_arr,legend=True,colo=plot_color,filename=plot_path+test_name+str(filecount))
	pt.line_trend(df_data_arr,ylab='Mean of '+ylab,legend=True,colo=plot_color,filename=plot_path+test_name+str(filecount))
	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1

	os.system("cp -rf %s/* %s"%(plot_path,nginx_plot_path))
	current_task.update_state(state='PROGRESS', meta={'process_percent': 100})
	#####Histogram

	#return file_list
	return random.random()

@shared_task
def survival_c(input_data, session_key):
	current_task.update_state(state='PROGRESS', meta={'process_percent1': 0})
	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+session_key+"/"
	nginx_plot_path= "/home/ubuntu/django_proj/pcta_updated/main/staticimages/"+session_key+"/"

	template_plot_path = "images/"+session_key+"/"
	test_name = "bcr_analysis"

	if not os.path.exists(plot_path):
		os.mkdir(plot_path)
		os.mkdir(nginx_plot_path)
	else:
		files = glob.glob(plot_path+"*")
		for f in files:
			os.remove(f)
		files = glob.glob(nginx_plot_path+"*")
		for f in files:
			os.remove(f)

	filecount = 0
	file_list = []

	pt=MY_PLOT()

	surv_dataset = ['TCGA_PRAD', 'GSE40272', 'GSE70769']
	expr_set_arr = []
	clinical_set_arr = []


	for surv in  surv_dataset:
		if surv!='TCGA_PRAD':
			expr_set = pd.read_csv('user_data/clinical_data/'+surv+'.csv',index_col=0)
			expr_set.index = expr_set.index.astype(int)
			expr_set.index = expr_set.index.astype(str)
			clinical_set = pd.read_excel('user_data/clinical_data/'+surv+'_clinical.xlsx',index_col=0)

		else:
			expr_set = tcga_expr_df
			clinical_set = pd.read_excel('user_data/clinical_data/'+surv+'_clinical.xlsx',index_col=0)

		expr_set_arr.append(expr_set)
		clinical_set_arr.append(clinical_set)

	current_task.update_state(state='PROGRESS', meta={'process_percent1': 20})

	test_set_arr = []
	for expr in expr_set_arr:
		if len(input_data)>1:
			df_data = gene_set_zscore_single_thr(expr, input_data, sample_status="multiple")
			test_set = pd.Series(data=df_data, index=expr.columns.tolist(), name='Input_Zscore')
			test_set_arr.append(test_set)
			surv_input = 'Input_Zscore'
		elif len(input_data)==1:
			test_set = expr.loc[input_data[0]]

			symbol = pcta_id.loc[input_data[0]]['Symbol']
			test_set = test_set.rename(symbol)

			test_set_arr.append(test_set)
			#surv_input = input_data[0]
			surv_input = symbol

	current_task.update_state(state='PROGRESS', meta={'process_percent1': 40})

	dv_arr = []
	for i, test in enumerate(test_set_arr):
		test_set_up = test[test > test.median()]
		test_set_down = test[test < test.median()]
		clinical_up = clinical_set_arr[i].loc[test_set_up.index.tolist()]
		clinical_down = clinical_set_arr[i].loc[test_set_down.index.tolist()]
		clinical_arr = [clinical_up, clinical_down]

		dv_arr.append(clinical_arr)

	current_task.update_state(state='PROGRESS', meta={'process_percent1': 60})
	for dv in dv_arr:
		pt.survival_plot_and_cox(dv, label=[surv_input+'_up', surv_input+'_down'], filename=plot_path+test_name+str(filecount))
		file_list.append(template_plot_path+test_name+str(filecount)+".png")
		filecount+=1

	current_task.update_state(state='PROGRESS', meta={'process_percent1': 80})

	table_arr = []
	cph = CoxPHFitter()
	for i, test in enumerate(test_set_arr):
		cph_data = test
		cph_data.loc[cph_data>test_set.median()]=1
		cph_data.loc[cph_data<test_set.median()]=0
		cph_data = pd.concat([cph_data, clinical_set_arr[i]], axis=1)

		table_data = []
		cph.fit(cph_data[[surv_input,'bcrstatus','bcrmonth']].dropna(), 'bcrmonth', event_col='bcrstatus')
		table_data = table_split(cph)

		cph.fit(cph_data[['gleason','bcrstatus','bcrmonth']].dropna(), 'bcrmonth', event_col='bcrstatus')
		table_data = table_data+table_split(cph)

		cph.fit(cph_data[[surv_input,'gleason','bcrstatus','bcrmonth']].dropna(), 'bcrmonth', event_col='bcrstatus')
		table_data = table_data+table_split(cph)

		#table_data = [['Variable','Coefficient','Hazard Ratio','P-value','C-index']]+table_data

		table_arr.append(table_data)
		rw = open(plot_path+surv_dataset[i]+'.cox.tsv',"w")
		rw.write("Variable\tCoefficient\tHR(95%Conf)\tP-value\tC-index\n")

		for x in table_data:
			rw.write('\t'.join(x)+'\n')
		rw.close()

	os.system("cp -rf %s/* %s"%(plot_path,nginx_plot_path))
	current_task.update_state(state='PROGRESS', meta={'process_percent1': 100})

	return random.random()

@shared_task
def correlation_c(input_data1,input_data2, dataset_input_data, group_info, group_samples, plot_color, session_key, name1='Input 1', name2='Input 2'):

	current_task.update_state(state='PROGRESS', meta={'process_percent2': 0})
	####ALL SAMPLES INSERTION####
	group_info = ['ALL']+group_info
	all_samples = [sample for samples in group_samples for sample in samples]
	group_samples.insert(0, all_samples)
	####ALL SAMPLES INSERTION####

	def input_data_process(df,input_d):
		if len(input_d)>1:
			df_data = [gene_set_zscore_single_thr(df[samples], input_d, sample_status="multiple") for samples in group_samples]

		elif len(input_d)==1:
			df_data = [df[samples].loc[input_d[0]].values.tolist() for samples in group_samples]
		return df_data

	current_task.update_state(state='PROGRESS', meta={'process_percent2': 20})
	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+session_key+"/"
	nginx_plot_path= "/home/ubuntu/django_proj/pcta_updated/main/staticimages/"+session_key+"/"

	template_plot_path = "images/"+session_key+"/"
	test_name = "correlation_test"

	if not os.path.exists(plot_path):
		os.mkdir(plot_path)
		os.mkdir(nginx_plot_path)
	else:
		files = glob.glob(plot_path+"*")
		for f in files:
			os.remove(f)
		files = glob.glob(nginx_plot_path+"*")
		for f in files:
			os.remove(f)

	filecount = 0
	file_list = []

	#####Scatter
	pt=MY_PLOT()

	if dataset_input_data=='PCTA':
		df = all_expr_df
		df = df.drop('Symbol',axis=1)
	else:
		df = tcga_expr_df

	df_data1 = input_data_process(df, input_data1)
	current_task.update_state(state='PROGRESS', meta={'process_percent2': 40})
	df_data2 = input_data_process(df, input_data2)
	current_task.update_state(state='PROGRESS', meta={'process_percent2': 60})
	group_df = []

	for a in range(len(group_info)):
		d = {name1:df_data1[a], name2:df_data2[a]}
		group_df.append(pd.DataFrame(data=d))
		current_task.update_state(state='PROGRESS', meta={'process_percent2': 80})

	pt.shared_scatter(group_df, x_dat=name1,y_dat=name2, title=group_info, color_arr=plot_color, filename=plot_path+test_name+str(filecount))
	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1

	os.system("cp -rf %s/* %s"%(plot_path,nginx_plot_path))
	current_task.update_state(state='PROGRESS', meta={'process_percent2': 100})

	return random.random()

@shared_task
def set_c(input_data, group_info, group_samples, session_key):

	def fet_f(a1,b1,total_gene_assum):
		a1_inter_b1 = list(set(a1).intersection(b1))
		a1_unique_fromb1 = list(set(a1)-set(a1_inter_b1))
		b1_unique_froma1 = list(set(b1)-set(a1_inter_b1))

		oddsratio, pvalue = stats.fisher_exact([[len(a1_inter_b1), len(b1_unique_froma1)], [len(a1_unique_fromb1), total_gene_assum-(len(a1_inter_b1)+len(b1_unique_froma1)+len(a1_unique_fromb1))]])
		return len(a1),len(a1_inter_b1),pvalue

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+session_key+"/"
	template_plot_path = "images/"+session_key+"/"
	test_name = "USER_SET"

	if not os.path.exists(plot_path):
		os.mkdir(plot_path)
	else:
		files = glob.glob(plot_path+"*")
		for f in files:
			os.remove(f)

	filecount = 0
	file_list = []

	pt=MY_PLOT()
	df_all = all_expr_df

	#############GSEA#############
	gmt_temp = 'USER_SET\tNA\t'+'\t'.join(input_data)
	#fixed_path_gmt = 'user_data/'+userID+'/user.gmt'
	fixed_path_gmt = plot_path+'user.gmt'
	rw = open(fixed_path_gmt,'w')
	rw.write(gmt_temp)
	rw.close()

	sample_list = group_samples
	class_vector = [[group_info[i]]*len(item) for i,item in enumerate(sample_list)]
	class_vector = [y for x in class_vector for y in x]
	class_vector = map(str,class_vector)

	current_task.update_state(state='PROGRESS', meta={'process_percent3': 20})

	df_user_s = [df_all[s] for s in sample_list]
	df_user_s = pd.concat(df_user_s,axis=1)

	df_user_s.columns = range(len(df_user_s.columns.tolist()))
	df_user_s = df_user_s.reset_index()

	gseapy.gsea(data=df_user_s, gene_sets=fixed_path_gmt, cls=class_vector, outdir=plot_path, min_size=2, max_size=1000, weighted_score_type=1, permutation_type = 'gene_set', method='signal_to_noise', ascending=False,figsize=(6.5,6), format='png')

	file_list.append(template_plot_path+"USER_SET.gsea.png")
	filecount += 1
	#############GSEA#############
	current_task.update_state(state='PROGRESS', meta={'process_percent3': 40})
	#############MRA#############
	mra_set_t = mra_set.T
	mra_list = list(set(mra_set.index.tolist()))
	mra_targets = [mra_set_t[x].values.tolist() for x in mra_list]
	mra_targets = [map(str,x[0]) if type(x[0])==list else [str(x[0])] for x in mra_targets]

	total_genes = len(list(set(mra_set.values.flatten())))

	pvals = [fet_f(mra_targets[a], input_data, total_genes) for a in range(len(mra_list))]
	pvals_list = [[mra_list[i],item[0],item[1],item[2]]for i,item in enumerate(pvals) if item[2] < 0.01 and item[0]>10]
	table_arr = pvals_list #####Table data

	current_task.update_state(state='PROGRESS', meta={'process_percent3': 60})

	rw = open(plot_path+'mra_candidates.tsv',"w")
	rw.write("Gene(EntrezID)\tTF_targets\tMapped_genes\tP-value\n")
	for x in table_arr:
		x = [str(y) for y in x]
		rw.write('\t'.join(x)+'\n')
	rw.close()

	network_data = pd.DataFrame(data=pvals_list, columns=['TF', 'targets','mapped','pval'])
	network_data = network_data.set_index('TF')
	network_data['prob_mapped'] = network_data['mapped']/network_data['targets']
	network_data = network_data.sort_values('prob_mapped',ascending=False)
	network_data = network_data.loc[network_data.index.tolist()[:10]]

	current_task.update_state(state='PROGRESS', meta={'process_percent': 80})

	selected_network_expr = df_all[group_samples[0]].loc[network_data.index.tolist()[:10]]

	pt.network_plot(network_data, selected_network_expr, tit=group_info[0],filename=plot_path+test_name+str(filecount))
	file_list.append(template_plot_path+test_name+str(filecount)+".png")

	current_task.update_state(state='PROGRESS', meta={'process_percent3': 100})

	#############MRA#############
	return random.random()
