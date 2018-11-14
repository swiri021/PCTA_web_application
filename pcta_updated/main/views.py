# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import unicode_literals
from django.shortcuts import render

from django.contrib.auth.decorators import login_required
from django.shortcuts import render, get_object_or_404, render_to_response, HttpResponseRedirect, HttpResponse, redirect, Http404
from django.contrib.auth import authenticate, login, logout
from django.db.models import Count

from django.conf import settings
from django.urls import reverse
from .models import *
from .forms import *
import StringIO
import zipfile
import glob

from PLOTS import MY_PLOT
import numpy as np
import pandas as pd
import gseapy
from scipy import stats
from celery.result import AsyncResult

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from multiprocessing import Process, Queue
from pcta_updated import all_expr_df, pcta_id, mra_set
from wand.image import Image
from wand.color import Color

import json
import os
import ast
import shutil
import os.path


from .tasks import association_c, set_c, survival_c, correlation_c
from celery.result import AsyncResult
import json

def poll_state(request):

	""" A view to report the progress to the user """
	data = 'Fail'
	if request.is_ajax():
		if 'task_id' in request.POST.keys() and request.POST['task_id']:
			task_id = request.POST['task_id']
			task = AsyncResult(task_id)
			data = task.result or task.state
		else:
			data = 'No task_id in the request;'
	else:
		data = 'This is not an ajax request'

	json_data = json.dumps(data)
	print json_data

	return HttpResponse(json_data, content_type='application/json')

################PAGE VIEW################
def initial_change(request):
	pathway_input = Pathway_option(request.POST)
	if request.method == 'POST':

		if pathway_input.is_valid():
			if request.POST['input_type']=='assoc':
				pathway_list = request.POST.getlist('pathway_select')
				for x in pathway_list:
					temp_genes = PathwayGenes.objects.filter(pathway_name=x).values_list('entrez_id',flat=True)
					temp_genes = [str(y) for y in temp_genes]
				temp_genes = "\n".join(temp_genes)
				request.session['pathway_assoc'] = temp_genes
				return redirect(reverse('main:main_template'))

			elif request.POST['input_type']=='corr1':
				pathway_list = request.POST.getlist('pathway_select')
				for x in pathway_list:
					temp_genes = PathwayGenes.objects.filter(pathway_name=x).values_list('entrez_id',flat=True)
					temp_genes = [str(y) for y in temp_genes]
				temp_genes = "\n".join(temp_genes)
				request.session['pathway_corr1'] = temp_genes
				request.session['pathway_name1'] = pathway_list[0]

				return redirect(reverse('main:main_template'))

			elif request.POST['input_type']=='corr2':
				pathway_list = request.POST.getlist('pathway_select')
				for x in pathway_list:
					temp_genes = PathwayGenes.objects.filter(pathway_name=x).values_list('entrez_id',flat=True)
					temp_genes = [str(y) for y in temp_genes]
				temp_genes = "\n".join(temp_genes)
				request.session['pathway_corr2'] = temp_genes
				request.session['pathway_name2'] = pathway_list[0]
				return redirect(reverse('main:main_template'))

			elif request.POST['input_type']=='set':
				pathway_list = request.POST.getlist('pathway_select')
				for x in pathway_list:
					temp_genes = PathwayGenes.objects.filter(pathway_name=x).values_list('entrez_id',flat=True)
					temp_genes = [str(y) for y in temp_genes]
				temp_genes = "\n".join(temp_genes)
				request.session['pathway_set'] = temp_genes
				request.session['pathway_set_name'] = pathway_list[0]
				return redirect(reverse('main:main_template'))

		else:
			msg = ''
			for key, value in pathway_input.errors.items():
				value = json.dumps(value)
				value = ast.literal_eval(str(value))
				msg = msg+" "+value[0]
			return render(request,'main/error.html',{'msg':msg})

def main_template(request):
	input_key = [key for key in request.session.keys() if key.find('pathway_')>-1]

	assoc_input = AssocInputText(request.POST)
	corr_input = CorrInputText(request.POST)
	set_input = SetInputText(request.POST)
	pathway_input = Pathway_option(request.POST)
	corr_name_input = CorrInputName(request.POST)
	set_name_input = SetInputName(request.POST)

	if request.method == 'POST':

		assoc_dat = request.POST['gene_input']
		corr_dat1 = request.POST['corr_input1']
		corr_dat2 = request.POST['corr_input2']
		set_dat = request.POST['set_input']
		request.session['analysis_type'] = request.POST['analysis_type']

		#print request.POST['pathway_select']

		if assoc_dat=='' and corr_dat1=='' and corr_dat2=='' and set_dat=='': ######no input check
			msg = "No input, try again"
			return render(request,'main/error.html',{'msg':msg})

		if request.POST['analysis_type']=='Association analysis start':

			input_g = input_process(assoc_dat)

			if len(input_g)>1000:
				msg = "Too many genes, Try it less than 1000 Genes"
				return render(request,'main/error.html',{'msg':msg})

			if input_scanning(input_g)==False:
				msg = "Gene does not exist in DB, Try another gene"
				return render(request,'main/error.html',{'msg':msg})

			request.session['gene_input']=input_g

		elif request.POST['analysis_type']=='Correlation analysis start':

			input_g1 = input_process(corr_dat1)
			input_g2 = input_process(corr_dat2)

			if len(input_g1)>1000 and len(input_g2)>1000:
				msg = "Too many genes, Try it less than 1000 Genes at both boxes"
				return render(request,'main/error.html',{'msg':msg})

			elif len(input_g1)>1000:
				msg = "Too many genes, Try it less than 1000 Genes at first box"
				return render(request,'main/error.html',{'msg':msg})

			elif len(input_g2)>1000:
				msg = "Too many genes, Try it less than 1000 Genes at second box"
				return render(request,'main/error.html',{'msg':msg})


			if len(input_g1)==0:
				msg = "At first box, Please type in gene or genelist"
				return render(request,'main/error.html',{'msg':msg})

			elif len(input_g2)==0:
				msg = "At second box, Please type in gene or genelist"
				return render(request,'main/error.html',{'msg':msg})


			if input_scanning(input_g1)==False and input_scanning(input_g2)==False:
				msg = "Gene in both boxes does not exist in DB, Try another gene"
				return render(request,'main/error.html',{'msg':msg})

			elif input_scanning(input_g1)==False:
				msg = "Gene in first box does not exist in DB, Try another gene"
				return render(request,'main/error.html',{'msg':msg})

			elif input_scanning(input_g2)==False:
				msg = "Gene in second box does not exist in DB, Try another gene"
				return render(request,'main/error.html',{'msg':msg})

			request.session['corr_input1']=input_g1
			request.session['corr_input2']=input_g2

			if request.POST['name_input1']==request.POST['name_input2']:
				request.session['corr_name_input1'] = request.POST['name_input1']
				request.session['corr_name_input2'] = request.POST['name_input2']+"_duplicated"
			else:
				request.session['corr_name_input1'] = request.POST['name_input1']
				request.session['corr_name_input2'] = request.POST['name_input2']

		elif request.POST['analysis_type']=='Set analysis start':
			input_g = input_process(set_dat)

			if len(input_g)>1000:
				msg = "Too many genes, Try it less than 1000 Genes"
				return render(request,'main/error.html',{'msg':msg})

			if len(input_g)<10:
				msg = "Input should be more than 10 genes"
				return render(request,'main/error.html',{'msg':msg})

			if input_scanning(input_g)==False:
				msg = "Gene does not exist, Try another gene"
				return render(request,'main/error.html',{'msg':msg})

			request.session['set_name_input'] = request.POST['set_name_input']
			request.session['set_input']=input_process(set_dat)

		return redirect(reverse('main:option'))

	else:

		if 'pathway_assoc' in input_key:
			assoc_input = AssocInputText(initial={'gene_input':request.session['pathway_assoc']})
		else:
			assoc_input = AssocInputText()

		if 'pathway_corr1' in input_key and 'pathway_corr2' in input_key:
			corr_name_input = CorrInputName(initial={'name_input1':request.session['pathway_name1'],'name_input2':request.session['pathway_name2']})
			corr_input = CorrInputText(initial={'corr_input1':request.session['pathway_corr1'], 'corr_input2':request.session['pathway_corr2']})
		elif ('pathway_corr1' in input_key) and ('pathway_corr2' not in input_key):
			corr_name_input = CorrInputName(initial={'name_input1':request.session['pathway_name1']})
			corr_input = CorrInputText(initial={'corr_input1':request.session['pathway_corr1']})
		elif ('pathway_corr2' in input_key) and ('pathway_corr1' not in input_key):
			corr_name_input = CorrInputName(initial={'name_input2':request.session['pathway_name2']})
			corr_input = CorrInputText(initial={'corr_input2':request.session['pathway_corr2']})
		else:
			corr_name_input = CorrInputName()
			corr_input = CorrInputText()

		if 'pathway_set' in input_key:
			set_name_input = SetInputName(initial={'set_name_input':request.session['pathway_set_name']})
			set_input = SetInputText(initial={'set_input':request.session['pathway_set']})
		else:
			set_input = SetInputText()
			set_name_input = SetInputName()

		pathway_input = Pathway_option()

	return render(request, 'main/input_page.html', {'assoc_input':assoc_input, 'corr_input':corr_input, 'set_input':set_input, 'pathway_input':pathway_input, 'corr_name_input':corr_name_input, 'set_name_input':set_name_input})

def option(request):
	input_key = [key for key in request.session.keys() if key.find('pathway_')>-1]
	for key in input_key:
		del request.session[key]

	analysis_type = request.session['analysis_type']
	#print request.session.session_key

	if analysis_type=='Association analysis start':
		optionform = AssocOptionForm(request.POST)
		if request.method=='POST':
			request.session['sample_option'] = request.POST['sample_option']
			return redirect(reverse('main:calculation'))
		else:
			optionform = AssocOptionForm()

	elif analysis_type=='Correlation analysis start':
		optionform = CorrOptionForm(request.POST)
		if request.method=='POST':
			request.session['sample_option'] = request.POST['sample_option']
			return redirect(reverse('main:calculation'))
		else:
			optionform = CorrOptionForm()

	elif analysis_type=='Set analysis start':
		optionform = SetOptionForm(request.POST)
		if request.method=='POST':
			request.session['sample_option'] = request.POST['sample_option']
			return redirect(reverse('main:calculation'))
		else:
			optionform = SetOptionForm()

	return render(request, 'main/option.html', {'optionform':optionform})

def about(request):
	return render(request, 'main/about.html', {})

def manual(request):
	return render(request, 'main/manual.html', {})

def download(request):
	return render(request, 'main/download.html', {})

def qna(request):
	return render(request, 'main/qna.html',{})

def plot_color(input_option):
	if input_option=='pcs':
		return ['#C95E5E','#57C757','#5377C6']
	elif input_option=='pam50':
		return ['#2A7E70','#BB733E','#BF404A']
	elif input_option=='n_category':
		return ['#449933','#040301','#BB733E','#BF404A']

def call_calculation(request):
	analysis_type = request.session['analysis_type']

	if analysis_type=='Association analysis start':
		input_data = request.session['gene_input']
		input_option = request.session['sample_option']
		if input_option!='bcr':

			if 'job' in request.GET:
				job_id = request.GET['job']
				job = AsyncResult(job_id)
				data = job.result or job.state
				context = {
					'data':data,
					'task_id':job_id,
					'pass_to':1
				}
				return render(request,"main/loading.html",context)

			else:
				plot_color_arr = plot_color(input_option)
				group_samples, group_info = grouping_samples(input_option)
				print plot_color_arr
				job = association_c.delay(input_data, group_info, group_samples, plot_color_arr, request.session.session_key)
				return HttpResponseRedirect(reverse('main:calculation') + '?job=' + job.id)

		else:

			if 'job' in request.GET:
				job_id = request.GET['job']
				job = AsyncResult(job_id)
				data = job.result or job.state
				context = {
					'data':data,
					'task_id':job_id,
					'pass_to':2
				}
				return render(request,"main/loading1.html",context)

			else:
				job = survival_c.delay(input_data, request.session.session_key)
				return HttpResponseRedirect(reverse('main:calculation') + '?job=' + job.id)
				#file_list, table_list = survival_calculator(input_data, request.session.session_key)
				#request.session['bcr_result_list'] = file_list
				#request.session['bcr_table_list'] = table_list
				#return redirect(reverse('main:bcr_result'))

	elif analysis_type=='Correlation analysis start':
		#return redirect(reverse('main:correlation_result'))
		#request.session['corr_result_list'] = file_list
		if 'job' in request.GET:
			job_id = request.GET['job']
			job = AsyncResult(job_id)
			data = job.result or job.state
			context = {
				'data':data,
				'task_id':job_id,
				'pass_to':3
			}
			return render(request,"main/loading2.html",context)

		else:
			input_data1 = request.session['corr_input1']
			input_data2 = request.session['corr_input2']
			input_option = request.session['sample_option']
			group_samples, group_info = grouping_samples(input_option)
			plot_color_arr = plot_color(input_option)

			job = correlation_c.delay(input_data1,input_data2, group_info, group_samples, plot_color_arr, request.session.session_key, name1=request.session['corr_name_input1'], name2=request.session['corr_name_input2'])
			return HttpResponseRedirect(reverse('main:calculation') + '?job=' + job.id)

	elif analysis_type=='Set analysis start':

		input_data = request.session['set_input']
		input_option = request.session['sample_option']

		selected_category = input_option.split(":")[0]
		selected_data = input_option.split(":")[1]

		group_samples, group_info = grouping_samples(selected_category)
		group_data_ext = group_samples.pop(group_info.index(selected_data))
		other_group_ext = [g for group in group_samples for g in group]

		group_samples = [group_data_ext, other_group_ext]
		group_info = [selected_data,'Other']

		file_list, table_list, gsea_mapping_rate = set_calculator(input_data, group_info, group_samples, request.session.session_key, set_name=request.session['set_name_input'])
		request.session['set_result_list'] = file_list
		request.session['set_table_list'] = table_list
		request.session['gsea_mapping_rate'] = gsea_mapping_rate
		return redirect(reverse('main:set_result'))

		"""
		if 'job' in request.GET:
			job_id = request.GET['job']

			job = AsyncResult(job_id)
			data = job.result or job.state
			print data
			context = {
				'data':data,
				'task_id':job_id,
				'pass_to':4
			}
			return render(request,"main/loading3.html",context)

		else:
			input_data = request.session['set_input']
			input_option = request.session['sample_option']

			selected_category = input_option.split(":")[0]
			selected_data = input_option.split(":")[1]

			group_samples, group_info = grouping_samples(selected_category)
			group_data_ext = group_samples.pop(group_info.index(selected_data))
			other_group_ext = [g for group in group_samples for g in group]

			group_samples = [group_data_ext, other_group_ext]
			group_info = [selected_data,'Other']

			job = set_c.delay(input_data, group_info, group_samples, request.session.session_key)
			return HttpResponseRedirect(reverse('main:calculation') + '?job=' + job.id)
			"""

def association_result(request):
	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+request.session.session_key+"/"
	template_plot_path = "images/"+request.session.session_key+"/"
	test_name = "association_test"

	result = [template_plot_path+test_name+"0.png",template_plot_path+test_name+"1.png",template_plot_path+test_name+"2.png"]

	table_result = []
	surv_dataset = ['onewayANOVA_result', 'ranksum_result']
	f1=open(plot_path+surv_dataset[0]+'.tsv')
	dat1 = f1.read().strip().split("\n")
	dat1 = dat1[1:]
	dat1 = [d.split('\t') for d in dat1]
	table_result.append(dat1)

	f1=open(plot_path+surv_dataset[1]+'.tsv')
	dat1 = f1.read().strip().split("\n")
	dat1 = dat1[1:]
	dat1 = [d.split('\t') for d in dat1]
	table_result.append(dat1)

	return render(request, 'main/association_result.html', {'result':result, 'table_result1':table_result[0], 'table_result2':table_result[1]})

def bcr_result(request):
	#result = request.session['bcr_result_list']
	#table_result = request.session['bcr_table_list']

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+request.session.session_key+"/"
	template_plot_path = "images/"+request.session.session_key+"/"
	test_name = "bcr_analysis"

	result = [template_plot_path+test_name+"0.png",template_plot_path+test_name+"1.png"]

	table_result = []
	surv_dataset = ['GSE40272', 'GSE70769']

	f1=open(plot_path+surv_dataset[0]+'.cox.tsv')
	dat1 = f1.read().strip().split("\n")
	dat1 = dat1[1:]
	dat1 = [d.split('\t') for d in dat1]
	table_result.append(dat1)

	f1=open(plot_path+surv_dataset[1]+'.cox.tsv')
	dat1 = f1.read().strip().split("\n")
	dat1 = dat1[1:]
	dat1 = [d.split('\t') for d in dat1]
	table_result.append(dat1)

	#print table_result[0]

	return render(request, 'main/bcr_result.html', {'result':result, 'table_result1':table_result[0], 'table_result2':table_result[1]})

def correlation_result(request):
	#result = request.session['corr_result_list']

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+request.session.session_key+"/"
	template_plot_path = "images/"+request.session.session_key+"/"
	test_name = "correlation_test"
	result = [template_plot_path+test_name+"0.png",template_plot_path+test_name+"1.png",template_plot_path+test_name+"2.png"]

	return render(request, 'main/correlation_result.html', {'result':result})

def set_result(request):
	result = request.session['set_result_list']
	table_result = request.session['set_table_list']
	add_info = request.session['gsea_mapping_rate']
	return render(request, 'main/set_result.html', {'result':result, 'table_result':table_result, 'add_info':add_info})
################PAGE VIEW################

################PAGE FUNCTION################

def zip_f(src, dst, fmat):

	memory_io = StringIO.StringIO()

	#zf = zipfile.ZipFile("%s.zip" % (dst), "w", zipfile.ZIP_DEFLATED)
	zf = zipfile.ZipFile(memory_io, "w")
	abs_src = os.path.abspath(src)

	for dirname, subdirs, files in os.walk(src):
		for filename in files:
			absname = os.path.abspath(os.path.join(dirname, filename))
			arcname = absname[len(abs_src) + 1:]
			#print 'zipping %s as %s' % (os.path.join(dirname, filename),arcname)

			if arcname.find(fmat)!=-1:
				print 'zipping %s as %s' % (os.path.join(dirname, filename),arcname)
				zf.write(absname, arcname)

	zf.close()
	return memory_io

def image_download(request):

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+request.session.session_key+"/"
	if request.method == 'POST':
		filename = request.POST['filetype']
		print filename
		if os.path.exists(plot_path):
			pdf_list = glob.glob(plot_path+"*.pdf")
			memory_io = zip_f(plot_path, plot_path+'file_dowload', filename)

		#zf = open(plot_path+'figures.zip',"rb")
		response = HttpResponse(memory_io.getvalue(),content_type='Content-Type: application/x-zip-compressed')
		response['Content-Disposition'] = 'attachment; filename=%s' % os.path.basename(plot_path+'file_dowload.zip')

		return response

def file_download(request):
	if request.method == 'POST':
		filename = request.POST['file_list']
		file_path = os.path.join(settings.MEDIA_ROOT, filename)

		print file_path
		if os.path.exists(file_path):
			with open(file_path, 'rb') as fh:
				response = HttpResponse(fh.read(), content_type='text/csv')
				response['Content-Disposition'] = 'attachment; filename=' + os.path.basename(file_path)
				#response['Content-Disposition'] = 'attachment; filename="file1.csv"'
				return response

			raise Http404

################PAGE FUNCTION################


###############Calculation Section###############

def id_searching(symbol):
	entrez = pcta_id.loc[pcta_id['Symbol']==symbol].index.tolist()
	if entrez==[]:
		"""
		syno = pcta_id['Synonyms'].values.tolist()
		count = 0
		while count < len(syno):
			matched = [syno[count] for a in str(syno[count]).split("|") if a==symbol]
			if matched!=[]:
				entrez = pcta_id.iloc[count].index.tolist()
				count = len(syno)
			count+=1

		if entrez==[]:
			return 'NaN'
		else:
			entrez[0]
		"""
		return 'NaN'

	else:
		return entrez[0]

def input_scanning(input_data):

	gene_list = pcta_id.index.tolist()
	check_list = [True for inp in input_data if inp in gene_list]
	true_count = check_list.count(True)

	if true_count==0:
		return False
	else:
		return True

def input_process(input_data):
	input_data = input_data.strip()
	input_data = input_data.replace("\r\n","\n")
	input_data = input_data.split("\n")
	input_data = [ str(i).encode('utf8') for i in input_data if i!='' ]
	input_data = list(set(input_data))
	input_data = [i.upper() for i in input_data]
	input_data = [ str(id_searching(i)) if str(i).isdigit()==False else i for i in input_data ]

	return input_data

def grouping_selection(grade, user_selection):
	if user_selection=='n_category':
		return MajorGrade.objects.filter(n_category=grade)
	elif user_selection=='disease_status':
		return MajorGrade.objects.filter(disease_status=grade)
	elif user_selection=='pcs':
		return MajorGrade.objects.filter(pcs=grade)
	elif user_selection=='pam50':
		return MajorGrade.objects.filter(pam50=grade)
	#elif user_selection=='all':
	#	return MajorGrade.objects.all()
	else:
		return MajorGrade.objects.filter(n_category=grade)#####default

def grouping_samples(input_option):
	m_grad = MajorGrade.objects.values(input_option).exclude(**{input_option : 'NULL'}).exclude(**{input_option : 'Normal'}).exclude(**{input_option : 'Benign'}).distinct().values_list(input_option, flat=True)
	m_grad = sorted(m_grad)
	group_samples = [grouping_selection(g, input_option) for g in m_grad]
	sample_list = [map(str,list(grade_orm.values('sample_id').values_list('sample_id', flat=True))) for grade_orm in group_samples]

	return sample_list, m_grad
"""
def gene_set_zscore( arr1, gene_set=[] ,sample_status="multiple"):

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
		zscore = [cal_z(arr1[x], inter) for x in arr1.columns.tolist()]

	return zscore
"""
def gene_set_zscore(df_sample, gene_set=[] ,sample_status="multiple"):

	def array_divide(piv, arr):
			res = []

			if len(arr)%piv==0:
				circ = int(round(len(arr)/piv))
			else:
				circ = int(round(len(arr)/piv) + 1)
			numb = 0

			for x in range(circ):
				try:
					numb = x*piv
					a = arr[numb:numb+piv]
					#numb = numb+100
				except IndexError:
					a = arr[numb:]

				res.append(a)

			return res

	def cal_z(temp_df, inter):

			selected = temp_df.loc[inter].tolist()
			selected_adj = [float(x) for x in selected]
			total = temp_df.tolist()
			total_adj = [float(x) for x in total]

			diff_mean = np.nanmean(selected_adj) - np.nanmean(total_adj)
			result = diff_mean*np.sqrt(len(selected_adj))/np.nanstd(total_adj,ddof=1)

			return result

	def multi_proc(df_head,order):
			result = [cal_z(df_sample[x], inter) for x in df_head]
			queue.put((result,order))

	ft_dat_index = gene_set
	arr1_index = df_sample.index.tolist()

	inter = list(set(arr1_index).intersection(ft_dat_index))

	if sample_status=="single":
		zscore = [cal_z(df_sample, inter)]

	elif sample_status=="multiple":

		queue = Queue()
		new_multi_head = array_divide(50,df_sample.columns.tolist())

		processes = [Process(target=multi_proc, args=(item,i)) for i,item in enumerate(new_multi_head)]
		for p in processes:
			p.start()

		for p in processes:
			p.join()
		results = [queue.get() for p in processes]
		results = sorted(results, key=lambda tup: tup[1])
		results = [x[0] for x in results]
		zscore = [y for x in results for y in x]

		#zscore = [cal_z(df_sample[x], inter) for x in df_sample.columns.tolist()]

	return zscore

def association_calculator(input_data, group_info, group_samples, session_key):

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+session_key+"/"
	template_plot_path = "images/"+session_key+"/"
	test_name = "association_test"

	if not os.path.exists(plot_path):
		os.mkdir(plot_path)
	else:
		files = glob.glob(plot_path+"*")
		print files
		for f in files:
			os.remove(f)

	filecount = 0
	file_list = []

	#####Lollipop
	pt=MY_PLOT()
	df = all_expr_df
	df = df.drop('Symbol',axis=1)

	if len(input_data)>1:
		df_data = [gene_set_zscore(df[samples], input_data, sample_status="multiple") for samples in group_samples]
		ylab = 'Z score'
	elif len(input_data)==1:
		df_data = [df[samples].loc[input_data[0]].values.tolist() for samples in group_samples]
		ylab = 'Expression'

	pt.lollipop([sorted(d) for d in df_data], ylab=ylab, label_data=group_info, filename=plot_path+test_name+str(filecount))
	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1
	#####Lollipop

	#####Violin
	df_data_arr = [ pd.DataFrame(data=d) for d in df_data]
	df_data_arr = pd.concat(df_data_arr, axis=1)
	df_data_arr.columns = group_info

	pt.violin_plt(df_data_arr,ylab=ylab,tit='',filename=plot_path+test_name+str(filecount))
	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1
	#####Violin

	#####Histogram
	pt.histogram_group(df_data_arr,legend=True,filename=plot_path+test_name+str(filecount))
	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1
	#####Histogram

	return file_list

def survival_calculator(input_data, session_key):

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+session_key+"/"

	template_plot_path = "images/"+session_key+"/"
	test_name = "bcr_analysis"

	if not os.path.exists(plot_path):
		os.mkdir(plot_path)
	else:
		files = glob.glob(plot_path+"*")
		for f in files:
			os.remove(f)

	filecount = 0
	file_list = []

	pt=MY_PLOT()

	surv_dataset = ['GSE40272', 'GSE70769']

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

				t_arr.append([testing+"(Multivariate)", coef, haz_c+"("+haz_c_int+")", pval, concor])

			return t_arr
		else:
			for x in cph_d_index:
				testing = x
				coef = format(cph_d['coef'].loc[x],".2f")
				haz_c = format(cph_d['exp(coef)'].loc[x],".2f")
				haz_c_int = str(format(np.exp(cph_d['lower 0.95'].loc[x]),".2f"))+"-"+str(format(np.exp(cph_d['upper 0.95'].loc[x]),".2f"))
				pval = format(cph_d['p'].loc[x],".3f")
				concor = format(concordance_index(cph.durations, -cph.predict_partial_hazard(cph.data).values.ravel(), cph.event_observed),".2f")

			return [[testing, coef, haz_c+"("+haz_c_int+")", pval, concor]]

	expr_set_arr = []
	clinical_set_arr = []

	for surv in  surv_dataset:
		expr_set = pd.read_csv('user_data/clinical_data/'+surv+'.csv',index_col=0)
		expr_set.index = expr_set.index.astype(int)
		expr_set.index = expr_set.index.astype(str)
		clinical_set = pd.read_excel('user_data/clinical_data/'+surv+'_clinical.xlsx',index_col=0)

		expr_set_arr.append(expr_set)
		clinical_set_arr.append(clinical_set)

	test_set_arr = []
	for expr in expr_set_arr:
		if len(input_data)>1:
			df_data = gene_set_zscore(expr, input_data, sample_status="multiple")
			test_set = pd.Series(data=df_data, index=expr.columns.tolist(), name='Input_Zscore')
			test_set_arr.append(test_set)
			surv_input = 'Input_Zscore'
		elif len(input_data)==1:
			test_set = expr.loc[input_data[0]]
			test_set_arr.append(test_set)
			surv_input = input_data[0]

	dv_arr = []
	for i, test in enumerate(test_set_arr):
		test_set_up = test[test > test.median()]
		test_set_down = test[test < test.median()]
		clinical_up = clinical_set_arr[i].loc[test_set_up.index.tolist()]
		clinical_down = clinical_set_arr[i].loc[test_set_down.index.tolist()]
		clinical_arr = [clinical_up, clinical_down]

		dv_arr.append(clinical_arr)

	for dv in dv_arr:
		pt.survival_plot_and_cox(dv, label=[surv_input+'_up', surv_input+'_down'], filename=plot_path+test_name+str(filecount))
		file_list.append(template_plot_path+test_name+str(filecount)+".png")
		filecount+=1

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

	return file_list, table_arr

def correlation_calculator(input_data1,input_data2, group_info, group_samples, session_key):

	####ALL SAMPLES INSERTION####
	group_info = ['ALL']+group_info
	all_samples = [sample for samples in group_samples for sample in samples]
	group_samples.insert(0, all_samples)
	####ALL SAMPLES INSERTION####

	def input_data_process(df,input_d):
		if len(input_d)>1:
			df_data = [gene_set_zscore(df[samples], input_d, sample_status="multiple") for samples in group_samples]
		elif len(input_d)==1:
			df_data = [df[samples].loc[input_d[0]].values.tolist() for samples in group_samples]
		return df_data

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

	filecount = 0
	file_list = []

	#####Scatter
	pt=MY_PLOT()
	df = all_expr_df
	df = df.drop('Symbol',axis=1)

	df_data1 = input_data_process(df, input_data1)
	df_data2 = input_data_process(df, input_data2)
	group_df = []

	for a in range(len(group_info)):
		d = {'Input1':df_data1[a], 'Input2':df_data2[a]}
		group_df.append(pd.DataFrame(data=d))

	pt.shared_scatter(group_df, x_dat='Input1',y_dat='Input2', title=group_info, filename=plot_path+test_name+str(filecount))
	file_list.append(template_plot_path+test_name+str(filecount)+".png")
	filecount += 1

	os.system("cp %s/* %s"%(plot_path,nginx_plot_path))

	return file_list

def set_calculator(input_data, group_info, group_samples, session_key, set_name='USER_GENE_SET'):

	def fet_f(a1,b1,total_gene_assum):
		a1_inter_b1 = list(set(a1).intersection(b1))
		a1_unique_fromb1 = list(set(a1)-set(a1_inter_b1))
		b1_unique_froma1 = list(set(b1)-set(a1_inter_b1))

		oddsratio, pvalue = stats.fisher_exact([[len(a1_inter_b1), len(b1_unique_froma1)], [len(a1_unique_fromb1), total_gene_assum-(len(a1_inter_b1)+len(b1_unique_froma1)+len(a1_unique_fromb1))]])
		return len(a1),len(a1_inter_b1),pvalue

	plot_path= "/home/ubuntu/django_proj/pcta_updated/main/static/images/"+session_key+"/"
	nginx_plot_path= "/home/ubuntu/django_proj/pcta_updated/main/staticimages/"+session_key+"/"

	template_plot_path = "images/"+session_key+"/"
	test_name = "USER_SET"

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
	df_all = all_expr_df

	gsea_mapping_rate = len(list(set(input_data).intersection(df_all.index.tolist())))
	list_numb = len(input_data)
	gsea_mapping_rate = float(gsea_mapping_rate)/float(list_numb)*100
	gsea_mapping_rate = "%.2f" % gsea_mapping_rate

	#############GSEA#############
	gmt_temp = set_name+'\tNA\t'+'\t'.join(input_data)
	#fixed_path_gmt = 'user_data/'+userID+'/user.gmt'
	fixed_path_gmt = plot_path+'user.gmt'
	rw = open(fixed_path_gmt,'w')
	rw.write(gmt_temp)
	rw.close()

	sample_list = group_samples
	class_vector = [[group_info[i]]*len(item) for i,item in enumerate(sample_list)]
	class_vector = [y for x in class_vector for y in x]
	class_vector = map(str,class_vector)

	df_user_s = [df_all[s] for s in sample_list]
	df_user_s = pd.concat(df_user_s,axis=1)

	df_user_s.columns = range(len(df_user_s.columns.tolist()))
	df_user_s = df_user_s.reset_index()

	gsea_result = gseapy.gsea(data=df_user_s, gene_sets=fixed_path_gmt, cls=class_vector, outdir=plot_path, min_size=2, max_size=1000, weighted_score_type=1, permutation_type = 'phenotype', method='signal_to_noise', ascending=False,figsize=(6.5,6), format='png')
	#ledge_genes = gsea_result.results[gsea_result.results.keys()[0]]['ledge_genes'].split(";")## leading edge subset
	pval_es = gsea_result.results[gsea_result.results.keys()[0]]['pval'] ###GSEA p-value

	#with Image(filename=plot_path+set_name+".gsea.pdf", resolution=300) as img:
	#	with Image(width=img.width, height=img.height, background=Color("white")) as bg:
	#		bg.composite(img,0,0)
	#		bg.save(filename=plot_path+set_name+".gsea.png")

	if pval_es<=0.05:

		file_list.append(template_plot_path+set_name+".gsea.png")
		filecount += 1
	else:
		file_list.append(" ")
		gsea_mapping_rate= 'Not Applicable'
		filecount += 1
	#############GSEA#############

	fold_change = all_expr_df[sample_list[0]].median(axis=1) - all_expr_df[sample_list[1]].median(axis=1)


	if pval_es<=0.05:
		#############MRA#############
		mra_set_t = mra_set.T
		mra_list = list(set(mra_set.index.tolist()))
		mra_targets = [mra_set_t[x].values.tolist() for x in mra_list]
		mra_targets = [map(str,x[0]) if type(x[0])==list else [str(x[0])] for x in mra_targets]

		total_genes = len(list(set(mra_set.values.flatten())))

		pvals = [fet_f(mra_targets[a], input_data, total_genes) for a in range(len(mra_list))]
		#pvals = [fet_f(mra_targets[a], ledge_genes, total_genes) for a in range(len(mra_list))] ## leading edge subset
		pvals_list = [[mra_list[i],int(item[0]),int(item[1]),float("{0:.4f}".format(item[2])), float("{0:.4f}".format(fold_change.loc[mra_list[i]]))]for i,item in enumerate(pvals) if item[2] < 0.01 and fold_change.loc[mra_list[i]] >= 0.1 and item[0] > 10]
		#table_arr = pvals_list #####Table data

		index_change = lambda x,y : [y]+x[1:]
		#table_arr = [index_change(x,pcta_id.loc[x[0]]['Symbol']) for x in table_arr]

		network_data = pd.DataFrame(data=pvals_list, columns=['TF', 'targets','mapped','pval','fc'])
		network_data = network_data.set_index('TF')
		network_data['Symbol'] = pcta_id.loc[map(str, network_data.index.tolist())]['Symbol']

		network_data['prob_mapped'] = network_data['mapped']/network_data['targets']
		network_data = network_data.sort_values('prob_mapped',ascending=False)
		network_data = network_data.loc[network_data.index.tolist()[:10]]
		network_data = network_data.round(4)
		network_data['targets'].astype(int)
		network_data['mapped'].astype(int)

		selected_network_expr = df_all[group_samples[0]].loc[network_data.index.tolist()[:10]]
		selected_network_expr['Symbol'] = pcta_id.loc[map(str, selected_network_expr.index.tolist())]['Symbol']

		network_data = network_data.set_index('Symbol')
		selected_network_expr = selected_network_expr.set_index('Symbol')

		pt.network_plot(network_data, selected_network_expr, tit=group_info[0],filename=plot_path+test_name+str(filecount))
		file_list.append(template_plot_path+test_name+str(filecount)+".png")

		network_data = network_data[['targets', 'mapped', 'prob_mapped', 'fc', 'pval']]

		int_ = lambda x : [int(x[0]), int(x[1])] + x[2:]
		table_arr = [[i]+int_(network_data.loc[i].tolist()) for i in network_data.index.tolist()]

		network_data.to_csv(plot_path+'mra_candidates.csv')

		os.system("cp -rf %s/* %s"%(plot_path,nginx_plot_path))

		#############MRA#############

	else:
		table_arr=[]
		file_list.append(" ")

	return file_list, table_arr, gsea_mapping_rate

###############Calculation Section###############
