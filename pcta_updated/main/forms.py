from django import forms
from .models import *
from django.contrib.auth.forms import UserCreationForm, AuthenticationForm
from django.contrib.auth.models import User


class AssocInputText(forms.Form):
	gene_input = forms.CharField(label='',required=False, widget=forms.Textarea(attrs={'placeholder': 'Copy and Paste Gene or Gene list(Entrez ID or Official Gene symbol)'}))
	def clean(self):
		value1 = self.cleaned_data['gene_input']
		value1_length = value1.split("\n")

		if len(value1_length) > 1000 and value1!='':
			raise forms.ValidationError("Input should be less than 1000")
		return value1

class CorrInputText(forms.Form):
	corr_input1 = forms.CharField(label='',required=False, widget=forms.Textarea(attrs={'placeholder': 'First Gene or Gene List (Entrez ID or Official Gene symbol)'}))
	corr_input2 = forms.CharField(label='',required=False, widget=forms.Textarea(attrs={'placeholder': 'Second Gene or Gene List (Entrez ID or Official Gene symbol)'}))
	def clean(self):
		value1 = self.cleaned_data['corr_input1']
		value1_length = value1.split("\n")
		value2 = self.cleaned_data['corr_input2']
		value2_length = value2.split("\n")

		if len(value1_length) > 1000 and value1!='':
			raise forms.ValidationError("Input should be less than 1000")
		if len(value2_length) > 1000 and value2!='':
			raise forms.ValidationError("Input should be less than 1000")

		if value1!='' and value2=='':
			raise forms.ValidationError("Please input second list for correlation")
		if value2!='' and value1=='':
			raise forms.ValidationError("Please input first list for correlation")
		return value1

class SetInputText(forms.Form):

	set_input = forms.CharField(label='',required=False, widget=forms.Textarea(attrs={'placeholder': 'Copy and Paste Gene list only (Entrez ID or Official Gene symbol)'}))
	def clean(self):
		value1 = self.cleaned_data['set_input']
		value1_length = value1.split("\n")

		if len(value1_length) > 1000 and value1!='':
			raise forms.ValidationError("Input should be less than 1000")
		if len(value1_length) < 3 and value1!='':
			raise forms.ValidationError("Input should be more than 3")
		return value1

class AssocOptionForm(forms.Form):
	choice_option = [('n_category', 'By Disease Course'), ('pcs', 'By PCS'), ('pam50', 'By PAM50'), ('bcr', 'By BCR')]
	sample_option = forms.ChoiceField(widget=forms.RadioSelect, label='', choices=choice_option, required=True, initial=['n_category'])

class CorrOptionForm(forms.Form):
	choice_option = [('n_category', 'By Disease Course'), ('pcs', 'By PCS'), ('pam50', 'By PAM50')]
	sample_option = forms.ChoiceField(widget=forms.RadioSelect, label='', choices=choice_option, required=True, initial=['n_category'])

class SetOptionForm(forms.Form):
	choice_option = (('By Disease Course', (('n_category:GS<7', 'GS<7 VS Others'), ('n_category:GS=7', 'GS=7 VS Others'), ('n_category:GS>7', 'GS>7 VS Others'), ('n_category:mCRPC', 'mCRPC VS Others'))), ('By PCS', (('pcs:PCS1', 'PCS1 VS Others'), ('pcs:PCS2', 'PCS2 VS Others'), ('pcs:PCS3', 'PCS3 VS Others'))), ('By PAM50', (('pam50:LumA', 'Luminal A VS Others'), ('pam50:LumB', 'Luminal B VS Others'), ('pam50:Basal', 'Basal VS Others'))))
	sample_option = forms.ChoiceField(label='',widget=forms.Select,required=True, choices=choice_option)


class Pathway_option(forms.Form):

	pathways = PathwayInfo.objects.all().values('pathway_name').distinct().values_list('pathway_name', flat=True)
	pathways = [str(x) for x in pathways if str(x) is not 'NULL']

	pathway_choice = [(i,i) for i in pathways]
	pathway_select = forms.MultipleChoiceField(widget=forms.SelectMultiple(attrs={'size':'20','style': 'font-size: 12px'}), label='Select Pathways ', choices=pathway_choice, required=True)
	#pathway_select = forms.ModelChoiceField(widget=forms.Select(attrs={'size':'20','style': 'font-size: 12px'}), label='Select Pathways ', choices=pathway_choice, required=True)

	def clean(self):
		value = self.cleaned_data['pathway_select']
		if len(value) != 1:
			raise forms.ValidationError("You need to select 1 item only.")
		return value

class CorrInputName(forms.Form):
	name_input1 = forms.CharField(label='Name 1 ',initial='Input 1', required=False)
	name_input2 = forms.CharField(label='Name 2 ',initial='Input 2', required=False)

class SetInputName(forms.Form):
	set_name_input = forms.CharField(label='Name ',initial='Input', required=False)