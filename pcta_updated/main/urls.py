from django.conf.urls import url
from . import views

from django.contrib.staticfiles.storage import staticfiles_storage
from django.views.generic.base import RedirectView


urlpatterns = [
	url(r'^$', views.main_template, name='main_template'),
	url(r'^option/$', views.option, name='option'),
	url(r'^calculation/$', views.call_calculation, name='calculation'),
	url(r'^about/$', views.about, name='about'),
	url(r'^manual/$', views.manual, name='manual'),
	url(r'^QnA/$', views.qna, name='qna'),
	url(r'^download/$', views.download, name='download'),
	url(r'^file_download/$', views.file_download, name='file_download'),
	url(r'^image_download/$', views.image_download, name='image_download'),
	url(r'^association_result/$', views.association_result, name='association_result'),
	url(r'^correlation_result/$', views.correlation_result, name='correlation_result'),
	url(r'^bcr_result/$', views.bcr_result, name='bcr_result'),
	url(r'^set_result/$', views.set_result, name='set_result'),
	url(r'^calculation/poll_state$', views.poll_state,name='poll_state'),
	url(r'^pathway_input$', views.initial_change, name='pathway_input'),

	url(r'^favicon.ico$', RedirectView.as_view( url=staticfiles_storage.url('images/favicon.ico'),), name="favicon"),
]