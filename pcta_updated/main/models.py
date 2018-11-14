# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from __future__ import unicode_literals
from django.db import models
from django.utils import timezone

class ExprInfo(models.Model):
    entrez_id = models.IntegerField()
    sample_id = models.IntegerField()
    value = models.FloatField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'expr_info'


class GeneInfo(models.Model):
    entrez_id = models.IntegerField()
    symbol = models.CharField(max_length=100)

    class Meta:
        managed = False
        db_table = 'gene_info'


class MajorGrade(models.Model):
    sample_id = models.IntegerField(primary_key=True)
    gleason = models.IntegerField(blank=True, null=True)
    disease_status = models.CharField(max_length=10, blank=True, null=True)
    hormone_response = models.CharField(max_length=10, blank=True, null=True)
    n_category = models.CharField(max_length=100, blank=True, null=True)
    pcs = models.CharField(db_column='PCS', max_length=100, blank=True, null=True)  # Field name made lowercase.
    pam50 = models.CharField(db_column='pam50',max_length=100, blank=True, null=True)
    sample_original_id = models.CharField(max_length=200, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'major_grade'


class MedInfo(models.Model):
    sample_id = models.IntegerField(primary_key=True)
    sample_original_id = models.CharField(max_length=200, blank=True, null=True)
    sampling_site = models.CharField(max_length=200, blank=True, null=True)
    histology = models.CharField(max_length=200, blank=True, null=True)
    invasiveness = models.CharField(max_length=10, blank=True, null=True)
    hormone_response = models.CharField(max_length=10, blank=True, null=True)
    tumor_content = models.IntegerField(blank=True, null=True)
    stromal_content = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'med_info'


class MinorGrade(models.Model):
    sample_id = models.IntegerField(primary_key=True)
    stage = models.IntegerField(blank=True, null=True)
    grade = models.IntegerField(blank=True, null=True)
    sample_original_id = models.CharField(max_length=200, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'minor_grade'


class MutStatusTable(models.Model):
    sample_id = models.IntegerField(primary_key=True)
    ets_expr = models.CharField(max_length=100, blank=True, null=True)
    egfr_mutation = models.CharField(max_length=100, blank=True, null=True)
    egfr_expr = models.CharField(max_length=100, blank=True, null=True)
    fusion = models.CharField(max_length=100, blank=True, null=True)
    sample_original_id = models.CharField(max_length=200, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'mut_status_table'


class PsaInfo(models.Model):
    sample_id = models.IntegerField(primary_key=True)
    psa_recurrence_month = models.IntegerField(blank=True, null=True)
    psa_recurrence_status = models.CharField(max_length=10, blank=True, null=True)
    psa_level = models.FloatField(blank=True, null=True)
    sample_original_id = models.CharField(max_length=200, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'psa_info'


class TimeInfo(models.Model):
    sample_id = models.IntegerField(primary_key=True)
    survival_month = models.IntegerField(blank=True, null=True)
    survival_status = models.CharField(max_length=10, blank=True, null=True)
    recurrence_month = models.IntegerField(blank=True, null=True)
    recurrence_status = models.CharField(max_length=10, blank=True, null=True)
    sample_original_id = models.CharField(max_length=200, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'time_info'


class UserInfo(models.Model):
    #user_id = models.CharField(db_column='user_ID', max_length=300, blank=True, null=False)  # Field name made lowercase.
    username = models.EmailField(db_column='user_ID', max_length=255, primary_key=True)
    org = models.CharField(max_length=300, blank=True, null=True)
    password1 = models.CharField(db_column='passwd',max_length=300, blank=True, null=False)

    class Meta:
        managed = False
        db_table = 'user_info'

class PathwayGenes(models.Model):
    pathway_name = models.CharField(max_length=300, blank=True, null=True)
    entrez_id = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'pathway_genes'


class PathwayInfo(models.Model):
    pathway_id = models.IntegerField(blank=True, null=True)
    pathway_name = models.CharField(max_length=300, blank=True, null=True)
    source1 = models.CharField(max_length=200, blank=True, null=True)
    source2 = models.CharField(max_length=100, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'pathway_info'

class News(models.Model):
    author = models.ForeignKey('auth.User')
    title = models.CharField(max_length=200)
    text = models.TextField()
    created_date = models.DateTimeField(default=timezone.now)
    published_date = models.DateTimeField(blank=True, null=True)

    def publish(self):
        self.published_date = timezone.now()
        self.save()

    def __str__(self):
        return self.title

