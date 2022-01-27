from statistics import mode
from django.db import models

# Create your models here.
class AtelSub(models.Model):
    subject = models.CharField(max_length=200)

    def __str__(self):
        return self.subject

class AtelSource(models.Model):
    srcra = models.FloatField('R.A. of the source')
    srcdec = models.FloatField('Decl. of the source')
    taskname = models.CharField(max_length=200)

    def __str__(self):
        return f'AtelSrc#{self.pk}'

class AtelMain(models.Model):
    atelID = models.IntegerField('ATel ID', primary_key=True)
    atelSub = models.ManyToManyField(AtelSub)
    atelSrc = models.ManyToManyField(AtelSource)

    def __str__(self):
        return f'ATel#{self.atelID}'

