from django.db import models

# Create your models here.
class SourceWeb(models.Model):
    taskname = models.CharField(max_length=200)
    taskdate = models.DateTimeField('Date task was submitted')
    srcra = models.FloatField('R.A. of the source')
    srcdec = models.FloatField('Decl. of the source')
    fileExist = models.BooleanField('Whether the file is still existed')

    def __str__(self):
        return self.taskname
