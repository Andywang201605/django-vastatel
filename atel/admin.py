from django.contrib import admin

from .models import AtelMain, AtelSource, AtelSub
# Register your models here.
admin.site.register(AtelMain)
admin.site.register(AtelSource)
admin.site.register(AtelSub)