# Generated by Django 3.1.2 on 2022-01-26 07:12

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('atel', '0002_atelsource_taskname'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='atelsource',
            name='atelID',
        ),
        migrations.AddField(
            model_name='atelmain',
            name='atelSrc',
            field=models.ManyToManyField(to='atel.AtelSource'),
        ),
    ]
