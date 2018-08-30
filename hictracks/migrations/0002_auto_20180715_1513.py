# Generated by Django 2.0.3 on 2018-07-15 15:13

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('hictracks', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='gene',
            name='chrom',
        ),
        migrations.RemoveField(
            model_name='gene',
            name='end',
        ),
        migrations.RemoveField(
            model_name='gene',
            name='id',
        ),
        migrations.RemoveField(
            model_name='gene',
            name='labels',
        ),
        migrations.RemoveField(
            model_name='gene',
            name='start',
        ),
        migrations.RemoveField(
            model_name='gene',
            name='strand',
        ),
        migrations.RemoveField(
            model_name='gene',
            name='track_file',
        ),
        migrations.AddField(
            model_name='gene',
            name='bedfileentry_ptr',
            field=models.OneToOneField(auto_created=True, default=1, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='hictracks.BedFileEntry'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='bedfileentry',
            name='track_file',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='file_entries', to='hictracks.TrackFile'),
        ),
        migrations.AlterField(
            model_name='genetophenotype',
            name='gene',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='hictracks.BedFileEntry'),
        ),
        migrations.AlterField(
            model_name='phenotype',
            name='genes',
            field=models.ManyToManyField(related_name='phenotypes', through='hictracks.GeneToPhenotype', to='hictracks.BedFileEntry'),
        ),
    ]