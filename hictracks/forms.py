from django import forms
from hictracks.models import Track, Plot, TrackFile, Eval, Assembly, BedFileEntry

class TrackForm(forms.ModelForm):

     class Meta:
        model = Track
        fields = ('height', 'colormap', 'min_value', 'max_value', 'style', 'column',
            'display', 'labels', 'inverted', 'x_labels', 'name_filter',
        'no', 'transform', 'color', 'edgecolor', 'title', 'bedgraph_style', 'bed_print_options', 'domains_file')


class CreateTrackForm(forms.ModelForm):

     class Meta:
        model = Track
        fields = ('no',  'title', 'height', 'track_file', 'column')


class PlotForm(forms.ModelForm):

     class Meta:
        model = Plot
        fields = ('title', 'name', 'assembly')


class TrackFileForm(forms.ModelForm):
    text = forms.CharField(widget=forms.Textarea)

    class Meta:
        model = TrackFile
        fields = ('assembly', 'name', 'text')

    def save(self, commit=True):
        return super(TrackFileForm, self).save(commit=commit)


class EvalForm(forms.ModelForm):
    text = forms.CharField(widget=forms.Textarea)

    class Meta:
        model = Eval
        fields = ('name', 'text', 'assembly')

    def save(self, commit=True):
        return super(EvalForm, self).save(commit=commit)

class EvalAddEntryForm(forms.ModelForm):
    class Meta:
        model = BedFileEntry
        fields = ('name', 'chrom', 'start', 'end')

    def clean(self):
        cleaned_data = super().clean()

        start = cleaned_data.get("start")
        end = cleaned_data.get("end")

        if start > end:
            raise forms.ValidationError(
                "The start coordinate should be less or equal to the end coordinate."
            )

        distance_limit =  10 * 1000 * 1000

        if end - start >  distance_limit:
            raise forms.ValidationError(
                "Distance from the start to the end coordinate must be less or equal to {:,}.".format(distance_limit)
            )