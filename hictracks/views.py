from django.shortcuts import render
from django.http import HttpResponse
from . import trackPlot
import matplotlib.pyplot as plt
from django.conf import settings
from django.template.response import TemplateResponse
from django.shortcuts import redirect
from django.core.exceptions import ObjectDoesNotExist
import re

from hictracks.models import Plot, Track, TrackFile, Assembly, Eval, Phenotype, Gene
from hictracks.misc import split_seq

from django.contrib import messages

from django.db.models import Q, ProtectedError

from io import StringIO, TextIOWrapper, BytesIO
from hictracks.readBed import BedOrBedGraphReader, ReadBedOrBedGraphException
from django.db import transaction

from .forms import CreateTrackForm, PlotForm, TrackForm, TrackFileForm, EvalForm, EvalAddEntryForm

from django_tables2 import RequestConfig
from .tables import PlotTable, TrackTable, TrackFileTable, EvalTable, EvalEntryTable, PhenotypeTable, GeneTable
from .tables import TrackFileFilter, EvalFilter, PlotFilter, PhenotypeFilter,  GeneFilter

from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin

from urllib.parse import urlencode, quote_plus

from collections import Counter
from scipy import stats


def is_object_readonly(request, obj):
    return not request.user.is_authenticated or request.user != obj.owner 

def only_public_or_user(request):
    user =  request.user if request.user.is_authenticated else None
    return Q(public = True) | Q(owner = user)


def get_chrom_length(plot_id, chrom):
    assembly = Plot.objects.get(id = plot_id).assembly
    return assembly.chromosomesassembly.chromosomes.get(name = chrom).size


def get_region_or_err_msg(search_text, plot_id):

    search_text  = re.sub(r'(\s|,)', '', search_text)

    if search_text == '':
        return 'Enter valid position query'

    sre = re.search('^chr(\d+|X|Y):(\d+)-(\d+)$',search_text)

    if not sre:
        return '"{search_text}" is not valid position query'.format(search_text = search_text)

    sre =  sre.groups()

    chrom = 'chr' + sre[0]
    start =  int(sre[1])
    end = int(sre[2])

    if start >= end:
        return 'End coordinate should be greater than start coordinate. Your position query: {search_text}'.format(search_text = search_text)


    try:    
        chrom_size = get_chrom_length(plot_id, chrom)
       
        if end > chrom_size:
            return 'End coordinate greater than chromosome size ({chrom_size}). Your position query: {search_text}'.format(chrom_size = chrom_size, search_text = search_text)

        if start > chrom_size:
            return 'Start coordinate greater then chromosome size ({chrom_size}). Your position query: {search_text}'.format(chrom_size = chrom_size, search_text = search_text)

    except ObjectDoesNotExist:
        return  'Chromosome "{chrom}" does not exists. Search text: {search_text}'.format(chrom = chrom, search_text = search_text)

    return  chrom, start, end

def move(start, end, perc, right):
    region_len =  end - start
    m =int(region_len  * perc)

    if right:
        return start + m, end + m, 100 * perc
    else:
        return max(start - m, 0), max(end - m, region_len), 100 * perc

def zoom(start, end, perc, zoom_in):

    if  zoom_in:
        perc = 1 / perc

    region_len =  end - start
    middle = start + region_len / 2

    new_region_len  = region_len * perc

    return max(int(middle - new_region_len / 2), 0), max(int(middle + new_region_len / 2), int(new_region_len ))


def setPlotCookie(request, p_id, p_chrom, p_start, p_end, p_interval_start, p_interval_end):
    request.session['plot_' + str(p_id)] = (p_chrom, p_start, p_end, p_interval_start, p_interval_end)
    request.session['plot'] = (p_id, p_chrom, p_start, p_end, p_interval_start, p_interval_end)

def getPlotCookie(request, p_id):
    if 'plot_' + str(p_id) in request.session:
        return request.session['plot_' + str(p_id)]
    return 'chr1', 30 * 1000 * 1000, 33 * 1000 * 1000 , None, None


def deletePlotCookie(request, p_id):
    if 'plot_' + str(p_id) in request.session:
        del request.session['plot_' + str(p_id)]

    if 'plot' in request.session and request.session['plot'][0] == str(p_id):
        del request.session['plot']

def printPlotCookie(request, p_id):
    if 'plot_' + str(p_id) in request.session:
        print(request.session['plot_' + str(p_id)])
    else:
        'Plot ' + str(p_id) + 'does not have cookie.'


def ranking(eval, p_chrom, p_interval_start, p_interval_end):
    p_interval_start, p_interval_end = int(p_interval_start), int(p_interval_end)

    region_start = max(p_interval_start - 3 * 1000 * 1000, 0)
    region_end = p_interval_end + 3 *1000 * 1000

    gene_file = TrackFile.objects.get(pk = 2)

    genes =  gene_file.get_entries(p_chrom, region_start,  region_end)

    genes = { gene.name : { 'gene': Gene.objects.get(pk =gene.id) } for gene in genes } 


    pLI_file = TrackFile.objects.get(pk = 6)
    pLI = { gene.name.upper() : gene.score for gene in pLI_file.get_entries(p_chrom, region_start, region_end) }

    clingen_file = TrackFile.objects.get(pk = 7)
    clingen = { gene.name : gene.score for gene in clingen_file.get_entries(p_chrom, region_start, region_end) }

    enh_prom_file = TrackFile.objects.get(pk = 40)
    enh_proms = [ enh_prom.name.upper() for enh_prom in enh_prom_file.get_entries(p_chrom, p_interval_start, p_interval_end)]
    enh_proms = Counter(enh_proms)


    for gene_name, d in genes.items():
        d['enh_prom'] = enh_proms[gene_name]

    enh_prom_scores = [ gene['enh_prom'] for gene in genes.values()]

    for gene_name, d in genes.items():
        gene = d['gene']
        d['gene_name'] = gene_name

        d['clingen'] = clingen.get(gene_name, None)
        d['clingen_score'] = 100 if clingen.get(gene_name, 0) in (3,2,30) else 0

        d['pLI'] = pLI.get(gene_name, None)

        d['distance'] = None

        d['distance'] = min( abs(p_interval_start - gene.start), 
                        abs(p_interval_end - gene.end))

        if d['distance'] < 1 * 1000 * 1000:
            d['distance_1Mb_score'] = 100
        else:
            d['distance_1Mb_score'] = 0

        d['phenotypes'] = gene.phenotypes.order_by('name', 'pheno_id')

        if d['phenotypes']:
            d['phenotype_score'] = 100
        else:
            d['phenotype_score'] = 0           

        if enh_proms[gene_name]:
            d['enh_prom_score'] = 100
        else:
            d['enh_prom_score'] = 0

        d['rank'] = d['clingen_score'] +  d['enh_prom_score'] +  d['phenotype_score'] + d['distance_1Mb_score'] 

    results = list(genes.items())
    results.sort(key = lambda x : (-x[1]['rank'], -x[1]['enh_prom'], x[1]['distance'], x[1]['gene_name']))

    """
    with open('/home/basia/CNVBrowser/latex/eval_' + str(eval.id) + '_' +  p_chrom + '.txt', 'w') as f_out:


        f_out.write(\\begin{center}
    \\begin{longtable}{ | l | l | l | l | l | l | l |}
  \\hline
   &  & &  Enhancer  &  & &\\\\ 
   &  & &  -promoter  & Distance  &  &\\\\ 
   Gene name  & pLI & Clingen & interactions & from & Phenotypes & Rank \\\\  
   &  &  & number & breakpoints & & \\\\  
  \\hline
"")
 
        for gene_name, gene in results:

            if gene['rank'] == 0:
                continue

            f_out.write(gene_name)
            f_out.write(' & ')
            f_out.write("{0:.4f}".format(gene['pLI']))
            f_out.write(' & ')
            f_out.write("{0:.0f}".format(gene['clingen']))
            f_out.write(' & ')
            f_out.write(str(gene['enh_prom']))
            f_out.write(' & ')
            f_out.write(str(gene['distance']))
            f_out.write(' & ')
            f_out.write('Yes' if gene['phenotypes'] else 'No')
            f_out.write(' & ')
            f_out.write(str(gene['rank']))
            f_out.write('\\\\\n')

        f_out.write('\\hline \\end{longtable} \\end{center}')

        f_out.write('\n\n')



        for gene_name, gene in results:

            phenotypes = gene['phenotypes']

            if not phenotypes:
                continue

            f_out.write('\\begin{center}\n')
            f_out.write('\t\\begin{longtable}{ | p{3cm} |  p{10cm} |}\n')

            f_out.write('\t\t\\hline\n')
            f_out.write('\\multicolumn{2}{ |c| }{' + gene_name + '} \\\\\n')

            f_out.write('\t\t\\hline\n')
            f_out.write('\t\t Name  & Definition and comments \\\\\n')
            f_out.write('\t\t\\hline\n')


            for phenotype in phenotypes:



                f_out.write(phenotype.name)
                f_out.write(' \href{' + phenotype.url + '}{' + phenotype.pheno_id +'}')
                f_out.write(' & ')

                info = []

                if phenotype.definition:
                    definition = 'Definition: ' + phenotype.definition.replace('%', '\%')
                    info.append(definition)
                if phenotype.comment:
                    comment = 'Comment: ' + phenotype.comment.replace('%', '\%')
                    info.append(comment)

                if phenotype.pheno_id in ('HP:0000006', 'HP:0000007'):
                    info = []

                if len(info) > 0:
                    f_out.write(info[0])

                f_out.write('\\\\\n')

                if len(info) == 2:
                    f_out.write(' & ')
                    f_out.write(info[1])                                         
                    f_out.write('\\\\\n')     

                f_out.write('\t\t\\hline\n')               
        
            f_out.write("\\end{longtable} \\end{center}")
"""
    return results    

def browser(request, p_id,  p_chrom = None, p_start = None, p_end = None):

    perc_move = (0.5, 0.25)
    perc_zoom = 1.25

    p_plot = Plot.objects.get(id = p_id)

    p_columns_dict = p_plot.getColumnsDict()
    
    c_chrom, c_start, c_end, c_interval_start, c_interval_end = getPlotCookie(request, p_id)

    if p_chrom is None:
        p_chrom = c_chrom
        p_start = c_start
        p_end = c_end


    p_name_filters = sorted(p_plot.getNameFilters(p_chrom, p_start, p_end))
    name_filter = ''
    get_url = ''
    results = ''

    if request.method == 'GET': 

        search_text = request.GET.get('search_text', '')
        name_filter = request.GET.get('name_filter', '')
        interval_start = request.GET.get('interval_start', c_interval_start or '')
        interval_end = request.GET.get('interval_end', c_interval_end or '')

        p_interval_start = request.GET.get('interval_start', c_interval_start)
        p_interval_end = request.GET.get('interval_end', c_interval_end)

        setPlotCookie(request, p_id, p_chrom, p_start, p_end, p_interval_start, p_interval_end)

        region_or_err_msg = get_region_or_err_msg(search_text, p_id)

        if len(region_or_err_msg) == 3:
            p_chrom, p_start, p_end = region_or_err_msg
        elif search_text  == '':
            pass
        else:
            messages.error(request, region_or_err_msg) 

        get_url = '?' +  urlencode({'name_filter':  name_filter, 
                                    'interval_start': interval_start, 
                                    'interval_end': interval_end}, quote_via=quote_plus)


        if p_plot.hasEval() and p_interval_start is not None and p_interval_end is not None:
            results = ranking(p_plot.eval, p_chrom, p_interval_start, p_interval_end)

    else:
        pass

    return TemplateResponse(request, 'hictracks/browser.html', {'p_id': p_id, 
                           'p_plot': p_plot,
                           'p_name_filters': p_name_filters,
                           'p_columns_dict': p_columns_dict,
                           'p_cols': len(p_columns_dict,),
                           'p_chrom': p_chrom, 
                           'p_start': p_start, 
                           'p_end': p_end,
                           'name_filter': name_filter,
                           'move': ( move( p_start, p_end, perc_move[0] , False), 
                                     move( p_start, p_end, perc_move[1] , False),
                                     move( p_start, p_end, perc_move[1] , True),
                                     move( p_start, p_end, perc_move[0] , True)
                                ),
                           'zoom_in': zoom( p_start, p_end, perc_zoom, True),
                           'zoom_out': zoom( p_start, p_end, perc_zoom, False),
                           'region_len': p_end - p_start,
                           'get_url':  get_url,
                           'results': results,
                           'interval_start': interval_start,
                           'interval_end': interval_end
            })

def image(request, p_cols, p_id, p_chrom, p_start, p_end):
    
    track = Track.objects.get(id = p_id)

    p_name_filter = None

    if request.method == 'GET': 
        name_filter = request.GET.get('name_filter', None)
        interval_start = request.GET.get('interval_start', '')

        if interval_start:
            interval_start = int(interval_start)

        interval_end = request.GET.get('interval_end', '')

        if interval_end:
            interval_end = int(interval_end)



    fig = track.draw_track(p_cols, chrom = p_chrom, start = p_start, end = p_end, 
        interval_start = interval_start, interval_end = interval_end, name_filter = name_filter)

    buf = BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)

    return HttpResponse(buf, content_type="image/png")

def plots(request):
    user =  request.user if request.user.is_authenticated else None

    plots = Plot.objects.filter(only_public_or_user(request))

   
    f = PlotFilter(request.GET, queryset=plots)
    table = PlotTable(f.qs)
    RequestConfig(request).configure(table)

    return render(request, 'hictracks/plots.html', {'table': table, 'filter': f})


def create_plot(request):

    if request.method == "POST":
        form = PlotForm(request.POST)

        if form.is_valid():
           plot = form.save(commit = False)
           plot.owner = request.user
           plot.save()

           p_id = plot.id

           return redirect(edit_plot, p_id= p_id)

    else:
        form = PlotForm()

    return render(request, 'hictracks/plot.html', {'form': form, 'assemblies': Assembly.objects.all() })


def edit_plot(request, p_id):

    plot = Plot.objects.get(pk =p_id)
    tracks = plot.tracks.all().order_by('column', 'no','id')
    table = TrackTable(tracks)
    RequestConfig(request).configure(table)

    if request.method == "POST":
        form = PlotForm(request.POST, instance = plot)

        if form.is_valid():
            plot = form.save()

            if 'add_track' in request.POST:
                return redirect(create_track, p_plot_id = p_id)

    else:
        form = PlotForm(instance = plot)

    return render(request, 'hictracks/plot.html', {'table': table,  'form': form,
                'p_id': p_id, 'assemblies': Assembly.objects.all() , 'has_tracks': len(tracks) > 0,
                'readonly': is_object_readonly(request, plot), 'plot' : plot})


def delete_plot(request, p_id):

    plot = Plot.objects.get(pk =p_id)

    try: 
     
        plot.delete()
        deletePlotCookie(request, p_id)

        messages.success(request, 'Plot "{}" successfully deleted.'.format(plot.name))
    except  ProtectedError:
        messages.error(request, 'You cannot delete this plot. It belongs to evaluation.')
        return redirect(edit_plot, p_id = p_id)

    return redirect(plots)


def datasources(request):

    track_files = TrackFile.objects.filter( ~Q(file_type= 'XA'),
                                            only_public_or_user(request),
                                            Q(eval__isnull=True) ).order_by('owner', 'assembly', 'file_type', 'name', 'id')

    f = TrackFileFilter(request.GET, queryset=track_files)
    table = TrackFileTable(f.qs)

    RequestConfig(request).configure(table)
    return render(request, 'hictracks/datasources.html', {'table': table, 'filter': f})



def edit_track(request, p_id):

    track = Track.objects.get(pk = p_id)

    if request.method == "POST":
        form = TrackForm(request.POST, instance = track)

        if form.is_valid():
            track = form.save()
        else:
            print(form.errors)
    else:
        
        form = TrackForm(instance = track)

    domains_files = None
    readonly = is_object_readonly(request, track.plot)

    if track.get_file_type() == 'HI':
        domains_files = TrackFile.objects.filter(file_type = 'BE').filter(only_public_or_user(request),  Q(eval__isnull=True)).order_by('name', 'id')


    return render(request, 'hictracks/track.html', {'form': form,
                   'p_id':  p_id, 
                   'plot_id': track.plot.id,
                   'file_type': track.get_file_type(), 
                   'file_sub_type': track.get_file_sub_type(),
                   'readonly': readonly,
                   'plot_br': track.plot.name,
                   'track_br': track.track_file.name,
                   'track_name': track.track_file.name,
                   'style_choices': track.get_style_choices(),
                   'domains_files': domains_files } )


def delete_track(request, p_id):

    track = Track.objects.get(pk =p_id)
    plot_id = track.plot.id
    track.delete()

    messages.success(request, 'Track successfully deleted.')
    return redirect(edit_plot, p_id = plot_id)

def create_track(request, p_plot_id):

    plot = Plot.objects.get(pk = p_plot_id)
    tracks = plot.tracks.all()

    track_numbers = [track.no for track in tracks] + [0]
    track_number = (max(track_numbers) + 10) // 10 * 10


    if request.method == "POST":
        form = CreateTrackForm(request.POST)

        if form.is_valid():
            track = form.save(commit = False)
            track.plot = plot
            track.save()

            messages.success(request, 'Track successfully created.')
            return redirect(edit_track, p_id = track.id)

    else:
        form = CreateTrackForm(initial={'no': track_number, 'column': 1})

    track_files = TrackFile.objects.filter(Q(assembly = plot.assembly), only_public_or_user(request), Q(eval__isnull=True)).order_by('name','id')

    return render(request, 'hictracks/add_track.html', {'form': form, 'p_plot_id': p_plot_id, 'track_files': track_files, 'plot': plot })

def edit_datasource(request, p_id):
    
    track_file = TrackFile.objects.get(pk = p_id)

    if request.method == 'POST':
        form = TrackFileForm(request.POST, instance = track_file)

        if form.is_valid():
            track_file = form.save()
            messages.success(request, 'Data source successfully saved.')
        else:
            messages.error(request, 'Data was NOT successfully saved.')            

            print(form['assembly'].errors)
            
    else:
        form = TrackFileForm(instance = track_file)

    return render(request, 'hictracks/datasource.html', {'form': form,  
                                                         'assemblies': Assembly.objects.all(),
                                                         'p_id': p_id,
                                                         'datasource_br': track_file.name,
                                                         'readonly': is_object_readonly(request, track_file )})


def delete_datasource(request, p_id):

    try: 
        track_file = TrackFile.objects.get(pk =p_id)
        track_file.delete()

        messages.success(request, 'Data source successfully deleted.')
    except ProtectedError:
        messages.error(request, 'You cannot delete this data source. It is used in a plot.')
        return redirect(edit_datasource, p_id= p_id)

    return redirect(datasources)


def save_datasource(track_file, file_handle):

    if file_handle:
        reader = BedOrBedGraphReader(file_handle = file_handle, track_file = track_file)

    track_file.save()

    if file_handle:
        for bed_entry in reader: 
            bed_entry.set_eval_pvalue()
            bed_entry.save()


@transaction.atomic
def save_datasource_atomic(track_file, file_handle):
    save_datasource(track_file, file_handle)


def  get_file_handle(p_type, form):

    if p_type == 'file':
        f = form.files['file']
        return  TextIOWrapper(f.file)            
    elif p_type == 'text' :
        text = form.data['text']
        if len(text) == 0:
            raise ReadBedOrBedGraphException('Paste in data in BED or BEDGraph format.')
        return StringIO(text)
    else:
        return None


def create_datasource(request, p_type):
    
    if request.method == 'POST':
        form = TrackFileForm(request.POST, request.FILES)

        if form.is_valid():

            try:

                track_file = form.save(commit = False)
                track_file.owner = request.user

                file_handle =  get_file_handle(p_type, form)
                save_datasource_atomic(track_file, file_handle)
            
                messages.success(request, "Data source '{}' successfully created.".format(track_file.name))

                return redirect(edit_datasource, p_id = track_file.id)

            except ReadBedOrBedGraphException as exp:
                messages.error(request, exp)   

    else:
        form = TrackFileForm()

    return render(request, 'hictracks/datasource.html', {'form': form,  'assemblies': Assembly.objects.all(), 'p_type': p_type, 'p_id': None})

def index(request):
    return render(request, 'hictracks/index.html')

def evals(request):
    table = None
    f = None

    if request.user.is_authenticated:
        evals = Eval.objects.filter(owner = request.user)
        table = EvalTable(evals)

        f = EvalFilter(request.GET, queryset=evals)
        table = EvalTable(f.qs)

        RequestConfig(request).configure(table)

    return render(request, 'hictracks/evals.html', {'table': table, 'filter':  f })

@transaction.atomic
def create_eval_atomic(request, form, p_type):
    track_file = TrackFile()
    track_file.owner = request.user
    assembly = form.instance.assembly
    track_file.assembly = assembly
    track_file.save()

    file_handle = get_file_handle(p_type, form)

    eval = form.save(commit = False)
    eval.owner = request.user
    eval.track_file =  track_file 


    columns = ( ( 1, 10001, 40, 6, ),)

    plot = Plot(assembly = assembly, owner = request.user)
    plot.title = 'Plot for evaluation \''  + form.cleaned_data['name'] + '\''
    plot.name = 'Plot for evaluation \''  + form.cleaned_data['name'] + '\''
    plot.save()

    
    for i, column in enumerate(columns):
        for j, track_id in enumerate(column): 

            if j == 0:
                 track = Track(plot = plot, 
                              track_file = TrackFile.objects.get(pk = track_id), 
                              no = (j + 1) * 10)               

            if j == 1:
                track = Track(plot = plot, 
                              track_file = TrackFile.objects.get(pk = track_id), 
                              no = (j + 1) * 10,
                              domains_file  = TrackFile.objects.get(pk = 10101))
            if j == 2:
                track = Track(plot = plot, 
                        track_file = TrackFile.objects.get(pk = track_id), 
                        name_filter = True, 
                        no = (j + 1) * 10, style = 'arcs')

            if j == 3:
                track = Track(plot = plot, 
                        track_file = TrackFile.objects.get(pk = track_id), 
                        no = (j + 1) * 10, style = 'tiles',
                        min_value = 0,
                        max_value = 1
                        )

            track.save()
    
    eval.plot = plot
    save_datasource(track_file, file_handle)    
    eval.save()

    return eval


def create_eval(request, p_type):

    if request.method == 'POST':
        form = EvalForm(request.POST, request.FILES)

        if form.is_valid():

            try:
                eval = create_eval_atomic(request, form, p_type)

                messages.success(request, "Successful evaluation.")

                return redirect(edit_eval, p_id = eval.id)

            except ReadBedOrBedGraphException as exp:
                messages.error(request, exp)   
    else:
        form = EvalForm()

    return render(request, 'hictracks/eval.html', {'form': form,  'assemblies': Assembly.objects.all(), 'p_type': p_type, 'p_id': None})

def edit_eval(request, p_id):
    
    eval = Eval.objects.get(pk = p_id)

    table = EvalEntryTable( eval.track_file.file_entries.all(), eval.plot.id)

    RequestConfig(request).configure(table)

    if request.method == 'POST':
        form = EvalForm(request.POST, instance = eval)

        if form.is_valid():
            eval = form.save()
            messages.success(request, 'Data source successfully saved.')
        else:
            print(form.errors)
            messages.error(request, 'Data was NOT successfully saved.')  
            
    else:
        form = TrackFileForm(instance = eval)

    return render(request, 'hictracks/eval.html', {'form': form,  
                                                 'p_id': p_id,
                                                 'assembly_name': eval.track_file.assembly.name , 
                                                 'eval_br': eval.name,
                                                 'table': table })

@transaction.atomic
def delete_eval(request, p_id):

    eval = Eval.objects.get(pk =p_id)

    eval.delete()

    eval.plot.delete()
    eval.track_file.delete()

    messages.success(request, "Evaluation of SVs '{}' successfully deleted.".format(eval.name))

    return redirect(evals)


def show_eval(request, p_id):

    eval = Eval.objects.get(pk =p_id) 

    tracks = eval.plot.getTracks()
 
    file_entries = split_seq(eval.track_file.file_entries.all(), 2)

    return render(request, 'hictracks/eval_show.html', {'p_id': p_id, 'tracks': tracks, 'file_entries' : file_entries} )

    

def add_entry_eval(request, p_id):

    eval = Eval.objects.get(pk =p_id) 

    chroms = eval.assembly.chromosomes.all()

    if request.method == 'POST':
        form = EvalAddEntryForm(request.POST)

        if form.is_valid():
            bed_file_entry = form.save(commit = False)
            bed_file_entry.track_file  = eval.track_file

            bed_file_entry.set_eval_pvalue()
            bed_file_entry.save()
            messages.success(request, "Breakpoint '{}' added to evaluation.".format(bed_file_entry.name))

            return redirect(edit_eval, p_id = eval.id)

        else:
            messages.error(request, 'Data was NOT successfully saved.')  
    else:
        form = EvalAddEntryForm()


    return render(request, 'hictracks/eval_add_entry.html', {'p_id': p_id, 'chroms': chroms , 'form': form} )





def phenotypes(request, p_db):

    phenotypes = Phenotype.objects.filter(db = p_db).order_by('pheno_id')


    f = PhenotypeFilter(request.GET, queryset=phenotypes)
    table = PhenotypeTable(f.qs)

    RequestConfig(request).configure(table)

    return render(request, 'hictracks/phenotypes.html', {'table': table, 'filter': f, 'db' : p_db })


def genes(request):

    genes = Gene.objects.all().order_by('name')

    f = GeneFilter(request.GET, queryset=genes)
    table = GeneTable(f.qs)

    RequestConfig(request).configure(table)

    return render(request, 'hictracks/genes.html', { 'table': table, 'filter': f})


def gene(request,p_id):

    gene = Gene.objects.get(pk=p_id)

    phenotypes = gene.phenotypes.all()
    f = PhenotypeFilter(request.GET, queryset=phenotypes)
    table = PhenotypeTable(f.qs)

    RequestConfig(request).configure(table)

    return render(request, 'hictracks/gene.html', {'table': table,'filter': f, 'gene':gene })
