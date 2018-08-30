from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('track_image/<int:p_cols>/<int:p_id>/<p_chrom>:<int:p_start>-<int:p_end>', views.image, name='image'),
    path('browser/<int:p_id>', views.browser, name='browser'),
    path('browser/<int:p_id>/<p_chrom>:<int:p_start>-<int:p_end>', views.browser, name='browser'),
    path('plots/', views.plots, name='plots'),
    path('plots/add', views.create_plot, name='create_plot'),
    path('plots/plot/<int:p_id>', views.edit_plot, name='edit_plot'),
    path('plots/plot/delete/<int:p_id>', views.delete_plot, name='delete_plot'),
    path('plots/plot/<int:p_plot_id>/add_track/', views.create_track, name='create_track'),
    path('datasources/', views.datasources, name='datasources'),
    path('datasources/datasource/add/<p_type>', views.create_datasource, name='add_datasource'),
    path('datasources/datasource/<int:p_id>', views.edit_datasource, name='edit_datasource'),
    path('datasources/datasource/delete/<int:p_id>', views.delete_datasource, name='delete_datasource'),
    path('track/<int:p_id>', views.edit_track, name='track'),
    path('track/delete/<int:p_id>', views.delete_track, name='delete_track'),
    path('evals', views.evals, name='evals'),
    path('evals/eval/add/<p_type>', views.create_eval, name='create_eval'),
    path('evals/eval/<int:p_id>', views.edit_eval, name='edit_eval'),
    path('evals/eval/delete/<int:p_id>', views.delete_eval, name='delete_eval'),
    path('evals/eval/<int:p_id>/show', views.show_eval, name='show_eval'),
    path('evals/eval/<int:p_id>/add', views.add_entry_eval, name='eval_add_entry'),
    path('ontologies/<p_db>', views.phenotypes, name='phenotypes'),
    path('ontologies/genes/', views.genes, name='genes'),
    path('ontologies/genes/<int:p_id>', views.gene, name='gene'),
]