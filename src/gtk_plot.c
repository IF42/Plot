#include "gtk_plot.h"

#include <gdk/gdk.h>


typedef float Coordinate[2];
typedef Vector(Coordinate)* Plot;


struct _GtkPlot
{
    GtkDrawingArea parrent;

    Vector(Plot) * scatter;
    Vector(Plot) * line;
};


G_DEFINE_TYPE(GtkPlot, gtk_plot, GTK_TYPE_DRAWING_AREA)


static void
gtk_plot_init(GtkPlot * widget)
{
    (void) widget;
}


static inline void
__plot_delete__(Vector(Plot) * plot)
{
    if(plot == NULL)
        return;

    for(size_t i = 0; i < VECTOR(plot)->length; i++)
        vector_delete(VECTOR(plot[i]));

    vector_delete(VECTOR(plot));
}


static void 
gtk_plot_dispose(GObject *object) 
{
    GtkPlot * widget = GTK_PLOT(object);
    
    if(widget->scatter != NULL)
    {
        __plot_delete__(widget->scatter);
        widget->scatter = NULL;
    }

    if(widget->line != NULL)
    {
        __plot_delete__(widget->line);
        widget->line = NULL;
    }

    printf("dispose\n");
    fflush(stdout);
    
    G_OBJECT_CLASS(gtk_plot_parent_class)->dispose(object);
}

static gboolean
gtk_plot_draw(
    GtkWidget * widget
    , cairo_t *cr) 
{
    int width = gtk_widget_get_allocated_width(widget);
    int height = gtk_widget_get_allocated_height(widget);

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    
    // Nastavíme tloušťku čáry
    cairo_set_line_width(cr, 1.0);
    
    cairo_rectangle(cr, 30, 30, width-60, height-60);   

    double x_step = (double)width / (VECTOR(GTK_PLOT(widget)->line[0])->length - 1);

    cairo_move_to(cr, 0, height / 2 - GTK_PLOT(widget)->line[0][0][0]);

    for(size_t i = 0; i < VECTOR(GTK_PLOT(widget)->line[0])->length; i++)
    {
        printf("%f\n", GTK_PLOT(widget)->line[0][i][0]);
        double x = i * x_step;
        double y = height / 2 - GTK_PLOT(widget)->line[0][i][0];
        cairo_line_to(cr, x, y);
    }

    // Provádíme kreslení
    cairo_stroke(cr);
    
    return FALSE;
}


static void
gtk_plot_class_init(GtkPlotClass * class)
{
    G_OBJECT_CLASS(class)->dispose = gtk_plot_dispose;
    GTK_WIDGET_CLASS(class)->draw = gtk_plot_draw;
}


GtkWidget *
gtk_plot_new(void)
{
	GtkWidget * widget = 
        GTK_WIDGET(g_object_new(GTK_PLOT_TYPE_WIDGET, NULL));

	GTK_PLOT(widget)->scatter = NULL;
	GTK_PLOT(widget)->line    = NULL;

	return widget;
}


void
gtk_plot_line(
    GtkPlot * self
    , size_t length
    , float * x
    , float * y)
{
    if(self->line == NULL)
        self->line = vector(Plot, 1);
    else
    {
        self->line = 
            vector_resize(
                VECTOR(self->line)
                , VECTOR(self->line)->length + 1);
    }

    Plot plot = vector(Coordinate, length);

    for(size_t i = 0; i < length; i++)
    {
        plot[i][0] = x[i];
        plot[i][1] = y[i];
    }

    self->line[VECTOR(self->line)->length - 1] = plot;
}


/*
static void
on_activate(GtkApplication * app)
{

}


Plotter * 
plotter_new(void)
{
     Plotter * self = malloc(sizeof(Plotter));

     self->scatter_plot = NULL;
     self->line_plot    = NULL;

     self->app = gtk_application_new(NULL, G_APPLICATION_DEFAULT_FLAGS);

     g_signal_connect(
         G_OBJECT(app)
         , "activate"
         , G_CALLBACK(view)
         , NULL);

     return self;
}


int
plotter_show(Plotter * self)
{
     int result = 
        g_application_run (
            G_APPLICATION (self->app)
            , 1, (char*[]){"plotter"});
    
     g_object_unref(self->app);

     return result;
}


void
plotter_line_plot(
    Plotter * self
    , size_t length
    , float * x
    , float * y)
{
    if(self->line_plot == NULL)
        self->line_plot = vector(Plot, 1);
    else
    {
        self->line_plot = 
            vector_resize(
                VECTOR(self->line_plot)
                , VECTOR(self->line_plot)->length + 1);
    }

    Plot plot = vector(Coordinate, length);

    for(size_t i = 0; i < length; i++)
    {
        plot[i][0] = x[i];
        plot[i][1] = y[i];
    }

    self->line_plot[VECTOR(self->line_plot)->length - 1] = plot;
}


void
plotter_line_plot_range(
    Plotter * self
    , size_t length
    , float begin
    , float * x
    , float * y)
{

}


void
plotter_scatter_plot(
    Plotter * self
    , size_t length
    , float * x
    , float * y)
{
    
}


void
plotter_scatter_plot_range(
    Plotter * self
    , size_t length
    , float begin
    , float * x
    , float * y)
{
    
}




void
plotter_delete(GtkPlot * self)
{
    if(self == NULL)
        return;

    __plot_delete__(self->scatter_plot);
    __plot_delete__(self->line_plot);

    free(self);
}

*/


