#include "plotter.h"

#include <gdk/gdk.h>


struct _GtkPlotter
{
    GtkDrawingArea parrent;

    Vector(Plot) * scatter_plot;
    Vector(Plot) * line_plot;
};


G_DEFINE_TYPE(GtkPlotter, gtk_plotter, GTK_TYPE_DRAWING_AREA)


static void
gtk_plotter_init(GtkPlotter * widget)
{
    (void) widget;
}


static void
gtk_plotter_class_init(GtkPlotterClass * class)
{
    (void)class;
}


static void
__draw_callback__(
	GtkPlotter * widget
	  , cairo_t * cr
	  , int width
	  , int height
	  , gpointer param)
{
	
}


GtkWidget *
gtk_plotter_new(void)
{
	GtkWidget * widget = 
        GTK_WIDGET(g_object_new(GTK_PLOTTER_TYPE_WIDGET, NULL));

	GTK_PLOTTER(widget)->scatter_plot = NULL;
	GTK_PLOTTER(widget)->line_plot    = NULL;

	gtk_drawing_area_set_draw_func(
            GTK_DRAWING_AREA(widget)
            , (GtkDrawingAreaDrawFunc) __draw_callback__
            , NULL
            , NULL);

	return widget;
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


static inline void
__plot_delete__(Vector(Plot) * plot)
{
    if(plot == NULL)
        return;

    for(size_t i = 0; i < VECTOR(plot)->length; i++)
        vector_delete(VECTOR(plot[i]));

    vector_delete(VECTOR(plot));
}


void
plotter_delete(GtkPlotter * self)
{
    if(self == NULL)
        return;

    __plot_delete__(self->scatter_plot);
    __plot_delete__(self->line_plot);

    free(self);
}

*/


