#ifndef _PLOTTER_H_
#define _PLOTTER_H_

#include <gtk/gtk.h>
#include <vector.h>


typedef float Coordinate[2];
typedef Vector(Coordinate)* Plot;


#define GTK_PLOTTER_TYPE_WIDGET (gtk_plotter_get_type())


struct _GtkPlotter;


G_DECLARE_FINAL_TYPE(GtkPlotter, gtk_plotter, GTK, PLOTTER, GtkDrawingArea)


GtkWidget *
gtk_plotter_new(void);



/*
void
gtk_plotter_line_plot(
    GtkPlotter * self
    , size_t length
    , float * x
    , float * y);


void
gtk_plotter_line_plot_range(
    GtkPlotter * self
    , size_t length
    , float begin
    , float * x
    , float * y);


void
gtk_plotter_scatter_plot(
    GtkPlotter * self
    , size_t length
    , float * x
    , float * y);


void
gtk_plotter_scatter_plot_range(
    GtkPlotter * self
    , size_t length
    , float begin
    , float * x
    , float * y);
*/


#endif
