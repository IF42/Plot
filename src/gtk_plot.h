#ifndef _GTK_PLOT_H_
#define _GTK_PLOT_H_

#include <gtk/gtk.h>
#include <vector.h>


#define GTK_PLOT_TYPE_WIDGET (gtk_plot_get_type())


struct _GtkPlot;



G_DECLARE_FINAL_TYPE(GtkPlot, gtk_plot, GTK, PLOT, GtkDrawingArea)


GtkWidget *
gtk_plot_new(void);




void
gtk_plot_line(
    GtkPlot * self
    , size_t length
    , float * x
    , float * y);

/*
void
gtk_plot_line_plot_range(
    GtkPlot * self
    , size_t length
    , float begin
    , float * x
    , float * y);


void
gtk_plot_scatter_plot(
    GtkPlot * self
    , size_t length
    , float * x
    , float * y);


void
gtk_plot_scatter_plot_range(
    GtkPlot * self
    , size_t length
    , float begin
    , float * x
    , float * y);
*/


#endif
