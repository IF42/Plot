#include <stdio.h>
#include <stdlib.h>

#include <gtk/gtk.h>
#include "../src/gtk_plot.h"


int
main(int argc, char ** argv)
{
    float x[] = {1, 2, 3, 4};
    float y[] = {2, 4, 6, 8};
    
    gtk_init(&argc, &argv);
    
    GtkWidget * window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    GtkWidget * plot   = gtk_plot_new();
   
    gtk_plot_line(GTK_PLOT(plot), 4, x, x);

    gtk_container_add(GTK_CONTAINER(window), plot);
    
    g_signal_connect(window, "destroy", G_CALLBACK(gtk_main_quit), NULL);
    
    gtk_widget_show_all(window);
    
    gtk_main();

    printf("Program exit..\n");

    return EXIT_SUCCESS;
}
