#ifndef __PLOT_H__
#define __PLOT_H__

#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <interface/vector.h>


typedef union {
    uint32_t dword;

    struct {
        unsigned char B;
        unsigned char R;
        unsigned char G;
        unsigned char A;
    };    
}RGBA;


#define RGBA(...) (RGBA){__VA_ARGS__}

#define RGBA_Gray RGBA(.R=0xB0, .G=0xB0, .B=0xB0, .A=0xFF)
#define RGBA_Black RGBA(.R=0x0, .G=0x0, .B=0x0, .A=0xFF)
#define RGBA_White RGBA(.R=0xFF, .G=0xFF, .B=0xFF, .A=0xFF)


typedef struct {
    size_t size;
    size_t width;
    size_t height;

    uint32_t pixel[];
}Plot;



typedef enum {
    Plot_LineType_Solid
    , Plot_LineType_Dotted
}Plot_LineType;


typedef struct
{
    char * legenda;
    Plot_LineType line_type; 
    double line_thickness;
    const vector * xs;
    const vector * ys;
    RGBA color;
}ScatterPlot_Series;


typedef enum {
    ScatterPlot_X_Axis_Auto
    , ScatterPlot_X_Axis_Top
    , ScatterPlot_X_Axis_Bottom
}ScatterPlot_X_Axis;


typedef enum {
    ScatterPlot_Y_Axis_Auto
    , ScatterPlot_Y_Axis_Left
    , ScatterPlot_Y_Axis_Right
}ScatterPlot_Y_Axis;


typedef struct {
    size_t width;
    size_t height;

    char * x_label;
    char * y_label;
    char * title;

    bool show_grid;
    RGBA grid_color;

    bool padding_auto;
    float padding_left;
    float padding_right;
    float padding_top;
    float padding_bottom;

    size_t serie_size;
    ScatterPlot_Series ** serie;
}ScatterPlot_Settings;


Plot * scatter_plot_draw(const vector * xs, const vector * ys);
Plot * scatter_plot_draw_from_settings(ScatterPlot_Settings * settings);


void plot_delete(Plot * self);


typedef struct {
    vector vector;
    double start;
    double end;
    double step;
    double value;
    double iterator;
}Range;


Range range(double start, double end, double step);


#define range_to_vector(T) _Generic((T), Range*: ((const vector*) (T)))


#endif
