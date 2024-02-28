#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <json.h>
#include <throw.h>

#include "../src/plot.h"


typedef struct {
    size_t size;
    double * x;
    double * y;
} Record;


Record record_init(size_t size) {
    return (Record) {size, malloc(sizeof(double) * size), malloc(sizeof(double) * size)};
}

void record_delete(Record * self) {
    if(self != NULL) {
        if(self->x != NULL)
            free(self->x);

        if(self->y != NULL)
            free(self->y);
    }
}


typedef struct {
    bool is_value;
    Record value;
}O_Record; 


O_Record __filter_record(Json * json, char * key, size_t start, size_t end) {
    Record record = record_init(end - start);

    for(size_t i = start; i < end && i < json->array.size; i++) {
        if(json_is_type(json_lookup(json->array.value[i], key), JsonFrac) == false) {
            debug("json dataset format error\n");
            record_delete(&record);
            return (O_Record) {.is_value = false};
        }

        record.y[i-start] = atof(json_lookup(json->array.value[i], key)->string);
        record.x[i-start] = i - start; 
    }

    return (O_Record) {.is_value = true, .value = record};
}


O_Record __load_price_record(char * filename) {
    FILE * f = fopen(filename, "r");

    if(f == NULL) {
        debug("Can't open dataset file: %s\n", filename);
        return (O_Record) {.is_value = false};
    }

    Json * json = json_parse(f);

    if(json_is_type(json, JsonArray) == false) {
        debug("Error reading json file: %s\n", filename);
        return (O_Record) {.is_value = false};
    }

	O_Record record = __filter_record(json, "open", 100, 500);

    json_delete(json);

    return record;
}


bool __show_candle_chart(void) {
	RGBABitmapImageReference * imageRef = CreateRGBABitmapImageReference();

    /*
     * loading chart data
     */ 
    O_Record record = __load_price_record("dataset/BITCOIN_M1_500.json");

    if(record.is_value == false) {
        return false;
	}

    /*
     * render price chart
     */
	ScatterPlotSeries open = GetDefaultScatterPlotSeriesSettings();
    open.xs = record.value.x;
    open.xsLength = record.value.size;
    open.ys = record.value.y;
    open.ysLength = record.value.size;
    open.linearInterpolation = true;
    open.lineType = Plot_LineType_Solid;
    open.lineThickness = 1;
    open.color = CreateRGBColor(0, 0.5, 0);

    ScatterPlotSettings settings = GetDefaultScatterPlotSettings();
    settings.width = 800;
    settings.height = 600;
    settings.autoBoundaries = false;
    settings.xMax = record.value.size;
    settings.xMin = 0;
    settings.yMax = 58000;
    settings.yMin = 55000;
    settings.autoPadding = true;
    settings.title = "BITCOIN";
    settings.xLabel = "time";
    settings.yLabel = "price";
    settings.showGrid = true;
    settings.scatterPlotSeries = (ScatterPlotSeries[]){open};
    settings.scatterPlotSeriesLength = 1;

    bool success = DrawScatterPlotFromSettings(imageRef, &settings); 

    /*
     * store plot as image
     */
    if(success == true) {
        ByteArray *pngdata = ConvertToPNG(imageRef->image);
        WriteToFile(pngdata, "plot.png");
		printf("chart rendered\n");
    } else {
        printf("error in chart rendering\n");
    }

    /*
     * clean memory
     */
    DeleteImage(imageRef->image);
	record_delete(&record.value);

    return true;
}


int
main(void)
{
    __show_candle_chart();
    printf("Programe exit..\n");
    return EXIT_SUCCESS;
}
