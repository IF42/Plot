#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <json.h>
#include <throw.h>
#include <math.h>

#include "../src/plot.h"


#pragma pack(push, 1)
typedef struct {
    uint16_t type;            // BM
    uint32_t size;            // Velikost souboru
    uint16_t reserved1;
    uint16_t reserved2;
    uint32_t offset;          // Offset dat pixelů
    uint32_t header_size;     // Velikost hlavičky BMP
    int32_t  width;           // Šířka obrázku
    int32_t  height;          // Výška obrázku
    uint16_t planes;          // Počet barevných rovin (vždy 1)
    uint16_t bits_per_pixel;  // Počet bitů na pixel (24 pro true color)
    uint32_t compression;     // Typ komprese (0 pro nekomprimovaný BMP)
    uint32_t image_size;      // Velikost dat pixelů
    int32_t  x_pixels_per_meter;
    int32_t  y_pixels_per_meter;
    uint32_t colors_used;     // Počet barev v palete (0 pro true color)
    uint32_t important_colors;
} BMPHeader;
#pragma pack(pop)


bool write_as_bmp(Plot * image, char * path) {
	FILE * f = fopen(path, "w");

	if(f == NULL) {
		printf("Can't open file :%s\n", path);
		return false;
	}

	BMPHeader header = {
		.type = 0x4D42,            // BM
		.size = sizeof(BMPHeader) + (image->size * sizeof(uint32_t)),  // Velikost souboru
		.reserved1 = 0,
		.reserved2 = 0,
		.offset = sizeof(BMPHeader),  // Offset dat pixelů
		.header_size = 40,           // Velikost hlavičky BMP
		.width = image->width,
		.height = image->height,
		.planes = 1,
		.bits_per_pixel = 32,        // 24 bitů na pixel pro true color
		.compression = 0,
		.image_size = image->size * sizeof(uint32_t),  // Velikost dat pixelů
		.x_pixels_per_meter = 0,
		.y_pixels_per_meter = 0,
		.colors_used = 0,            // Pro true color není paleta
		.important_colors = 0
	 };

	fwrite(&header, 1, sizeof(BMPHeader), f);
	fwrite(image->pixel, 1, image->size * sizeof(uint32_t), f);
	fclose(f);

	return true;
}


typedef struct {
    size_t size;
    double * x;
    double * y;
} Batch;


Batch batch_init(size_t size) {
    return (Batch) {size, malloc(sizeof(double) * size), malloc(sizeof(double) * size)};
}

void batch_delete(Batch * self) {
    if(self != NULL) {
        if(self->x != NULL)
            free(self->x);

        if(self->y != NULL)
            free(self->y);
    }
}


typedef struct {
    bool is_value;
    Batch value;
}O_Batch; 


O_Batch __filter_batch(Json * json, char * key) {
    Batch record = batch_init(json->array.size);

    for(size_t i = 0; i < json->array.size; i++) {
        if(json_is_type(json_lookup(json->array.value[i], key), JsonFrac) == false) {
            debug("json dataset format error\n");
            batch_delete(&record);
            return (O_Batch) {.is_value = false};
        }

        record.y[i] = atof(json_lookup(json->array.value[i], key)->string);
        record.x[i] = i; 
    }

    return (O_Batch) {.is_value = true, .value = record};
}


O_Batch __load_price_batch(char * filename) {
    FILE * f = fopen(filename, "r");

    if(f == NULL) {
        debug("Can't open dataset file: %s\n", filename);
        return (O_Batch) {.is_value = false};
    }

    Json * json = json_parse(f);

    if(json_is_type(json, JsonArray) == false) {
        debug("Error reading json file: %s\n", filename);
        return (O_Batch) {.is_value = false};
    }

	O_Batch record = __filter_batch(json, "open");

    json_delete(json);

    return record;
}


Json * __load_dataset(char * filename) {
    FILE * f = fopen(filename, "r");

    if(f == NULL) {
        debug("Can't open dataset file\n");
		return NULL;
    }

    Json * json = json_parse(f);

    if(json_is_type(json, JsonArray) == false) {
        debug("Error reading json file\n");
    }

	return json;
}


bool __show_candle_chart(void) {
	Json * json = __load_dataset("dataset/BITCOIN_M1_500.json");

    if(json == NULL) {
        return false;
    }

	O_Batch batch_open = __filter_batch(json, "high");
	O_Batch batch_close = __filter_batch(json, "low");

    if(batch_open.is_value == false || batch_close.is_value == false) {
        return true;
	}

    ScatterPlot_Series serie_open = {
        .line_type = Plot_LineType_Solid
        , .size = batch_open.value.size
        , .xs = batch_open.value.x
        , .ys = batch_open.value.y
        , .legenda = "High price" 
        , .line_thickness = 1
        , .color = RGBA(.R=0xFF, .G=0x00, .B=0x20, .A=0xF0)};

    ScatterPlot_Series serie_close = {
        .line_type = Plot_LineType_Solid
        , .size = batch_close.value.size
        , .xs = batch_close.value.x
        , .ys = batch_close.value.y
        , .legenda = "Low price" 
        , .line_thickness = 1
        , .color = RGBA(.R=0x00, .G=0xFF, .B=0x20, .A=0xF0)};

    ScatterPlot_Settings settings = {
        .width = 1200
        , .height = 600
        , .x_label = "time"
        , .y_label = "price"
        , .title = "BITCOIN"
        , .padding_auto = true
        , .show_grid = true
        , .grid_color = RGBA_Gray
        , .serie_size = 2
        , .serie = (ScatterPlot_Series *[]) {&serie_open, &serie_close}
    };

    Plot * image = scatter_plot_draw_from_settings(&settings);

    /*
     * store plot as image
     */
    if(image != NULL) {
		write_as_bmp(image, "plot.bmp");
        plot_delete(image);
    } else {
        printf("error in chart rendering\n");
    }

	batch_delete(&batch_open.value);
	batch_delete(&batch_close.value);
    json_delete(json);

    return true;
}


int main(void) {
	__show_candle_chart();

    printf("Programe exit..\n");
    return EXIT_SUCCESS;
}





