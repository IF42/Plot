#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "../src/plot.h"
#include "../src/supportLib.h"

#define points 50



int
main(void)
{
	bool success;
	
	double xs[] = {-2, -1, 0, 1, 2};
	double ys[] = {-2, 1, 2, 1, -2};


	RGBABitmapImageReference *imageRef = CreateRGBABitmapImageReference();
	success = DrawScatterPlot(imageRef, 800, 600, xs, 5, ys, 5);

	double x1 = MapXCoordinateAutoSettings(0, imageRef->image, xs, 5);
  double y1 = MapYCoordinateAutoSettings(0, imageRef->image, ys, 5);

  double x2 = MapXCoordinateAutoSettings(1, imageRef->image, xs, 5);
  double y2 = MapYCoordinateAutoSettings(1, imageRef->image, ys, 5);

  DrawLine(imageRef->image, x1, y1, x2, y2, 2, GetGray(0.3));	


	if(success)
	{
		ByteArray *pngdata = ConvertToPNG(imageRef->image);
		DeleteImage(imageRef->image);

		WriteToFile(pngdata, "example3.png");
	}

    printf("Programe exit..\n");
    return EXIT_SUCCESS;
}
