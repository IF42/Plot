#include "plot.h"

#include <stdio.h>
#include <ctype.h>


static inline void aSwapElementsOfNumberArray(double *A, double ai, double bi){
    double tmp = A[(int)(ai)];

    A[(int)(ai)] = A[(int)(bi)];
    A[(int)(bi)] = tmp;
}


void WriteToFile(ByteArray *data, char *filename) {
	FILE* file = fopen(filename, "wb");
	fwrite(data->bytes, 1, data->bytesLength, file);
	fclose(file);
}


bool Loess(
        double *xs
        , size_t xsLength
        , double *ys
        , size_t ysLength
        , double bandwidth
        , double robustnessIters
        , double accuracy
        , NumberArrayReference *resultXs)
{
	double *weights;
	size_t weightsLength;

	weights = (double*)calloc(sizeof(double) * ((double)xsLength), 1);
	weightsLength = (double)xsLength;
	aFillNumberArray(weights, weightsLength, 1.0);

	return Lowess(
				xs
				, xsLength
				, ys
				, ysLength
				, weights
				, weightsLength
				, bandwidth
				, robustnessIters
				, accuracy
				, resultXs);
}


double Truncate(double x){
  double t;

  if(x >= 0.0){
    t = floor(x);
  }else{
    t = ceil(x);
  }

  return t;
}


static void RearrangeArray(double *as, size_t asLength, double *indexes, size_t indexesLength){
  double *bs;
  size_t bsLength;
  double i;

  bs = (double*)calloc(sizeof(double) * ((double)asLength), 1);
  bsLength = (double)asLength;

  AssignNumberArray(bs, bsLength, as, asLength);

  for(i = 0.0; i < (double)indexesLength; i = i + 1.0)
    as[(int)(i)] = bs[(int)(indexes[(int)(i)])];

  free(bs);
}


bool Lowess(
	double *xs
	, size_t xsLength
	, double *ys
	, size_t ysLength
	, double *weights
	, size_t weightsLength
	, double bandwidth
	, double robustnessIters
	, double accuracy
	, NumberArrayReference *resultXs)
{
  double *res, *residuals, *sortedResiduals, *robustnessWeights, *indexes;
  size_t resLength, residualsLength, sortedResidualsLength, robustnessWeightsLength, indexesLength;
  double n, i, k;
  double x, sumWeights, sumX, sumXSquared, sumY, sumXY, denom;
  double xk, yk, dist, w, xkw;
  double meanX, meanY, meanXY, meanXSquared;
  double alpha, beta;
  double arg, iter, medianResidual;
  double *bandwidthInterval;
  double ileft, iright, edge;
  double left, right, nextRight, nextLeft, bandwidthInPoints;
  bool success, done;

  /* Sort arrays */
  indexes = QuickSortNumbersWithIndexes(&indexesLength, xs, xsLength);
  RearrangeArray(ys, ysLength, indexes, indexesLength);

  if((double)xsLength == (double)ysLength && (double)xsLength != 0.0){
    n = (double)xsLength;

    if(n == 1.0 || n == 2.0){
      if(n == 1.0){
        res = (double*)calloc(sizeof(double) * (1.0), 1);
        resLength = 1.0;
        res[0] = ys[0];
      }else{
        res = (double*)calloc(sizeof(double) * (2.0), 1);
        resLength = 2.0;
        res[0] = ys[0];
        res[1] = ys[1];
      }

      resultXs->numberArray = res;
      resultXs->numberArrayLength = resLength;
      success = true;
    }else{
      bandwidthInPoints = Truncate(bandwidth*n);

      if(bandwidthInPoints >= 2.0){
        res = (double*)calloc(sizeof(double) * (n), 1);
        resLength = n;
        residuals = (double*)calloc(sizeof(double) * (n), 1);
        residualsLength = n;

        robustnessWeights = (double*)calloc(sizeof(double) * (n), 1);
        robustnessWeightsLength = n;
        aFillNumberArray(robustnessWeights, robustnessWeightsLength, 1.0);

        done = false;
        for(iter = 0.0; iter <= robustnessIters &&  !done ; iter = iter + 1.0){
          bandwidthInterval = (double*)calloc(sizeof(double) * (2.0), 1);
          bandwidthInterval[0] = 0.0;
          bandwidthInterval[1] = bandwidthInPoints - 1.0;

          for(i = 0.0; i < n; i = i + 1.0){
            x = xs[(int)(i)];

            if(i > 0.0){
              left = bandwidthInterval[0];
              right = bandwidthInterval[1];

              nextRight = FindNextNonZeroElement(weights, weightsLength, right);
              nextLeft = left;
              for(; nextRight < (double)xsLength && xs[(int)(nextRight)] - xs[(int)(i)] < xs[(int)(i)] - xs[(int)(nextLeft)]; ){
                nextLeft = FindNextNonZeroElement(weights, weightsLength, bandwidthInterval[0]);
                bandwidthInterval[0] = nextLeft;
                bandwidthInterval[1] = nextRight;
                nextRight = FindNextNonZeroElement(weights, weightsLength, nextRight);
              }
            }

            ileft = bandwidthInterval[0];
            iright = bandwidthInterval[1];

            if(xs[(int)(i)] - xs[(int)(ileft)] > xs[(int)(iright)] - xs[(int)(i)]){
              edge = ileft;
            }else{
              edge = iright;
            }

            sumWeights = 0.0;
            sumX = 0.0;
            sumXSquared = 0.0;
            sumY = 0.0;
            sumXY = 0.0;
            denom = fabs(1.0/(xs[(int)(edge)] - x));
            for(k = ileft; k <= iright; k = k + 1.0){
              xk = xs[(int)(k)];
              yk = ys[(int)(k)];

              if(k < i){
                dist = x - xk;
              }else{
                dist = xk - x;
              }

              w = Tricube(dist*denom)*robustnessWeights[(int)(k)]*weights[(int)(k)];
              xkw = xk*w;
              sumWeights = sumWeights + w;
              sumX = sumX + xkw;
              sumXSquared = sumXSquared + xk*xkw;
              sumY = sumY + yk*w;
              sumXY = sumXY + yk*xkw;
            }

            meanX = sumX/sumWeights;
            meanY = sumY/sumWeights;
            meanXY = sumXY/sumWeights;
            meanXSquared = sumXSquared/sumWeights;

            if(sqrt(fabs(meanXSquared - meanX*meanX)) < accuracy){
              beta = 0.0;
            }else{
              beta = (meanXY - meanX*meanY)/(meanXSquared - meanX*meanX);
            }

            alpha = meanY - beta*meanX;

            res[(int)(i)] = beta*x + alpha;

            residuals[(int)(i)] = fabs(ys[(int)(i)] - res[(int)(i)]);
          }

          if(iter == robustnessIters){
            done = true;
          }

          if( !done ){
            sortedResiduals = aCopyNumberArray(&sortedResidualsLength, residuals, residualsLength);
            QuickSortNumbers(sortedResiduals, sortedResidualsLength);

            medianResidual = sortedResiduals[(int)(n/2.0)];

            if(fabs(medianResidual) < accuracy){
              done = true;
            }

            if( !done ){
              for(i = 0.0; i < n; i = i + 1.0){
                arg = residuals[(int)(i)]/(6.0*medianResidual);
                if(arg >= 1.0){
                  robustnessWeights[(int)(i)] = 0.0;
                }else{
                  w = 1.0 - arg*arg;
                  robustnessWeights[(int)(i)] = w*w;
                }
              }
            }
          }
        }

        resultXs->numberArray = res;
        resultXs->numberArrayLength = resLength;
        success = true;
      }else{
        success = false;
        fprintf(stderr, "There must be at least two points.");
      }
    }
  }else{
    success = false;
    fprintf(stderr, "There must be equal number of points, and over zero.");
  }

  return success;
}


void AssignNumberArray(double *as, size_t asLength, double *bs, size_t bsLength){
  double i;

  for(i = 0.0; i < fmin((double)asLength, (double)bsLength); i = i + 1.0){
    as[(int)(i)] = bs[(int)(i)];
  }
}


double FindNextNonZeroElement(double *array, size_t arrayLength, double offset){
  double position;
  bool done;

  done = false;
  for(position = offset + 1.0; position < (double)arrayLength &&  !done ; position = position + 1.0){
    if(array[(int)(position)] != 0.0){
      done = true;
    }
  }

  return position;
}


double Tricube(double x){
  double ax, result;

  ax = fabs(x);

  if(ax >= 1.0){
    result = 0.0;
  }else{
    result = 1.0 - ax*ax*ax;
    result = result*result*result;
  }

  return result;
}


bool CropLineWithinBoundary(NumberReference *x1Ref, NumberReference *y1Ref, NumberReference *x2Ref, NumberReference *y2Ref, double xMin, double xMax, double yMin, double yMax){
  double x1, y1, x2, y2;
  bool success, p1In, p2In;
  double dx, dy, f1, f2, f3, f4, f;

  x1 = x1Ref->numberValue;
  y1 = y1Ref->numberValue;
  x2 = x2Ref->numberValue;
  y2 = y2Ref->numberValue;

  p1In = x1 >= xMin && x1 <= xMax && y1 >= yMin && y1 <= yMax;
  p2In = x2 >= xMin && x2 <= xMax && y2 >= yMin && y2 <= yMax;

  if(p1In && p2In){
    success = true;
  }else if( !p1In  && p2In){
    dx = x1 - x2;
    dy = y1 - y2;

    if(dx != 0.0){
      f1 = (xMin - x2)/dx;
      f2 = (xMax - x2)/dx;
    }else{
      f1 = 1.0;
      f2 = 1.0;
    }
    if(dy != 0.0){
      f3 = (yMin - y2)/dy;
      f4 = (yMax - y2)/dy;
    }else{
      f3 = 1.0;
      f4 = 1.0;
    }

    if(f1 < 0.0){
      f1 = 1.0;
    }
    if(f2 < 0.0){
      f2 = 1.0;
    }
    if(f3 < 0.0){
      f3 = 1.0;
    }
    if(f4 < 0.0){
      f4 = 1.0;
    }

    f = fmin(f1, fmin(f2, fmin(f3, f4)));

    x1 = x2 + f*dx;
    y1 = y2 + f*dy;

    success = true;
  }else if(p1In &&  !p2In ){
    dx = x2 - x1;
    dy = y2 - y1;

    if(dx != 0.0){
      f1 = (xMin - x1)/dx;
      f2 = (xMax - x1)/dx;
    }else{
      f1 = 1.0;
      f2 = 1.0;
    }
    if(dy != 0.0){
      f3 = (yMin - y1)/dy;
      f4 = (yMax - y1)/dy;
    }else{
      f3 = 1.0;
      f4 = 1.0;
    }

    if(f1 < 0.0){
      f1 = 1.0;
    }
    if(f2 < 0.0){
      f2 = 1.0;
    }
    if(f3 < 0.0){
      f3 = 1.0;
    }
    if(f4 < 0.0){
      f4 = 1.0;
    }

    f = fmin(f1, fmin(f2, fmin(f3, f4)));

    x2 = x1 + f*dx;
    y2 = y1 + f*dy;

    success = true;
  }else{
    success = false;
  }

  x1Ref->numberValue = x1;
  y1Ref->numberValue = y1;
  x2Ref->numberValue = x2;
  y2Ref->numberValue = y2;

  return success;
}


RGBA **Get8HighContrastColors(size_t *returnArrayLength){
  RGBA **colors;
  size_t colorsLength;
  colors = (RGBA**)calloc(sizeof(RGBA) * (8.0), 1);
  colorsLength = 8.0;
  colors[0] = CreateRGBColor(3.0/256.0, 146.0/256.0, 206.0/256.0);
  colors[1] = CreateRGBColor(253.0/256.0, 83.0/256.0, 8.0/256.0);
  colors[2] = CreateRGBColor(102.0/256.0, 176.0/256.0, 50.0/256.0);
  colors[3] = CreateRGBColor(208.0/256.0, 234.0/256.0, 43.0/256.0);
  colors[4] = CreateRGBColor(167.0/256.0, 25.0/256.0, 75.0/256.0);
  colors[5] = CreateRGBColor(254.0/256.0, 254.0/256.0, 51.0/256.0);
  colors[6] = CreateRGBColor(134.0/256.0, 1.0/256.0, 175.0/256.0);
  colors[7] = CreateRGBColor(251.0/256.0, 153.0/256.0, 2.0/256.0);
  *returnArrayLength = colorsLength;
  return colors;
}


void DrawFilledRectangleWithBorder(RGBABitmapImage *image, double x, double y, double w, double h, RGBA *borderColor, RGBA *fillColor){
  if(h > 0.0 && w > 0.0){
    DrawFilledRectangle(image, x, y, w, h, fillColor);
    DrawRectangle1px(image, x, y, w, h, borderColor);
  }
}


RGBABitmapImageReference *CreateRGBABitmapImageReference(){
  RGBABitmapImageReference *reference;

  reference = (RGBABitmapImageReference *)calloc(sizeof(RGBABitmapImageReference), 1);
  reference->image = (RGBABitmapImage *)calloc(sizeof(RGBABitmapImage), 1);
  reference->image->pixels = (uint32_t*)calloc(sizeof(uint32_t), 1);
  reference->image->xLength = 0.0;
  reference->image->yLength = 0.0;

  return reference;
}


bool RectanglesOverlap(Rectangle *r1, Rectangle *r2){
  bool overlap;

  overlap = false;

  overlap = overlap || (r2->x1 >= r1->x1 && r2->x1 <= r1->x2 && r2->y1 >= r1->y1 && r2->y1 <= r1->y2);
  overlap = overlap || (r2->x2 >= r1->x1 && r2->x2 <= r1->x2 && r2->y1 >= r1->y1 && r2->y1 <= r1->y2);
  overlap = overlap || (r2->x1 >= r1->x1 && r2->x1 <= r1->x2 && r2->y2 >= r1->y1 && r2->y2 <= r1->y2);
  overlap = overlap || (r2->x2 >= r1->x1 && r2->x2 <= r1->x2 && r2->y2 >= r1->y1 && r2->y2 <= r1->y2);

  return overlap;
}


Rectangle *CreateRectangle(double x1, double y1, double x2, double y2){
  Rectangle *r;
  r = (Rectangle *)calloc(sizeof(Rectangle), 1);
  r->x1 = x1;
  r->y1 = y1;
  r->x2 = x2;
  r->y2 = y2;
  return r;
}


void CopyRectangleValues(Rectangle *rd, Rectangle *rs){
  rd->x1 = rs->x1;
  rd->y1 = rs->y1;
  rd->x2 = rs->x2;
  rd->y2 = rs->y2;
}


void DrawXLabelsForPriority(double p, double xMin, double oy, double xMax, double xPixelMin, double xPixelMax, NumberReference *nextRectangle, RGBA *gridLabelColor, RGBABitmapImage *canvas, double *xGridPositions, size_t xGridPositionsLength, StringArrayReference *xLabels, NumberArrayReference *xLabelPriorities, Rectangle **occupied, bool textOnBottom){
  bool overlap, currentOverlaps;
  double i, j, x, px, padding;
  char *text;
  size_t textLength;
  Rectangle *r;

  r = (Rectangle *)calloc(sizeof(Rectangle), 1);
  padding = 10.0;

  overlap = false;
  for(i = 0.0; i < (double)xLabels->stringArrayLength; i = i + 1.0){
    if(xLabelPriorities->numberArray[(int)(i)] == p){

      x = xGridPositions[(int)(i)];
      px = MapXCoordinate(x, xMin, xMax, xPixelMin, xPixelMax);
      text = xLabels->stringArray[(int)(i)]->string;
      textLength = xLabels->stringArray[(int)(i)]->stringLength;

      r->x1 = floor(px - GetTextWidth(textLength)/2.0);
      if(textOnBottom){
        r->y1 = floor(oy + 5.0);
      }else{
        r->y1 = floor(oy - 20.0);
      }
      r->x2 = r->x1 + GetTextWidth(textLength);
      r->y2 = r->y1 + GetTextHeight;

      /* Add padding */
      r->x1 = r->x1 - padding;
      r->y1 = r->y1 - padding;
      r->x2 = r->x2 + padding;
      r->y2 = r->y2 + padding;

      currentOverlaps = false;

      for(j = 0.0; j < nextRectangle->numberValue; j = j + 1.0){
        currentOverlaps = currentOverlaps || RectanglesOverlap(r, occupied[(int)(j)]);
      }

      if( !currentOverlaps  && p == 1.0){
        DrawText(canvas, r->x1 + padding, r->y1 + padding, text, gridLabelColor);

        CopyRectangleValues(occupied[(int)(nextRectangle->numberValue)], r);
        nextRectangle->numberValue = nextRectangle->numberValue + 1.0;
      }

      overlap = overlap || currentOverlaps;
    }
  }
  if( !overlap  && p != 1.0){
    for(i = 0.0; i < (double)xGridPositionsLength; i = i + 1.0){
      x = xGridPositions[(int)(i)];
      px = MapXCoordinate(x, xMin, xMax, xPixelMin, xPixelMax);

      if(xLabelPriorities->numberArray[(int)(i)] == p){
        text = xLabels->stringArray[(int)(i)]->string;
        textLength = xLabels->stringArray[(int)(i)]->stringLength;

        r->x1 = floor(px - GetTextWidth(textLength)/2.0);
        if(textOnBottom){
          r->y1 = floor(oy + 5.0);
        }else{
          r->y1 = floor(oy - 20.0);
        }
        r->x2 = r->x1 + GetTextWidth(textLength);
        r->y2 = r->y1 + GetTextHeight;

        DrawText(canvas, r->x1, r->y1, text, gridLabelColor);

        CopyRectangleValues(occupied[(int)(nextRectangle->numberValue)], r);
        nextRectangle->numberValue = nextRectangle->numberValue + 1.0;
      }
    }
  }
}


void DrawYLabelsForPriority(double p, double yMin, double ox, double yMax, double yPixelMin, double yPixelMax, NumberReference *nextRectangle, RGBA *gridLabelColor, RGBABitmapImage *canvas, double *yGridPositions, size_t yGridPositionsLength, StringArrayReference *yLabels, NumberArrayReference *yLabelPriorities, Rectangle **occupied, bool textOnLeft){
  bool overlap, currentOverlaps;
  double i, j, y, py, padding;
  char *text;
  size_t textLength;
  Rectangle *r;

  r = (Rectangle *)calloc(sizeof(Rectangle), 1);
  padding = 10.0;
  overlap = false;

  for(i = 0.0; i < (double)yLabels->stringArrayLength; i = i + 1.0){
    if(yLabelPriorities->numberArray[(int)(i)] == p){

      y = yGridPositions[(int)(i)];
      py = MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax);
      text = yLabels->stringArray[(int)(i)]->string;
      textLength = yLabels->stringArray[(int)(i)]->stringLength;

      if(textOnLeft){
        r->x1 = floor(ox - GetTextWidth(textLength) - 10.0);
      }else{
        r->x1 = floor(ox + 10.0);
      }
      r->y1 = floor(py - 6.0);
      r->x2 = r->x1 + GetTextWidth(textLength);
      r->y2 = r->y1 + GetTextHeight;

      /* Add padding */
      r->x1 = r->x1 - padding;
      r->y1 = r->y1 - padding;
      r->x2 = r->x2 + padding;
      r->y2 = r->y2 + padding;

      currentOverlaps = false;

      for(j = 0.0; j < nextRectangle->numberValue; j = j + 1.0){
        currentOverlaps = currentOverlaps || RectanglesOverlap(r, occupied[(int)(j)]);
      }

      /* Draw labels with priority 1 if they do not overlap anything else. */
      if( !currentOverlaps  && p == 1.0){
        DrawText(canvas, r->x1 + padding, r->y1 + padding, text, gridLabelColor);

        CopyRectangleValues(occupied[(int)(nextRectangle->numberValue)], r);
        nextRectangle->numberValue = nextRectangle->numberValue + 1.0;
      }

      overlap = overlap || currentOverlaps;
    }
  }
  if( !overlap  && p != 1.0){
    for(i = 0.0; i < (double)yGridPositionsLength; i = i + 1.0){
      y = yGridPositions[(int)(i)];
      py = MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax);

      if(yLabelPriorities->numberArray[(int)(i)] == p){
        text = yLabels->stringArray[(int)(i)]->string;
        textLength = yLabels->stringArray[(int)(i)]->stringLength;

        if(textOnLeft){
          r->x1 = floor(ox - GetTextWidth(textLength) - 10.0);
        }else{
          r->x1 = floor(ox + 10.0);
        }
        r->y1 = floor(py - 6.0);
        r->x2 = r->x1 + GetTextWidth(textLength);
        r->y2 = r->y1 + GetTextHeight;

        DrawText(canvas, r->x1, r->y1, text, gridLabelColor);

        CopyRectangleValues(occupied[(int)(nextRectangle->numberValue)], r);
        nextRectangle->numberValue = nextRectangle->numberValue + 1.0;
      }
    }
  }
}


bool EpsilonCompare(double a, double b, double epsilon){
  return fabs(a - b) < epsilon;
}


double *ComputeGridLinePositions(size_t *returnArrayLength, double cMin, double cMax, StringArrayReference *labels, NumberArrayReference *priorities){
  double *positions;
  size_t positionsLength;
  double cLength, p, pMin, pMax, pInterval, pNum, i, num, rem, priority, mode;

  cLength = cMax - cMin;

  p = floor(log10(cLength));
  pInterval = pow(10.0, p);
  /* gives 10-1 lines for 100-10 diff */
  pMin = ceil(cMin/pInterval)*pInterval;
  pMax = floor(cMax/pInterval)*pInterval;
  pNum = round((pMax - pMin)/pInterval + 1.0);

  mode = 1.0;

  if(pNum <= 3.0){
    p = floor(log10(cLength) - 1.0);
    /* gives 100-10 lines for 100-10 diff */
    pInterval = pow(10.0, p);
    pMin = ceil(cMin/pInterval)*pInterval;
    pMax = floor(cMax/pInterval)*pInterval;
    pNum = round((pMax - pMin)/pInterval + 1.0);

    mode = 4.0;
  }else if(pNum <= 6.0){
    p = floor(log10(cLength));
    pInterval = pow(10.0, p)/4.0;
    /* gives 40-5 lines for 100-10 diff */
    pMin = ceil(cMin/pInterval)*pInterval;
    pMax = floor(cMax/pInterval)*pInterval;
    pNum = round((pMax - pMin)/pInterval + 1.0);

    mode = 3.0;
  }else if(pNum <= 10.0){
    p = floor(log10(cLength));
    pInterval = pow(10.0, p)/2.0;
    /* gives 20-3 lines for 100-10 diff */
    pMin = ceil(cMin/pInterval)*pInterval;
    pMax = floor(cMax/pInterval)*pInterval;
    pNum = round((pMax - pMin)/pInterval + 1.0);

    mode = 2.0;
  }

  positions = (double*)calloc(sizeof(double) * (pNum), 1);
  positionsLength = pNum;
  labels->stringArray = (StringReference**)calloc(sizeof(StringReference) * (pNum), 1);
  labels->stringArrayLength = pNum;
  priorities->numberArray = (double*)calloc(sizeof(double) * (pNum), 1);
  priorities->numberArrayLength = pNum;

  for(i = 0.0; i < pNum; i = i + 1.0){
    num = pMin + pInterval*i;
    positions[(int)(i)] = num;

    /* Always print priority 1 labels. Only draw priority 2 if they can all be drawn. Then, only draw priority 3 if they can all be drawn. */
    priority = 1.0;

    /* Prioritize x.25, x.5 and x.75 lower. */
    if(mode == 2.0 || mode == 3.0){
      rem = fmod(fabs(round(num/pow(10.0, p - 2.0))), 100.0);

      priority = 1.0;
      if(rem == 50.0){
        priority = 2.0;
      }else if(rem == 25.0 || rem == 75.0){
        priority = 3.0;
      }
    }

    /* Prioritize x.1-x.4 and x.6-x.9 lower */
    if(mode == 4.0){
      rem = fmod(fabs(round(num/pow(10.0, p))), 10.0);

      priority = 1.0;
      if(rem == 1.0 || rem == 2.0 || rem == 3.0 || rem == 4.0 || rem == 6.0 || rem == 7.0 || rem == 8.0 || rem == 9.0){
        priority = 2.0;
      }
    }

    /* 0 has lowest priority. */
    if(EpsilonCompare(num, 0.0, pow(10.0, p - 5.0))){
      priority = 3.0;
    }

    priorities->numberArray[(int)(i)] = priority;

    /* The label itself. */
    labels->stringArray[(int)(i)] = (StringReference *)calloc(sizeof(StringReference), 1);
    if(p < 0.0){
      if(mode == 2.0 || mode == 3.0){
        num = RoundToDigits(num,  -(p - 1.0));
      }else{
        num = RoundToDigits(num,  -p);
      }
    }
    labels->stringArray[(int)(i)]->string = CreateStringDecimalFromNumber(&labels->stringArray[(int)(i)]->stringLength, num);
  }

  *returnArrayLength = positionsLength;

  return positions;
}


double MapYCoordinate(double y, double yMin, double yMax, double yPixelMin, double yPixelMax){
  double yLength, yPixelLength;

  yLength = yMax - yMin;
  yPixelLength = yPixelMax - yPixelMin;

  y = y - yMin;
  y = y*yPixelLength/yLength;
  y = yPixelLength - y;
  y = y + yPixelMin;
  return y;
}


double MapXCoordinate(double x, double xMin, double xMax, double xPixelMin, double xPixelMax){
  double xLength, xPixelLength;

  xLength = xMax - xMin;
  xPixelLength = xPixelMax - xPixelMin;

  x = x - xMin;
  x = x*xPixelLength/xLength;
  x = x + xPixelMin;
  return x;
}


double MapXCoordinateAutoSettings(double x, RGBABitmapImage *image, double *xs, size_t xsLength){
  return MapXCoordinate(x, GetMinimum(xs, xsLength), GetMaximum(xs, xsLength), GetDefaultPaddingPercentage()*ImageWidth(image), (1.0 - GetDefaultPaddingPercentage())*ImageWidth(image));
}


double MapYCoordinateAutoSettings(double y, RGBABitmapImage *image, double *ys, size_t ysLength){
  return MapYCoordinate(y, GetMinimum(ys, ysLength), GetMaximum(ys, ysLength), GetDefaultPaddingPercentage()*ImageHeight(image), (1.0 - GetDefaultPaddingPercentage())*ImageHeight(image));
}


double MapXCoordinateBasedOnSettings(double x, ScatterPlotSettings *settings){
  double xMin, xMax, xPadding, xPixelMin, xPixelMax;
  Rectangle *boundaries;

  boundaries = (Rectangle *)calloc(sizeof(Rectangle), 1);
  ComputeBoundariesBasedOnSettings(settings, boundaries);
  xMin = boundaries->x1;
  xMax = boundaries->x2;

  if(settings->autoPadding){
    xPadding = floor(GetDefaultPaddingPercentage()*settings->width);
  }else{
    xPadding = settings->xPadding;
  }

  xPixelMin = xPadding;
  xPixelMax = settings->width - xPadding;

  return MapXCoordinate(x, xMin, xMax, xPixelMin, xPixelMax);
}


double MapYCoordinateBasedOnSettings(double y, ScatterPlotSettings *settings){
  double yMin, yMax, yPadding, yPixelMin, yPixelMax;
  Rectangle *boundaries;

  boundaries = (Rectangle *)calloc(sizeof(Rectangle), 1);
  ComputeBoundariesBasedOnSettings(settings, boundaries);
  yMin = boundaries->y1;
  yMax = boundaries->y2;

  if(settings->autoPadding){
    yPadding = floor(GetDefaultPaddingPercentage()*settings->height);
  }else{
    yPadding = settings->yPadding;
  }

  yPixelMin = yPadding;
  yPixelMax = settings->height - yPadding;

  return MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax);
}


double GetDefaultPaddingPercentage(){
  return 0.10;
}


void DrawText(RGBABitmapImage *canvas, double x, double y, char *text, RGBA *color){
  double i, charWidth, spacing;

  charWidth = 8.0;
  spacing = 2.0;

  for(i = 0.0; text[(int)i] != 0; i = i + 1.0){
    DrawAsciiCharacter(canvas, x + i*(charWidth + spacing), y, text[(int)(i)], color);
  }
}


void DrawTextUpwards(RGBABitmapImage *canvas, double x, double y, char *text, size_t textLength, RGBA *color){
  RGBABitmapImage *buffer, *rotated;

  buffer = CreateImage(GetTextWidth(textLength), GetTextHeight, GetTransparent());
  DrawText(buffer, 0.0, 0.0, text, color);
  rotated = RotateAntiClockwise90Degrees(buffer);
  DrawImageOnImage(canvas, rotated, x, y);
  DeleteImage(buffer);
  DeleteImage(rotated);
}


bool DrawScatterPlot(RGBABitmapImageReference *canvasReference, double width, double height, double *xs, size_t xsLength, double *ys, size_t ysLength){
  ScatterPlotSettings settings = GetDefaultScatterPlotSettings();

  settings.width = width;
  settings.height = height;
  settings.scatterPlotSeries = (ScatterPlotSeries[1]) {0};//(ScatterPlotSeries**) calloc(sizeof(ScatterPlotSeries) * (1.0), 1);
  settings.scatterPlotSeriesLength = 1.0;
  settings.scatterPlotSeries[0] = GetDefaultScatterPlotSeriesSettings();
  settings.scatterPlotSeries[0].xs = xs;
  settings.scatterPlotSeries[0].xsLength = xsLength;
  settings.scatterPlotSeries[0].ys = ys;
  settings.scatterPlotSeries[0].ysLength = ysLength;

  return DrawScatterPlotFromSettings(canvasReference, &settings);
}


bool DrawScatterPlotFromSettings(RGBABitmapImageReference *canvasReference, ScatterPlotSettings *settings)
{
  double xMin, xMax, yMin, yMax, xLength, yLength, i, x, y, xPrev, yPrev, px, py, pxPrev, pyPrev, originX, originY, p, l, plot;
  Rectangle *boundaries;
  double xPadding, yPadding, originXPixels, originYPixels;
  double xPixelMin, yPixelMin, xPixelMax, yPixelMax, xLengthPixels, yLengthPixels, axisLabelPadding;
  NumberReference *nextRectangle, *x1Ref, *y1Ref, *x2Ref, *y2Ref, *patternOffset;
  bool prevSet, success;
  RGBA *gridLabelColor;
  RGBABitmapImage *canvas;
  double *xs, *ys;
  size_t xsLength;
  bool linearInterpolation;
  ScatterPlotSeries *sp;
  double *xGridPositions, *yGridPositions;
  size_t xGridPositionsLength, yGridPositionsLength;
  StringArrayReference *xLabels, *yLabels;
  NumberArrayReference *xLabelPriorities, *yLabelPriorities;
  Rectangle **occupied;
  size_t occupiedLength;
  bool *linePattern;
  size_t linePatternLength;
  bool originXInside, originYInside, textOnLeft, textOnBottom;
  double originTextX, originTextY, originTextXPixels, originTextYPixels, side, yaxis;

  canvas = CreateImage(settings->width, settings->height, GetWhite());
  patternOffset = CreateNumberReference(0.0);

  success = ScatterPlotFromSettingsValid(settings);

  if(success){

    boundaries = (Rectangle *)calloc(sizeof(Rectangle), 1);
    ComputeBoundariesBasedOnSettings(settings, boundaries);
    xMin = boundaries->x1;
    yMin = boundaries->y1;
    xMax = boundaries->x2;
    yMax = boundaries->y2;

    /* If zero, set to defaults. */
    if(xMin - xMax == 0.0){
      xMin = 0.0;
      xMax = 10.0;
    }

    if(yMin - yMax == 0.0){
      yMin = 0.0;
      yMax = 10.0;
    }

    xLength = xMax - xMin;
    yLength = yMax - yMin;

    if(settings->autoPadding){
      xPadding = floor(GetDefaultPaddingPercentage()*settings->width);
      yPadding = floor(GetDefaultPaddingPercentage()*settings->height);
    }else{
      xPadding = settings->xPadding;
      yPadding = settings->yPadding;
    }

    /* Draw title */
    DrawText(canvas, floor(settings->width/2.0 - GetTextWidth(strlen(settings->title))/2.0), floor(yPadding/3.0), settings->title, GetBlack());

    /* Draw grid */
    xPixelMin = xPadding;
    yPixelMin = yPadding;
    xPixelMax = settings->width - xPadding;
    yPixelMax = settings->height - yPadding;
    xLengthPixels = xPixelMax - xPixelMin;
    yLengthPixels = yPixelMax - yPixelMin;
    DrawRectangle1px(canvas, xPixelMin, yPixelMin, xLengthPixels, yLengthPixels, settings->gridColor);

    gridLabelColor = GetGray(0.5);

    xLabels = (StringArrayReference *)calloc(sizeof(StringArrayReference), 1);
    xLabelPriorities = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);
    yLabels = (StringArrayReference *)calloc(sizeof(StringArrayReference), 1);
    yLabelPriorities = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);
    xGridPositions = ComputeGridLinePositions(&xGridPositionsLength, xMin, xMax, xLabels, xLabelPriorities);
    yGridPositions = ComputeGridLinePositions(&yGridPositionsLength, yMin, yMax, yLabels, yLabelPriorities);


    if(settings->showGrid){
      /* X-grid */
      for(i = 0.0; i < (double)xGridPositionsLength; i = i + 1.0){
        x = xGridPositions[(int)(i)];
        px = MapXCoordinate(x, xMin, xMax, xPixelMin, xPixelMax);
        DrawLine1px(canvas, px, yPixelMin, px, yPixelMax, settings->gridColor);
      }

      /* Y-grid */
      for(i = 0.0; i < (double)yGridPositionsLength; i = i + 1.0){
        y = yGridPositions[(int)(i)];
        py = MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax);
        DrawLine1px(canvas, xPixelMin, py, xPixelMax, py, settings->gridColor);
      }
    }

    /* Compute origin information. */
    originYInside = yMin < 0.0 && yMax > 0.0;
    originY = 0.0;
    if(settings->xAxisAuto){
      if(originYInside){
        originY = 0.0;
      }else{
        originY = yMin;
      }
    }else{
if(settings->xAxisTop){
        originY = yMax;
      }
      if(settings->xAxisBottom){
        originY = yMin;
      }
    }
    originYPixels = MapYCoordinate(originY, yMin, yMax, yPixelMin, yPixelMax);

    originXInside = xMin < 0.0 && xMax > 0.0;
    originX = 0.0;
    if(settings->yAxisAuto){
      if(originXInside){
        originX = 0.0;
      }else{
        originX = xMin;
      }
    }else{
if(settings->yAxisLeft){
        originX = xMin;
      }
      if(settings->yAxisRight){
        originX = xMax;
      }
    }
    originXPixels = MapXCoordinate(originX, xMin, xMax, xPixelMin, xPixelMax);

    if(originYInside){
      originTextY = 0.0;
    }else{
      originTextY = yMin + yLength/2.0;
    }
    originTextYPixels = MapYCoordinate(originTextY, yMin, yMax, yPixelMin, yPixelMax);

    if(originXInside){
      originTextX = 0.0;
    }else{
      originTextX = xMin + xLength/2.0;
    }
    originTextXPixels = MapXCoordinate(originTextX, xMin, xMax, xPixelMin, xPixelMax);

    /* Labels */
    occupied = (Rectangle**)calloc(sizeof(Rectangle) * ((double)xLabels->stringArrayLength + (double)yLabels->stringArrayLength), 1);
    occupiedLength = (double)xLabels->stringArrayLength + (double)yLabels->stringArrayLength;

    for(i = 0.0; i < (double)occupiedLength; i = i + 1.0){
      occupied[(int)(i)] = CreateRectangle(0.0, 0.0, 0.0, 0.0);
    }
    nextRectangle = CreateNumberReference(0.0);

    /* x labels */
    for(i = 1.0; i <= 5.0; i = i + 1.0){
      textOnBottom = true;
      if( !settings->xAxisAuto  && settings->xAxisTop){
        textOnBottom = false;
      }
      DrawXLabelsForPriority(i, xMin, originYPixels, xMax, xPixelMin, xPixelMax, nextRectangle, gridLabelColor, canvas, xGridPositions, xGridPositionsLength, xLabels, xLabelPriorities, occupied, textOnBottom);
    }

    /* y labels */
    for(i = 1.0; i <= 5.0; i = i + 1.0){
      textOnLeft = true;
      if( !settings->yAxisAuto  && settings->yAxisRight){
        textOnLeft = false;
      }

      DrawYLabelsForPriority(i, yMin, originXPixels, yMax, yPixelMin, yPixelMax, nextRectangle, gridLabelColor, canvas, yGridPositions, yGridPositionsLength, yLabels, yLabelPriorities, occupied, textOnLeft);
    }

    /* Draw origin line axis titles. */
    axisLabelPadding = 20.0;

    /* x origin line */
    if(originYInside){
      DrawLine1px(canvas, round(xPixelMin), round(originYPixels), round(xPixelMax), round(originYPixels), GetBlack());
    }

    /* y origin line */
    if(originXInside){
      DrawLine1px(canvas, round(originXPixels), round(yPixelMin), round(originXPixels), round(yPixelMax), GetBlack());
    }


    /* Draw origin axis titles. */
    DrawTextUpwards(canvas, 10.0, floor(originTextYPixels - GetTextWidth(strlen(settings->yLabel))/2.0), settings->yLabel, strlen(settings->yLabel), GetBlack());
    DrawText(canvas, floor(originTextXPixels - GetTextWidth(strlen(settings->xLabel))/2.0), yPixelMax + axisLabelPadding, settings->xLabel, GetBlack());

    /* X-grid-markers */
    for(i = 0.0; i < (double)xGridPositionsLength; i = i + 1.0){
      x = xGridPositions[(int)(i)];
      px = MapXCoordinate(x, xMin, xMax, xPixelMin, xPixelMax);
      p = xLabelPriorities->numberArray[(int)(i)];
      l = 1.0;
      if(p == 1.0){
        l = 8.0;
      }else if(p == 2.0){
        l = 3.0;
      }
      side =  -1.0;
      if( !settings->xAxisAuto  && settings->xAxisTop){
        side = 1.0;
      }
      DrawLine1px(canvas, px, originYPixels, px, originYPixels + side*l, GetBlack());
    }

    /* Y-grid-markers */
    for(i = 0.0; i < (double)yGridPositionsLength; i = i + 1.0){
      y = yGridPositions[(int)(i)];
      py = MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax);
      p = yLabelPriorities->numberArray[(int)(i)];
      l = 1.0;
      if(p == 1.0){
        l = 8.0;
      }else if(p == 2.0){
        l = 3.0;
      }
      side = 1.0;
      if( !settings->yAxisAuto  && settings->yAxisRight){
        side =  -1.0;
      }
      DrawLine1px(canvas, originXPixels, py, originXPixels + side*l, py, GetBlack());
    }

    /* Draw points */
    for(plot = 0.0; plot < (double)settings->scatterPlotSeriesLength; plot = plot + 1.0){
      sp = &settings->scatterPlotSeries[(int)(plot)];

      xs = sp->xs;
      xsLength = sp->xsLength;
      ys = sp->ys;
      linearInterpolation = sp->linearInterpolation;

      x1Ref = (NumberReference *)calloc(sizeof(NumberReference), 1);
      y1Ref = (NumberReference *)calloc(sizeof(NumberReference), 1);
      x2Ref = (NumberReference *)calloc(sizeof(NumberReference), 1);
      y2Ref = (NumberReference *)calloc(sizeof(NumberReference), 1);
      if(linearInterpolation){
        prevSet = false;
        xPrev = 0.0;
        yPrev = 0.0;
        for(i = 0.0; i < (double)xsLength; i = i + 1.0){
          x = xs[(int)(i)];
          y = ys[(int)(i)];

          if(prevSet){
            x1Ref->numberValue = xPrev;
            y1Ref->numberValue = yPrev;
            x2Ref->numberValue = x;
            y2Ref->numberValue = y;

            success = CropLineWithinBoundary(x1Ref, y1Ref, x2Ref, y2Ref, xMin, xMax, yMin, yMax);

            if(success){
              pxPrev = floor(MapXCoordinate(x1Ref->numberValue, xMin, xMax, xPixelMin, xPixelMax));
              pyPrev = floor(MapYCoordinate(y1Ref->numberValue, yMin, yMax, yPixelMin, yPixelMax));
              px = floor(MapXCoordinate(x2Ref->numberValue, xMin, xMax, xPixelMin, xPixelMax));
              py = floor(MapYCoordinate(y2Ref->numberValue, yMin, yMax, yPixelMin, yPixelMax));

              switch(sp->lineType) {
                  case Plot_LineType_Solid:
                      if(sp->lineThickness == 1.0)
                          DrawLine1px(canvas, pxPrev, pyPrev, px, py, sp->color);
                      else 
                          DrawLine(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, sp->color);
                      break;
                  case Plot_LineType_Dashed:
                      linePattern = GetLinePattern1(&linePatternLength);
                      DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, linePatternLength, patternOffset, sp->color);
                      break;
                  case Plot_LineType_Dotted:
                      linePattern = GetLinePattern2(&linePatternLength);
                      DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, linePatternLength, patternOffset, sp->color);
                      break;
                  case Plot_LineType_Dotdash:
                      linePattern = GetLinePattern3(&linePatternLength);
                      DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, linePatternLength, patternOffset, sp->color);
                      break;
                  case Plot_LineType_Longdash:
                      linePattern = GetLinePattern4(&linePatternLength);
                      DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, linePatternLength, patternOffset, sp->color);
                      break;
                  case Plot_LineType_Twodash:
                      linePattern = GetLinePattern5(&linePatternLength);
                      DrawLineBresenhamsAlgorithmThickPatterned(canvas, pxPrev, pyPrev, px, py, sp->lineThickness, linePattern, linePatternLength, patternOffset, sp->color);
                      break;
              }
            }
          }

          prevSet = true;
          xPrev = x;
          yPrev = y;
        }
      }else{
        for(i = 0.0; i < (double)xsLength; i = i + 1.0){
          x = xs[(int)(i)];
          y = ys[(int)(i)];

          if(x > xMin && x < xMax && y > yMin && y < yMax){

            x = floor(MapXCoordinate(x, xMin, xMax, xPixelMin, xPixelMax));
            y = floor(MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax));

            switch(sp->pointType) {
                case Plot_PointType_Crosses:
                    DrawPixel(canvas, x, y, sp->color);
                    DrawPixel(canvas, x + 1.0, y, sp->color);
                    DrawPixel(canvas, x + 2.0, y, sp->color);
                    DrawPixel(canvas, x - 1.0, y, sp->color);
                    DrawPixel(canvas, x - 2.0, y, sp->color);
                    DrawPixel(canvas, x, y + 1.0, sp->color);
                    DrawPixel(canvas, x, y + 2.0, sp->color);
                    DrawPixel(canvas, x, y - 1.0, sp->color);
                    DrawPixel(canvas, x, y - 2.0, sp->color);
                    break;
                case Plot_PointType_Circles:
                    DrawCircle(canvas, x, y, 3.0, sp->color);
                    break;
                case Plot_PointType_Dots:
                    DrawFilledCircle(canvas, x, y, 3.0, sp->color);
                    break;
                case Plot_PointType_Triangles:
                    DrawTriangle(canvas, x, y, 3.0, sp->color);
                    break;
                case Plot_PointType_FilledTriangles:
                    DrawFilledTriangle(canvas, x, y, 3.0, sp->color);
                    break;
                case Plot_PointType_Pixels:
                    DrawPixel(canvas, x, y, sp->color);
                    break;
                case Plot_PointType_DotlineToXAxis:
                    DrawFilledCircle(canvas, x, y, 3.0, sp->color);
                    yaxis = floor(MapYCoordinate(0.0, yMin, yMax, yPixelMin, yPixelMax));
                    yaxis = fmin(fmax(yaxis, yPixelMin), yPixelMax);
                    DrawLine(canvas, x, y, x, yaxis, sp->lineThickness, sp->color);
                    break;
            }
         }
        }
      }
    }

    canvasReference->image = canvas;
  }

  return success;
}


void ComputeBoundariesBasedOnSettings(ScatterPlotSettings *settings, Rectangle *boundaries){
  ScatterPlotSeries *sp;
  double plot, xMin, xMax, yMin, yMax;

  if((double)settings->scatterPlotSeriesLength >= 1.0){
    xMin = GetMinimum(settings->scatterPlotSeries[0].xs, settings->scatterPlotSeries[0].xsLength)*1.05;
    xMax = GetMaximum(settings->scatterPlotSeries[0].xs, settings->scatterPlotSeries[0].xsLength)*1.05;
    yMin = GetMinimum(settings->scatterPlotSeries[0].ys, settings->scatterPlotSeries[0].ysLength)*1.05;
    yMax = GetMaximum(settings->scatterPlotSeries[0].ys, settings->scatterPlotSeries[0].ysLength)*1.05;
  }else{
    xMin =  -10.0;
    xMax = 10.0;
    yMin =  -10.0;
    yMax = 10.0;
  }

  if( !settings->autoBoundaries ){
    xMin = settings->xMin;
    xMax = settings->xMax;
    yMin = settings->yMin;
    yMax = settings->yMax;
  }else{
    for(plot = 1.0; plot < (double)settings->scatterPlotSeriesLength; plot = plot + 1.0){
      sp = &settings->scatterPlotSeries[(int)(plot)];

      xMin = fmin(xMin, GetMinimum(sp->xs, sp->xsLength));
      xMax = fmax(xMax, GetMaximum(sp->xs, sp->xsLength));
      yMin = fmin(yMin, GetMinimum(sp->ys, sp->ysLength));
      yMax = fmax(yMax, GetMaximum(sp->ys, sp->ysLength));
    }
  }

  boundaries->x1 = xMin;
  boundaries->y1 = yMin;
  boundaries->x2 = xMax;
  boundaries->y2 = yMax;
}


bool ScatterPlotFromSettingsValid(ScatterPlotSettings *settings){
  bool success, found;
  ScatterPlotSeries *series;
  double i;

  success = true;

  /* Check axis placement. */
  if( !settings->xAxisAuto ){
    if(settings->xAxisTop && settings->xAxisBottom){
      success = false;
      fprintf(stderr, "x-axis not automatic and configured to be both on top and on bottom.\n");
    }
    if( !settings->xAxisTop  &&  !settings->xAxisBottom ){
      success = false;
      fprintf(stderr, "x-axis not automatic and configured to be neither on top nor on bottom.\n");
    }
  }

  if( !settings->yAxisAuto ){
    if(settings->yAxisLeft && settings->yAxisRight){
      success = false;
      fprintf(stderr, "y-axis not automatic and configured to be both on top and on bottom.\n");
    }
    if( !settings->yAxisLeft  &&  !settings->yAxisRight ){
      success = false;
      fprintf(stderr, "y-axis not automatic and configured to be neither on top nor on bottom.\n");
    }
  }

  /* Check series lengths. */
  for(i = 0.0; i < (double)settings->scatterPlotSeriesLength; i = i + 1.0){
    series = &settings->scatterPlotSeries[(int)(i)];
    if((double)series->xsLength != (double)series->ysLength){
      success = false;
      fprintf(stderr, "x and y series must be of the same length.\n");
    }
    if((double)series->xsLength == 0.0){
      success = false;
      fprintf(stderr, "There must be data in the series to be plotted.\n");
    }
    if(series->linearInterpolation && (double)series->xsLength == 1.0){
      success = false;
      fprintf(stderr, "Linear interpolation requires at least two data points to be plotted.\n");
    }
  }

  /* Check bounds. */
  if( !settings->autoBoundaries ){
    if(settings->xMin >= settings->xMax){
      success = false;
      fprintf(stderr, "x min is higher than or equal to x max.\n");
    }
    if(settings->yMin >= settings->yMax){
      success = false;
      fprintf(stderr, "y min is higher than or equal to y max.\n");
    }
  }

  /* Check padding. */
  if( !settings->autoPadding ){
    if(2.0*settings->xPadding >= settings->width){
      success = false;
      fprintf(stderr, "The x padding is more then the width.\n");
    }
    if(2.0*settings->yPadding >= settings->height){
      success = false;
      fprintf(stderr, "The y padding is more then the height.");
    }
  }

  /* Check width and height. */
  if(settings->width < 0.0){
    success = false;
    fprintf(stderr, "The width is less than 0.\n");
  }
  if(settings->height < 0.0){
    success = false;
    fprintf(stderr, "The height is less than 0.\n");
  }

  /* Check point types. */
  for(i = 0.0; i < (double)settings->scatterPlotSeriesLength; i = i + 1.0){
    series = &settings->scatterPlotSeries[(int)(i)];

    if(series->lineThickness < 0.0){
      success = false;
      fprintf(stderr, "The line thickness is less than 0.\n");
    }

    if( !series->linearInterpolation ){
      /* Point type. */
      switch(series->pointType) {
          case Plot_PointType_Crosses:
          case Plot_PointType_Circles:
          case Plot_PointType_Dots:
          case Plot_PointType_Triangles:
          case Plot_PointType_FilledTriangles:
          case Plot_PointType_Pixels:
          case Plot_PointType_DotlineToXAxis:
              found = true;
              break;
          default:
              found = false;
      }

      if( !found ){
        success = false;
        fprintf(stderr, "The point type is unknown.\n");
      }
    }else{
      /* Line type. */
      switch(series->lineType) {
          case Plot_LineType_Solid:
          case Plot_LineType_Dashed:
          case Plot_LineType_Dotted:
          case Plot_LineType_Dotdash:
          case Plot_LineType_Longdash:
          case Plot_LineType_Twodash:
              found = true;
              break;
          default:
              found = false;
      }

      if( !found ){
        success = false;
        fprintf(stderr, "The line type is unknown.\n");
      }
    }
  }

  return success;
}


RGBABitmapImage *DrawBarPlotNoErrorCheck(double width, double height, double *ys, size_t ysLength){
  RGBABitmapImageReference *canvasReference = CreateRGBABitmapImageReference();
  DrawBarPlot(canvasReference, width, height, ys, ysLength);

  return canvasReference->image;
}


bool DrawBarPlot(RGBABitmapImageReference *canvasReference, double width, double height, double *ys, size_t ysLength){
  BarPlotSettings settings;

  settings = GetDefaultBarPlotSettings();
  settings.barPlotSeries = (BarPlotSeries[1]){0}; 
  settings.barPlotSeriesLength = 1.0;
  settings.barPlotSeries[0] = GetDefaultBarPlotSeriesSettings();
  settings.barPlotSeries[0].ys = ys;
  settings.barPlotSeries[0].ysLength = ysLength;
  settings.width = width;
  settings.height = height;

  return DrawBarPlotFromSettings(canvasReference, &settings);
}


bool DrawBarPlotFromSettings(RGBABitmapImageReference *canvasReference, BarPlotSettings *settings){
  double xPadding, yPadding;
  double xPixelMin, yPixelMin, yPixelMax, xPixelMax;
  double xLengthPixels, yLengthPixels;
  double s, n, y, x, w, h, yMin, yMax, b, i, py, yValue;
  RGBA **colors;
  size_t colorsLength;
  double *ys, *yGridPositions;
  size_t yGridPositionsLength;
  double yTop, yBottom, ss, bs;
  double groupSeparation, barSeparation, barWidth, textwidth;
  StringArrayReference *yLabels;
  NumberArrayReference *yLabelPriorities;
  Rectangle **occupied;
  size_t occupiedLength;
  NumberReference *nextRectangle;
  RGBA *gridLabelColor, *barColor;
  char *label;
  size_t labelLength;
  bool success;
  RGBABitmapImage *canvas;

  success = BarPlotSettingsIsValid(settings);

  if(success){
    canvas = CreateImage(settings->width, settings->height, GetWhite());

    ss = (double)settings->barPlotSeriesLength;
    gridLabelColor = GetGray(0.5);

    /* padding */
    if(settings->autoPadding){
      xPadding = floor(GetDefaultPaddingPercentage()*ImageWidth(canvas));
      yPadding = floor(GetDefaultPaddingPercentage()*ImageHeight(canvas));
    }else{
      xPadding = settings->xPadding;
      yPadding = settings->yPadding;
    }

    /* Draw title */
    DrawText(canvas, floor(ImageWidth(canvas)/2.0 - GetTextWidth(strlen(settings->title))/2.0), floor(yPadding/3.0), settings->title, GetBlack());
    DrawTextUpwards(canvas, 10.0, floor(ImageHeight(canvas)/2.0 - GetTextWidth(strlen(settings->yLabel))/2.0), settings->yLabel, strlen(settings->yLabel), GetBlack());

    /* min and max */
    if(settings->autoBoundaries){
      if(ss >= 1.0){
        yMax = GetMaximum(settings->barPlotSeries[0].ys, settings->barPlotSeries[0].ysLength)*1.05;
        yMin = fmin(0.0, GetMinimum(settings->barPlotSeries[0].ys, settings->barPlotSeries[0].ysLength))*1.05;

        for(s = 0.0; s < ss; s = s + 1.0){
          yMax = fmax(yMax, GetMaximum(settings->barPlotSeries[(int)(s)].ys, settings->barPlotSeries[(int)(s)].ysLength));
          yMin = fmin(yMin, GetMinimum(settings->barPlotSeries[(int)(s)].ys, settings->barPlotSeries[(int)(s)].ysLength));
        }
      }else{
        yMax = 10.0;
        yMin = 0.0;
      }
    }else{
      yMin = settings->yMin;
      yMax = settings->yMax;
    }

    /* boundaries */
    xPixelMin = xPadding;
    yPixelMin = yPadding;
    xPixelMax = ImageWidth(canvas) - xPadding;
    yPixelMax = ImageHeight(canvas) - yPadding;
    xLengthPixels = xPixelMax - xPixelMin;
    yLengthPixels = yPixelMax - yPixelMin;

    /* Draw boundary. */
    DrawRectangle1px(canvas, xPixelMin, yPixelMin, xLengthPixels, yLengthPixels, settings->gridColor);

    /* Draw grid lines. */
    yLabels = (StringArrayReference *)calloc(sizeof(StringArrayReference), 1);
    yLabelPriorities = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);
    yGridPositions = ComputeGridLinePositions(&yGridPositionsLength, yMin, yMax, yLabels, yLabelPriorities);

    if(settings->showGrid){
      /* Y-grid */
      for(i = 0.0; i < (double)yGridPositionsLength; i = i + 1.0){
        y = yGridPositions[(int)(i)];
        py = MapYCoordinate(y, yMin, yMax, yPixelMin, yPixelMax);
        DrawLine1px(canvas, xPixelMin, py, xPixelMax, py, settings->gridColor);
      }
    }

    /* Draw origin. */
    if(yMin < 0.0 && yMax > 0.0){
      py = MapYCoordinate(0.0, yMin, yMax, yPixelMin, yPixelMax);
      DrawLine1px(canvas, xPixelMin, py, xPixelMax, py, settings->gridColor);
    }

    /* Labels */
    occupied = (Rectangle**)calloc(sizeof(Rectangle) * ((double)yLabels->stringArrayLength), 1);
    occupiedLength = (double)yLabels->stringArrayLength;
    for(i = 0.0; i < (double)occupiedLength; i = i + 1.0){
      occupied[(int)(i)] = CreateRectangle(0.0, 0.0, 0.0, 0.0);
    }
    nextRectangle = CreateNumberReference(0.0);

    for(i = 1.0; i <= 5.0; i = i + 1.0){
      DrawYLabelsForPriority(i, yMin, xPixelMin, yMax, yPixelMin, yPixelMax, nextRectangle, gridLabelColor, canvas, yGridPositions, yGridPositionsLength, yLabels, yLabelPriorities, occupied, true);
    }

    /* Draw bars. */
    if(settings->autoColor){
      if( !settings->grayscaleAutoColor ){
        colors = Get8HighContrastColors(&colorsLength);
      }else{
        colors = (RGBA**)calloc(sizeof(RGBA) * (ss), 1);
        colorsLength = ss;
        if(ss > 1.0){
          for(i = 0.0; i < ss; i = i + 1.0){
            colors[(int)(i)] = GetGray(0.7 - (i/ss)*0.7);
          }
        }else{
          colors[0] = GetGray(0.5);
        }
      }
    }else{
      colors = (RGBA**)calloc(sizeof(RGBA) * (0.0), 1);
      colorsLength = 0.0;
    }

    /* distances */
    bs = (double)settings->barPlotSeries[0].ysLength;

    if(settings->autoSpacing){
      groupSeparation = ImageWidth(canvas)*0.05;
      barSeparation = ImageWidth(canvas)*0.005;
    }else{
      groupSeparation = settings->groupSeparation;
      barSeparation = settings->barSeparation;
    }

    barWidth = (xLengthPixels - groupSeparation*(bs - 1.0) - barSeparation*(bs*(ss - 1.0)))/(bs*ss);

    /* Draw bars. */
    b = 0.0;
    for(n = 0.0; n < bs; n = n + 1.0){
      for(s = 0.0; s < ss; s = s + 1.0){
        ys = settings->barPlotSeries[(int)(s)].ys;

        yValue = ys[(int)(n)];

        yBottom = MapYCoordinate(yValue, yMin, yMax, yPixelMin, yPixelMax);
        yTop = MapYCoordinate(0.0, yMin, yMax, yPixelMin, yPixelMax);

        x = xPixelMin + n*(groupSeparation + ss*barWidth) + s*(barWidth) + b*barSeparation;
        w = barWidth;

        if(yValue >= 0.0){
          y = yBottom;
          h = yTop - y;
        }else{
          y = yTop;
          h = yBottom - yTop;
        }

        /* Cut at boundaries. */
        if(y < yPixelMin && y + h > yPixelMax){
          y = yPixelMin;
          h = yPixelMax - yPixelMin;
        }else if(y < yPixelMin){
          y = yPixelMin;
          if(yValue >= 0.0){
            h = yTop - y;
          }else{
            h = yBottom - y;
          }
        }else if(y + h > yPixelMax){
          h = yPixelMax - y;
        }

        /* Get color */
        if(settings->autoColor){
          barColor = colors[(int)(s)];
        }else{
          barColor = settings->barPlotSeries[(int)(s)].color;
        }

        /* Draw */
        if(settings->barBorder){
          DrawFilledRectangleWithBorder(canvas, round(x), round(y), round(w), round(h), GetBlack(), barColor);
        }else{
          DrawFilledRectangle(canvas, round(x), round(y), round(w), round(h), barColor);
        }

        b = b + 1.0;
      }
      b = b - 1.0;
    }

    /* x-labels */
    for(n = 0.0; n < bs; n = n + 1.0){
      if(settings->autoLabels){
        label = CreateStringDecimalFromNumber(&labelLength, n + 1.0);
      }else{
        label = settings->xLabels[(int)(n)];
        labelLength = strlen(settings->xLabels[(int)(n)]);
      }

      textwidth = GetTextWidth(labelLength);

      x = xPixelMin + (n + 0.5)*(ss*barWidth + (ss - 1.0)*barSeparation) + n*groupSeparation - textwidth/2.0;

      DrawText(canvas, floor(x), ImageHeight(canvas) - yPadding + 20.0, label, gridLabelColor);

      b = b + 1.0;
    }

    canvasReference->image = canvas;
  }

  return success;
}


bool BarPlotSettingsIsValid(BarPlotSettings *settings){
  bool success, lengthSet;
  BarPlotSeries *series;
  double i, length;

  success = true;

  /* Check series lengths. */
  lengthSet = false;
  length = 0.0;
  for(i = 0.0; i < (double)settings->barPlotSeriesLength; i = i + 1.0){
    series = &settings->barPlotSeries[(int)(i)];

    if( !lengthSet ){
      length = (double)series->ysLength;
      lengthSet = true;
    }else if(length != (double)series->ysLength){
      success = false;
      fprintf(stderr, "The number of data points must be equal for all series.");
    }
  }

  /* Check bounds. */
  if( !settings->autoBoundaries ){
    if(settings->yMin >= settings->yMax){
      success = false;
      fprintf(stderr, "Minimum y lower than maximum y.");
    }
  }

  /* Check padding. */
  if( !settings->autoPadding ){
    if(2.0*settings->xPadding >= settings->width){
      success = false;
      fprintf(stderr, "Double the horizontal padding is larger than or equal to the width.");
    }
    if(2.0*settings->yPadding >= settings->height){
      success = false;
      fprintf(stderr, "Double the vertical padding is larger than or equal to the height.");
    }
  }

  /* Check width and height. */
  if(settings->width < 0.0){
    success = false;
    fprintf(stderr,  "Width lower than zero.");
  }
  if(settings->height < 0.0){
    success = false;
    fprintf(stderr, "Height lower than zero.");
  }

  /* Check spacing */
  if( !settings->autoSpacing ){
    if(settings->groupSeparation < 0.0){
      success = false;
      fprintf(stderr, "Group separation lower than zero.");
    }
    if(settings->barSeparation < 0.0){
      success = false;
      fprintf(stderr, "Bar separation lower than zero.");
    }
  }

  return success;
}


double GetMinimum(double *data, size_t dataLength){
  double i, minimum;

  minimum = data[0];
  for(i = 0.0; i < (double)dataLength; i = i + 1.0){
    minimum = fmin(minimum, data[(int)(i)]);
  }

  return minimum;
}


double GetMaximum(double *data, size_t dataLength){
  double i, maximum;

  maximum = data[0];
  for(i = 0.0; i < (double)dataLength; i = i + 1.0){
    maximum = fmax(maximum, data[(int)(i)]);
  }

  return maximum;
}


double RoundToDigits(double element, double digitsAfterPoint){
  return round(element*pow(10.0, digitsAfterPoint))/pow(10.0, digitsAfterPoint);
}


RGBA *GetBlack(){
  RGBA *black;
  black = (RGBA *)calloc(sizeof(RGBA), 1);
  black->a = 1.0;
  black->r = 0.0;
  black->g = 0.0;
  black->b = 0.0;
  return black;
}


RGBA *GetWhite(){
  RGBA *white;
  white = (RGBA *)calloc(sizeof(RGBA), 1);
  white->a = 1.0;
  white->r = 1.0;
  white->g = 1.0;
  white->b = 1.0;
  return white;
}


RGBA *GetTransparent(){
  RGBA *transparent;
  transparent = (RGBA *)calloc(sizeof(RGBA), 1);
  transparent->a = 0.0;
  transparent->r = 0.0;
  transparent->g = 0.0;
  transparent->b = 0.0;
  return transparent;
}


RGBA *GetGray(double percentage){
  RGBA *black;
  black = (RGBA *)calloc(sizeof(RGBA), 1);
  black->a = 1.0;
  black->r = 1.0 - percentage;
  black->g = 1.0 - percentage;
  black->b = 1.0 - percentage;
  return black;
}


RGBA *CreateRGBColor(double r, double g, double b){
  RGBA *color;
  color = (RGBA *)calloc(sizeof(RGBA), 1);
  color->a = 1.0;
  color->r = r;
  color->g = g;
  color->b = b;
  return color;
}


RGBA *CreateRGBAColor(double r, double g, double b, double a){
  RGBA *color;
  color = (RGBA *)calloc(sizeof(RGBA), 1);
  color->a = a;
  color->r = r;
  color->g = g;
  color->b = b;
  return color;
}


RGBABitmapImage *CreateImage(double w, double h, RGBA *color){
  RGBABitmapImage *image;
  double i, j;

  image = (RGBABitmapImage *)calloc(sizeof(RGBABitmapImage), 1);
  image->pixels = (uint32_t*)calloc(sizeof(uint32_t) * w * h, 1);
  image->xLength = w;
  image->yLength = h;
  for(i = 0.0; i < w; i = i + 1.0){
    for(j = 0.0; j < h; j = j + 1.0){
      SetPixel(image, i, j, color);
    }
  }

  return image;
}


void DeleteImage(RGBABitmapImage *image){
  free(image->pixels);
  free(image);
}


double ImageWidth(RGBABitmapImage *image){
  return image->xLength;
}


double ImageHeight(RGBABitmapImage *image){
  double height;

  if(ImageWidth(image) == 0.0){
    height = 0.0;
  }else{
    height = image->yLength;
  }

  return height;
}


void SetPixel(RGBABitmapImage *image, double x, double y, RGBA *color){
  if(x >= 0.0 && x < ImageWidth(image) && y >= 0.0 && y < ImageHeight(image)){
		int pixel = x + y * ImageWidth(image);
		int r = round(color->r * 255);
		int g = round(color->g * 255);
		int b = round(color->b * 255);
		int a = round(color->a * 255);

		image->pixels[pixel] = (r << 24) | (g << 16) | (b << 8) | a;
  }
}


void DrawPixel(RGBABitmapImage *image, double x, double y, RGBA *color){
  double ra, ga, ba, aa;
  double rb, gb, bb, ab;
  double ro, go, bo, ao;
	RGBA newColor;
	uint32_t oldColor;
	int pixel;

  if(x >= 0.0 && x < ImageWidth(image) && y >= 0.0 && y < ImageHeight(image)){
    ra = color->r;
    ga = color->g;
    ba = color->b;
    aa = color->a;

		pixel = x + y * ImageWidth(image);
		oldColor	= image->pixels[pixel];
		rb = ((oldColor >> 24) & 0xFF) / 255.0;
		gb = ((oldColor >> 16) & 0xFF) / 255.0;
		bb = ((oldColor >> 8) & 0xFF) / 255.0;
		ab = ((oldColor >> 0) & 0xFF) / 255.0;

    ao = CombineAlpha(aa, ab);

    ro = AlphaBlend(ra, aa, rb, ab, ao);
    go = AlphaBlend(ga, aa, gb, ab, ao);
    bo = AlphaBlend(ba, aa, bb, ab, ao);
		
    newColor.r = ro;
    newColor.g = go;
    newColor.b = bo;
    newColor.a = ao;

    SetPixel(image, x, y, &newColor);
  }
}


double CombineAlpha(double as, double ad){
  return as + ad*(1.0 - as);
}


double AlphaBlend(double cs, double as, double cd, double ad, double ao){
  return (cs*as + cd*ad*(1.0 - as))/ao;
}


void DrawHorizontalLine1px(RGBABitmapImage *image, double x, double y, double length, RGBA *color){
  double i;

  for(i = 0.0; i < length; i = i + 1.0){
    DrawPixel(image, x + i, y, color);
  }
}


void DrawVerticalLine1px(RGBABitmapImage *image, double x, double y, double height, RGBA *color){
  double i;

  for(i = 0.0; i < height; i = i + 1.0){
    DrawPixel(image, x, y + i, color);
  }
}


void DrawRectangle1px(RGBABitmapImage *image, double x, double y, double width, double height, RGBA *color){
  DrawHorizontalLine1px(image, x, y, width + 1.0, color);
  DrawVerticalLine1px(image, x, y + 1.0, height + 1.0 - 1.0, color);
  DrawVerticalLine1px(image, x + width, y + 1.0, height + 1.0 - 1.0, color);
  DrawHorizontalLine1px(image, x + 1.0, y + height, width + 1.0 - 2.0, color);
}


void DrawImageOnImage(RGBABitmapImage *dst, RGBABitmapImage *src, double topx, double topy){
  double y, x;

  for(y = 0.0; y < ImageHeight(src); y = y + 1.0){
    for(x = 0.0; x < ImageWidth(src); x = x + 1.0){
      if(topx + x >= 0.0 && topx + x < ImageWidth(dst) && topy + y >= 0.0 && topy + y < ImageHeight(dst)){
				RGBA o = GetImagePixelStruct(src, x, y);
        DrawPixel(dst, topx + x, topy + y, &o);
      }
    }
  }
}


void DrawLine1px(RGBABitmapImage *image, double x0, double y0, double x1, double y1, RGBA *color){
  XiaolinWusLineAlgorithm(image, x0, y0, x1, y1, color);
}


void XiaolinWusLineAlgorithm(RGBABitmapImage *image, double x0, double y0, double x1, double y1, RGBA *color){
  bool steep;
  double x, t, dx, dy, g, xEnd, yEnd, xGap, xpxl1, ypxl1, intery, xpxl2, ypxl2, olda;

  olda = color->a;

  steep = fabs(y1 - y0) > fabs(x1 - x0);

  if(steep){
    t = x0;
    x0 = y0;
    y0 = t;

    t = x1;
    x1 = y1;
    y1 = t;
  }
  if(x0 > x1){
    t = x0;
    x0 = x1;
    x1 = t;

    t = y0;
    y0 = y1;
    y1 = t;
  }

  dx = x1 - x0;
  dy = y1 - y0;
  g = dy/dx;

  if(dx == 0.0){
    g = 1.0;
  }

  xEnd = round(x0);
  yEnd = y0 + g*(xEnd - x0);
  xGap = OneMinusFractionalPart(x0 + 0.5);
  xpxl1 = xEnd;
  ypxl1 = floor(yEnd);
  if(steep){
    DrawPixel(image, ypxl1, xpxl1, SetBrightness(color, OneMinusFractionalPart(yEnd)*xGap));
    DrawPixel(image, ypxl1 + 1.0, xpxl1, SetBrightness(color, FractionalPart(yEnd)*xGap));
  }else{
    DrawPixel(image, xpxl1, ypxl1, SetBrightness(color, OneMinusFractionalPart(yEnd)*xGap));
    DrawPixel(image, xpxl1, ypxl1 + 1.0, SetBrightness(color, FractionalPart(yEnd)*xGap));
  }
  intery = yEnd + g;

  xEnd = round(x1);
  yEnd = y1 + g*(xEnd - x1);
  xGap = FractionalPart(x1 + 0.5);
  xpxl2 = xEnd;
  ypxl2 = floor(yEnd);
  if(steep){
    DrawPixel(image, ypxl2, xpxl2, SetBrightness(color, OneMinusFractionalPart(yEnd)*xGap));
    DrawPixel(image, ypxl2 + 1.0, xpxl2, SetBrightness(color, FractionalPart(yEnd)*xGap));
  }else{
    DrawPixel(image, xpxl2, ypxl2, SetBrightness(color, OneMinusFractionalPart(yEnd)*xGap));
    DrawPixel(image, xpxl2, ypxl2 + 1.0, SetBrightness(color, FractionalPart(yEnd)*xGap));
  }

  if(steep){
    for(x = xpxl1 + 1.0; x <= xpxl2 - 1.0; x = x + 1.0){
      DrawPixel(image, floor(intery), x, SetBrightness(color, OneMinusFractionalPart(intery)));
      DrawPixel(image, floor(intery) + 1.0, x, SetBrightness(color, FractionalPart(intery)));
      intery = intery + g;
    }
  }else{
    for(x = xpxl1 + 1.0; x <= xpxl2 - 1.0; x = x + 1.0){
      DrawPixel(image, x, floor(intery), SetBrightness(color, OneMinusFractionalPart(intery)));
      DrawPixel(image, x, floor(intery) + 1.0, SetBrightness(color, FractionalPart(intery)));
      intery = intery + g;
    }
  }

  color->a = olda;
}


double OneMinusFractionalPart(double x){
  return 1.0 - FractionalPart(x);
}


double FractionalPart(double x){
  return x - floor(x);
}


RGBA *SetBrightness(RGBA *color, double newBrightness){
  color->a = newBrightness;
  return color;
}


void DrawQuadraticBezierCurve(RGBABitmapImage *image, double x0, double y0, double cx, double cy, double x1, double y1, RGBA *color){
  double t, dt, dx, dy;
  NumberReference *xs, *ys, *xe, *ye;

  dx = fabs(x0 - x1);
  dy = fabs(y0 - y1);

  dt = 1.0/sqrt(pow(dx, 2.0) + pow(dy, 2.0));

  xs = (NumberReference *)calloc(sizeof(NumberReference), 1);
  ys = (NumberReference *)calloc(sizeof(NumberReference), 1);
  xe = (NumberReference *)calloc(sizeof(NumberReference), 1);
  ye = (NumberReference *)calloc(sizeof(NumberReference), 1);

  QuadraticBezierPoint(x0, y0, cx, cy, x1, y1, 0.0, xs, ys);
  for(t = dt; t <= 1.0; t = t + dt){
    QuadraticBezierPoint(x0, y0, cx, cy, x1, y1, t, xe, ye);
    DrawLine1px(image, xs->numberValue, ys->numberValue, xe->numberValue, ye->numberValue, color);
    xs->numberValue = xe->numberValue;
    ys->numberValue = ye->numberValue;
  }

  free(xs);
  free(ys);
  free(xe);
  free(ye);
}


void QuadraticBezierPoint(double x0, double y0, double cx, double cy, double x1, double y1, double t, NumberReference *x, NumberReference *y){
  x->numberValue = pow(1.0 - t, 2.0)*x0 + (1.0 - t)*2.0*t*cx + pow(t, 2.0)*x1;
  y->numberValue = pow(1.0 - t, 2.0)*y0 + (1.0 - t)*2.0*t*cy + pow(t, 2.0)*y1;
}


void DrawCubicBezierCurve(RGBABitmapImage *image, double x0, double y0, double c0x, double c0y, double c1x, double c1y, double x1, double y1, RGBA *color){
  double t, dt, dx, dy;
  NumberReference *xs, *ys, *xe, *ye;

  dx = fabs(x0 - x1);
  dy = fabs(y0 - y1);

  dt = 1.0/sqrt(pow(dx, 2.0) + pow(dy, 2.0));

  xs = (NumberReference *)calloc(sizeof(NumberReference), 1);
  ys = (NumberReference *)calloc(sizeof(NumberReference), 1);
  xe = (NumberReference *)calloc(sizeof(NumberReference), 1);
  ye = (NumberReference *)calloc(sizeof(NumberReference), 1);

  CubicBezierPoint(x0, y0, c0x, c0y, c1x, c1y, x1, y1, 0.0, xs, ys);
  for(t = dt; t <= 1.0; t = t + dt){
    CubicBezierPoint(x0, y0, c0x, c0y, c1x, c1y, x1, y1, t, xe, ye);
    DrawLine1px(image, xs->numberValue, ys->numberValue, xe->numberValue, ye->numberValue, color);
    xs->numberValue = xe->numberValue;
    ys->numberValue = ye->numberValue;
  }

  free(xs);
  free(ys);
  free(xe);
  free(ye);
}


void CubicBezierPoint(double x0, double y0, double c0x, double c0y, double c1x, double c1y, double x1, double y1, double t, NumberReference *x, NumberReference *y){
  x->numberValue = pow(1.0 - t, 3.0)*x0 + pow(1.0 - t, 2.0)*3.0*t*c0x + (1.0 - t)*3.0*pow(t, 2.0)*c1x + pow(t, 3.0)*x1;

  y->numberValue = pow(1.0 - t, 3.0)*y0 + pow(1.0 - t, 2.0)*3.0*t*c0y + (1.0 - t)*3.0*pow(t, 2.0)*c1y + pow(t, 3.0)*y1;
}


RGBABitmapImage *CopyImage(RGBABitmapImage *image){
  RGBABitmapImage *copy;

  copy = CreateImage(ImageWidth(image), ImageHeight(image), GetTransparent());

  memcpy(copy->pixels, image->pixels, image->xLength * image->yLength);

  return copy;
}


RGBA* GetImagePixel(RGBABitmapImage *image, double x, double y){
	RGBA *rgba;
	int pixel;
	uint32_t color;

	rgba = (RGBA *)calloc(sizeof(RGBA), 1);

	pixel = x + y * ImageWidth(image);

	color	= image->pixels[pixel];

	rgba->r = ((color >> 24) & 0xFF) / 255.0;
	rgba->g = ((color >> 16) & 0xFF) / 255.0;
	rgba->b = ((color >> 8) & 0xFF) / 255.0;
	rgba->a = ((color >> 0) & 0xFF) / 255.0;

  return rgba;
}


RGBA GetImagePixelStruct(RGBABitmapImage *image, double x, double y){
	RGBA rgba;
	int pixel;
	uint32_t color;

	pixel = x + y * ImageWidth(image);

	color	= image->pixels[pixel];

	rgba.r = (color >> 24) & 0xFF;
	rgba.g = (color >> 16) & 0xFF;
	rgba.b = (color >> 8) & 0xFF;
	rgba.a = (color >> 0) & 0xFF;

  return rgba;
}


void HorizontalFlip(RGBABitmapImage *img){
  double y, x;
  double tmp;
  RGBA c1, c2;

  for(y = 0.0; y < ImageHeight(img); y = y + 1.0){
    for(x = 0.0; x < ImageWidth(img)/2.0; x = x + 1.0){
      c1 = GetImagePixelStruct(img, x, y);
      c2 = GetImagePixelStruct(img, ImageWidth(img) - 1.0 - x, y);

      tmp = c1.a;
      c1.a = c2.a;
      c2.a = tmp;

      tmp = c1.r;
      c1.r = c2.r;
      c2.r = tmp;

      tmp = c1.g;
      c1.g = c2.g;
      c2.g = tmp;

      tmp = c1.b;
      c1.b = c2.b;
      c2.b = tmp;
    }
  }
}


void DrawFilledRectangle(RGBABitmapImage *image, double x, double y, double w, double h, RGBA *color){
  double i, j;

  for(i = 0.0; i < w; i = i + 1.0){
    for(j = 0.0; j < h; j = j + 1.0){
      SetPixel(image, x + i, y + j, color);
    }
  }
}


RGBABitmapImage *RotateAntiClockwise90Degrees(RGBABitmapImage *image){
  RGBABitmapImage *rotated;
  double x, y;

  rotated = CreateImage(ImageHeight(image), ImageWidth(image), GetBlack());

  for(y = 0.0; y < ImageHeight(image); y = y + 1.0){
    for(x = 0.0; x < ImageWidth(image); x = x + 1.0){
      SetPixel(rotated, y, ImageWidth(image) - 1.0 - x, GetImagePixel(image, x, y));
    }
  }

  return rotated;
}


void DrawCircle(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color){
  DrawCircleBasicAlgorithm(canvas, xCenter, yCenter, radius, color);
}


void BresenhamsCircleDrawingAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color){
  double x, y, delta;

  y = radius;
  x = 0.0;

  delta = 3.0 - 2.0*radius;
  for(; y >= x; x = x + 1.0){
    DrawLine1px(canvas, xCenter + x, yCenter + y, xCenter + x, yCenter + y, color);
    DrawLine1px(canvas, xCenter + x, yCenter - y, xCenter + x, yCenter - y, color);
    DrawLine1px(canvas, xCenter - x, yCenter + y, xCenter - x, yCenter + y, color);
    DrawLine1px(canvas, xCenter - x, yCenter - y, xCenter - x, yCenter - y, color);

    DrawLine1px(canvas, xCenter - y, yCenter + x, xCenter - y, yCenter + x, color);
    DrawLine1px(canvas, xCenter - y, yCenter - x, xCenter - y, yCenter - x, color);
    DrawLine1px(canvas, xCenter + y, yCenter + x, xCenter + y, yCenter + x, color);
    DrawLine1px(canvas, xCenter + y, yCenter - x, xCenter + y, yCenter - x, color);

    if(delta < 0.0){
      delta = delta + 4.0*x + 6.0;
    }else{
      delta = delta + 4.0*(x - y) + 10.0;
      y = y - 1.0;
    }
  }
}


void DrawCircleMidpointAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color){
  double d, x, y;

  d = floor((5.0 - radius*4.0)/4.0);
  x = 0.0;
  y = radius;

  for(; x <= y; x = x + 1.0){
    DrawPixel(canvas, xCenter + x, yCenter + y, color);
    DrawPixel(canvas, xCenter + x, yCenter - y, color);
    DrawPixel(canvas, xCenter - x, yCenter + y, color);
    DrawPixel(canvas, xCenter - x, yCenter - y, color);
    DrawPixel(canvas, xCenter + y, yCenter + x, color);
    DrawPixel(canvas, xCenter + y, yCenter - x, color);
    DrawPixel(canvas, xCenter - y, yCenter + x, color);
    DrawPixel(canvas, xCenter - y, yCenter - x, color);

    if(d < 0.0){
      d = d + 2.0*x + 1.0;
    }else{
      d = d + 2.0*(x - y) + 1.0;
      y = y - 1.0;
    }
  }
}


void DrawCircleBasicAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color){
  double pixels, a, da, dx, dy;

  /* Place the circle in the center of the pixel. */
  xCenter = floor(xCenter) + 0.5;
  yCenter = floor(yCenter) + 0.5;

  pixels = 2.0*M_PI*radius;

  /* Below a radius of 10 pixels, over-compensate to get a smoother circle. */
  if(radius < 10.0){
    pixels = pixels*10.0;
  }

  da = 2.0*M_PI/pixels;

  for(a = 0.0; a < 2.0*M_PI; a = a + da){
    dx = cos(a)*radius;
    dy = sin(a)*radius;

    /* Floor to get the pixel coordinate. */
    DrawPixel(canvas, floor(xCenter + dx), floor(yCenter + dy), color);
  }
}


void DrawFilledCircle(RGBABitmapImage *canvas, double x, double y, double r, RGBA *color){
  DrawFilledCircleBasicAlgorithm(canvas, x, y, r, color);
}


void DrawFilledCircleMidpointAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color){
  double d, x, y;

  d = floor((5.0 - radius*4.0)/4.0);
  x = 0.0;
  y = radius;

  for(; x <= y; x = x + 1.0){
    DrawLineBresenhamsAlgorithm(canvas, xCenter + x, yCenter + y, xCenter - x, yCenter + y, color);
    DrawLineBresenhamsAlgorithm(canvas, xCenter + x, yCenter - y, xCenter - x, yCenter - y, color);
    DrawLineBresenhamsAlgorithm(canvas, xCenter + y, yCenter + x, xCenter - y, yCenter + x, color);
    DrawLineBresenhamsAlgorithm(canvas, xCenter + y, yCenter - x, xCenter - y, yCenter - x, color);

    if(d < 0.0){
      d = d + 2.0*x + 1.0;
    }else{
      d = d + 2.0*(x - y) + 1.0;
      y = y - 1.0;
    }
  }
}


void DrawFilledCircleBasicAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color){
  double pixels, a, da, dx, dy;

  /* Place the circle in the center of the pixel. */
  xCenter = floor(xCenter) + 0.5;
  yCenter = floor(yCenter) + 0.5;

  pixels = 2.0*M_PI*radius;

  /* Below a radius of 10 pixels, over-compensate to get a smoother circle. */
  if(radius < 10.0){
    pixels = pixels*10.0;
  }

  da = 2.0*M_PI/pixels;

  /* Draw lines for a half-circle to fill an entire circle. */
  for(a = 0.0; a < M_PI; a = a + da){
    dx = cos(a)*radius;
    dy = sin(a)*radius;

    /* Floor to get the pixel coordinate. */
    DrawVerticalLine1px(canvas, floor(xCenter - dx), floor(yCenter - dy), floor(2.0*dy) + 1.0, color);
  }
}


void DrawTriangle(RGBABitmapImage *canvas, double xCenter, double yCenter, double height, RGBA *color){
  double x1, y1, x2, y2, x3, y3;

  x1 = floor(xCenter + 0.5);
  y1 = floor(floor(yCenter + 0.5) - height);
  x2 = x1 - 2.0*height*tan(M_PI/6.0);
  y2 = floor(y1 + 2.0*height);
  x3 = x1 + 2.0*height*tan(M_PI/6.0);
  y3 = floor(y1 + 2.0*height);

  DrawLine1px(canvas, x1, y1, x2, y2, color);
  DrawLine1px(canvas, x1, y1, x3, y3, color);
  DrawLine1px(canvas, x2, y2, x3, y3, color);
}


void DrawFilledTriangle(RGBABitmapImage *canvas, double xCenter, double yCenter, double height, RGBA *color){
  double i, offset, x1, y1;

  x1 = floor(xCenter + 0.5);
  y1 = floor(floor(yCenter + 0.5) - height);

  for(i = 0.0; i <= 2.0*height; i = i + 1.0){
    offset = floor(i*tan(M_PI/6.0));
    DrawHorizontalLine1px(canvas, x1 - offset, y1 + i, 2.0*offset, color);
  }
}


void DrawLine(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, double thickness, RGBA *color){
  DrawLineBresenhamsAlgorithmThick(canvas, x1, y1, x2, y2, thickness, color);
}


#define sign(a)     \
    a > 0 ? 1  :    \
    a < 0 ? -1 :    \
             0


void DrawLineBresenhamsAlgorithmThick(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, double thickness, RGBA *color){
  double x, y, dx, dy, incX, incY, pdx, pdy, es, el, err, t, r;

  dx = x2 - x1;
  dy = y2 - y1;

  incX = sign(dx);
  incY = sign(dy);

  dx = fabs(dx);
  dy = fabs(dy);

  if(dx > dy){
    pdx = incX;
    pdy = 0.0;
    es = dy;
    el = dx;
  }else{
    pdx = 0.0;
    pdy = incY;
    es = dx;
    el = dy;
  }

  x = x1;
  y = y1;
  err = el/2.0;

  if(thickness >= 3.0){
    r = thickness/2.0;
    DrawCircle(canvas, x, y, r, color);
  }else if(floor(thickness) == 2.0){
    DrawFilledRectangle(canvas, x, y, 2.0, 2.0, color);
  }else if(floor(thickness) == 1.0){
    DrawPixel(canvas, x, y, color);
  }

  for(t = 0.0; t < el; t = t + 1.0){
    err = err - es;
    if(err < 0.0){
      err = err + el;
      x = x + incX;
      y = y + incY;
    }else{
      x = x + pdx;
      y = y + pdy;
    }

    if(thickness >= 3.0){
      r = thickness/2.0;
      DrawCircle(canvas, x, y, r, color);
    }else if(floor(thickness) == 2.0){
      DrawFilledRectangle(canvas, x, y, 2.0, 2.0, color);
    }else if(floor(thickness) == 1.0){
      DrawPixel(canvas, x, y, color);
    }
  }
}


void DrawLineBresenhamsAlgorithm(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, RGBA *color){
  double x, y, dx, dy, incX, incY, pdx, pdy, es, el, err, t;

  dx = x2 - x1;
  dy = y2 - y1;

  incX = sign(dx);
  incY = sign(dy);

  dx = fabs(dx);
  dy = fabs(dy);

  if(dx > dy){
    pdx = incX;
    pdy = 0.0;
    es = dy;
    el = dx;
  }else{
    pdx = 0.0;
    pdy = incY;
    es = dx;
    el = dy;
  }

  x = x1;
  y = y1;
  err = el/2.0;
  DrawPixel(canvas, x, y, color);

  for(t = 0.0; t < el; t = t + 1.0){
    err = err - es;
    if(err < 0.0){
      err = err + el;
      x = x + incX;
      y = y + incY;
    }else{
      x = x + pdx;
      y = y + pdy;
    }

    DrawPixel(canvas, x, y, color);
  }
}


void DrawLineBresenhamsAlgorithmThickPatterned(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, double thickness, bool *pattern, size_t patternLength, NumberReference *offset, RGBA *color){
  double x, y, dx, dy, incX, incY, pdx, pdy, es, el, err, t, r;

  dx = x2 - x1;
  dy = y2 - y1;

  incX = sign(dx);
  incY = sign(dy);

  dx = fabs(dx);
  dy = fabs(dy);

  if(dx > dy){
    pdx = incX;
    pdy = 0.0;
    es = dy;
    el = dx;
  }else{
    pdx = 0.0;
    pdy = incY;
    es = dx;
    el = dy;
  }

  x = x1;
  y = y1;
  err = el/2.0;

  offset->numberValue = fmod(offset->numberValue + 1.0, (double)patternLength*thickness);

  if(pattern[(int)(floor(offset->numberValue/thickness))]){
    if(thickness >= 3.0){
      r = thickness/2.0;
      DrawCircle(canvas, x, y, r, color);
    }else if(floor(thickness) == 2.0){
      DrawFilledRectangle(canvas, x, y, 2.0, 2.0, color);
    }else if(floor(thickness) == 1.0){
      DrawPixel(canvas, x, y, color);
    }
  }

  for(t = 0.0; t < el; t = t + 1.0){
    err = err - es;
    if(err < 0.0){
      err = err + el;
      x = x + incX;
      y = y + incY;
    }else{
      x = x + pdx;
      y = y + pdy;
    }

    offset->numberValue = fmod(offset->numberValue + 1.0, (double)patternLength*thickness);

    if(pattern[(int)(floor(offset->numberValue/thickness))]){
      if(thickness >= 3.0){
        r = thickness/2.0;
        DrawCircle(canvas, x, y, r, color);
      }else if(floor(thickness) == 2.0){
        DrawFilledRectangle(canvas, x, y, 2.0, 2.0, color);
      }else if(floor(thickness) == 1.0){
        DrawPixel(canvas, x, y, color);
      }
    }
  }
}


bool *GetLinePattern5(size_t *returnArrayLength){
  bool *pattern;
  size_t patternLength;

  pattern = (bool*)calloc(sizeof(bool) * (19.0), 1);
  patternLength = 19.0;

  pattern[0] = true;
  pattern[1] = true;
  pattern[2] = true;
  pattern[3] = true;
  pattern[4] = true;
  pattern[5] = true;
  pattern[6] = true;
  pattern[7] = true;
  pattern[8] = true;
  pattern[9] = true;
  pattern[10] = false;
  pattern[11] = false;
  pattern[12] = false;
  pattern[13] = true;
  pattern[14] = true;
  pattern[15] = true;
  pattern[16] = false;
  pattern[17] = false;
  pattern[18] = false;

  *returnArrayLength = patternLength;
  return pattern;
}


bool *GetLinePattern4(size_t *returnArrayLength){
  bool *pattern;
  size_t patternLength;

  pattern = (bool*)calloc(sizeof(bool) * (13.0), 1);
  patternLength = 13.0;

  pattern[0] = true;
  pattern[1] = true;
  pattern[2] = true;
  pattern[3] = true;
  pattern[4] = true;
  pattern[5] = true;
  pattern[6] = true;
  pattern[7] = true;
  pattern[8] = true;
  pattern[9] = true;
  pattern[10] = false;
  pattern[11] = false;
  pattern[12] = false;

  *returnArrayLength = patternLength;
  return pattern;
}


bool *GetLinePattern3(size_t *returnArrayLength){
  bool *pattern;
  size_t patternLength;

  pattern = (bool*)calloc(sizeof(bool) * (13.0), 1);
  patternLength = 13.0;

  pattern[0] = true;
  pattern[1] = true;
  pattern[2] = true;
  pattern[3] = true;
  pattern[4] = true;
  pattern[5] = true;
  pattern[6] = false;
  pattern[7] = false;
  pattern[8] = false;
  pattern[9] = true;
  pattern[10] = true;
  pattern[11] = false;
  pattern[12] = false;

  *returnArrayLength = patternLength;
  return pattern;
}


bool *GetLinePattern2(size_t *returnArrayLength){
  bool *pattern;
  size_t patternLength;

  pattern = (bool*)calloc(sizeof(bool) * (4.0), 1);
  patternLength = 4.0;

  pattern[0] = true;
  pattern[1] = true;
  pattern[2] = false;
  pattern[3] = false;

  *returnArrayLength = patternLength;
  return pattern;
}


bool *GetLinePattern1(size_t *returnArrayLength){
  bool *pattern;
  size_t patternLength;

  pattern = (bool*)calloc(sizeof(bool) * (8.0), 1);
  patternLength = 8.0;

  pattern[0] = true;
  pattern[1] = true;
  pattern[2] = true;
  pattern[3] = true;
  pattern[4] = true;
  pattern[5] = false;
  pattern[6] = false;
  pattern[7] = false;

  *returnArrayLength = patternLength;
  return pattern;
}


char *CreateStringDecimalFromNumber(size_t *returnArrayLength, double decimal){
  StringReference *stringReference = (StringReference *)calloc(sizeof(StringReference), 1);

    stringReference->stringLength = snprintf(NULL, 0, "%.1f", decimal);
    stringReference->string = malloc(sizeof(char) * (stringReference->stringLength + 1));
    snprintf(stringReference->string, stringReference->stringLength + 1, "%.1f", decimal);

  *returnArrayLength = stringReference->stringLength;
  return stringReference->string;
}


char *GetDigitCharacterTable(size_t *returnArrayLength){
  char *numberTable;
  size_t numberTableLength;

  numberTable = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  numberTableLength = strlen(numberTable);

  *returnArrayLength = numberTableLength;
  return numberTable;
}


bool CreateNumberFromDecimalStringWithCheck(char *string, size_t stringLength, NumberReference *decimalReference){
  return CreateNumberFromStringWithCheck(string, stringLength, 10.0, decimalReference);
}


double CreateNumberFromDecimalString(char *string, size_t stringLength){
  NumberReference *doubleReference;
  StringReference *stringReference;
  double number;

  doubleReference = CreateNumberReference(0.0);
  stringReference = CreateStringReference(strparam(""));
  CreateNumberFromStringWithCheck(string, stringLength, 10.0, doubleReference);
  number = doubleReference->numberValue;

  free(doubleReference);
  free(stringReference);

  return number;
}


bool CreateNumberFromStringWithCheck(char *string, size_t stringLength, double base, NumberReference *numberReference){
  bool success;
  BooleanReference *numberIsPositive, *exponentIsPositive;
  NumberArrayReference *beforePoint, *afterPoint, *exponent;

  numberIsPositive = CreateBooleanReference(true);
  exponentIsPositive = CreateBooleanReference(true);
  beforePoint = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);
  afterPoint = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);
  exponent = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);

  if(base >= 2.0 && base <= 36.0){
    success = ExtractPartsFromNumberString(string, stringLength, base, numberIsPositive, beforePoint, afterPoint, exponentIsPositive, exponent);

    if(success){
      numberReference->numberValue = CreateNumberFromParts(base, numberIsPositive->booleanValue, beforePoint->numberArray, beforePoint->numberArrayLength, afterPoint->numberArray, afterPoint->numberArrayLength, exponentIsPositive->booleanValue, exponent->numberArray, exponent->numberArrayLength);
    }
  }else{
    success = false;
    fprintf(stderr, "Base must be from 2 to 36.");
  }

  return success;
}


double CreateNumberFromParts(double base, bool numberIsPositive, double *beforePoint, size_t beforePointLength, double *afterPoint, size_t afterPointLength, bool exponentIsPositive, double *exponent, size_t exponentLength){
  double n, i, p, e;

  n = 0.0;

  for(i = 0.0; i < (double)beforePointLength; i = i + 1.0){
    p = beforePoint[(int)((double)beforePointLength - i - 1.0)];

    n = n + p*pow(base, i);
  }

  for(i = 0.0; i < (double)afterPointLength; i = i + 1.0){
    p = afterPoint[(int)(i)];

    n = n + p*pow(base,  -(i + 1.0));
  }

  if((double)exponentLength > 0.0){
    e = 0.0;
    for(i = 0.0; i < (double)exponentLength; i = i + 1.0){
      p = exponent[(int)((double)exponentLength - i - 1.0)];

      e = e + p*pow(base, i);
    }

    if( !exponentIsPositive ){
      e =  -e;
    }

    n = n*pow(base, e);
  }

  if( !numberIsPositive ){
    n =  -n;
  }

  return n;
}


bool ExtractPartsFromNumberString(char *n, size_t nLength, double base, BooleanReference *numberIsPositive, NumberArrayReference *beforePoint, NumberArrayReference *afterPoint, BooleanReference *exponentIsPositive, NumberArrayReference *exponent){
  double i, j, count;
  bool success, done, complete;

  i = 0.0;
  complete = false;

  if(i < (double)nLength){
    if(n[(int)(i)] == '-'){
      numberIsPositive->booleanValue = false;
      i = i + 1.0;
    }else if(n[(int)(i)] == '+'){
      numberIsPositive->booleanValue = true;
      i = i + 1.0;
    }

    success = true;
  }else{
    success = false;
    fprintf(stderr,  "Number cannot have length zero.");
  }

  if(success){
    done = false;
    count = 0.0;
    for(; i + count < (double)nLength &&  !done ; ){
      if(CharacterIsNumberCharacterInBase(n[(int)(i + count)], base)){
        count = count + 1.0;
      }else{
        done = true;
      }
    }

    if(count >= 1.0){
      beforePoint->numberArray = (double*)calloc(sizeof(double) * (count), 1);
      beforePoint->numberArrayLength = count;

      for(j = 0.0; j < count; j = j + 1.0){
        beforePoint->numberArray[(int)(j)] = GetNumberFromNumberCharacterForBase(n[(int)(i + j)], base);
      }

      i = i + count;

      if(i < (double)nLength){
        success = true;
      }else{
        afterPoint->numberArray = (double*)calloc(sizeof(double) * (0.0), 1);
        afterPoint->numberArrayLength = 0.0;
        exponent->numberArray = (double*)calloc(sizeof(double) * (0.0), 1);
        exponent->numberArrayLength = 0.0;
        success = true;
        complete = true;
      }
    }else{
      success = false;
      fprintf(stderr, "Number must have at least one number after the optional sign.");
    }
  }

  if(success &&  !complete ){
    if(n[(int)(i)] == '.'){
      i = i + 1.0;

      if(i < (double)nLength){
        done = false;
        count = 0.0;
        for(; i + count < (double)nLength &&  !done ; ){
          if(CharacterIsNumberCharacterInBase(n[(int)(i + count)], base)){
            count = count + 1.0;
          }else{
            done = true;
          }
        }

        if(count >= 1.0){
          afterPoint->numberArray = (double*)calloc(sizeof(double) * (count), 1);
          afterPoint->numberArrayLength = count;

          for(j = 0.0; j < count; j = j + 1.0){
            afterPoint->numberArray[(int)(j)] = GetNumberFromNumberCharacterForBase(n[(int)(i + j)], base);
          }

          i = i + count;

          if(i < (double)nLength){
            success = true;
          }else{
            exponent->numberArray = (double*)calloc(sizeof(double) * (0.0), 1);
            exponent->numberArrayLength = 0.0;
            success = true;
            complete = true;
          }
        }else{
          success = false;
          fprintf(stderr, "There must be at least one digit after the decimal point.");
        }
      }else{
        success = false;
        fprintf(stderr, "There must be at least one digit after the decimal point.");
      }
    }else if(base <= 14.0 && (n[(int)(i)] == 'e' || n[(int)(i)] == 'E')){
      if(i < (double)nLength){
        success = true;
        afterPoint->numberArray = (double*)calloc(sizeof(double) * (0.0), 1);
        afterPoint->numberArrayLength = 0.0;
      }else{
        success = false;
        fprintf(stderr, "There must be at least one digit after the exponent.");
      }
    }else{
      success = false;
      fprintf(stderr, "Expected decimal point or exponent symbol.");
    }
  }

  if(success &&  !complete ){
    if(base <= 14.0 && (n[(int)(i)] == 'e' || n[(int)(i)] == 'E')){
      i = i + 1.0;

      if(i < (double)nLength){
        if(n[(int)(i)] == '-'){
          exponentIsPositive->booleanValue = false;
          i = i + 1.0;
        }else if(n[(int)(i)] == '+'){
          exponentIsPositive->booleanValue = true;
          i = i + 1.0;
        }

        if(i < (double)nLength){
          done = false;
          count = 0.0;
          for(; i + count < (double)nLength &&  !done ; ){
            if(CharacterIsNumberCharacterInBase(n[(int)(i + count)], base)){
              count = count + 1.0;
            }else{
              done = true;
            }
          }

          if(count >= 1.0){
            exponent->numberArray = (double*)calloc(sizeof(double) * (count), 1);
            exponent->numberArrayLength = count;

            for(j = 0.0; j < count; j = j + 1.0){
              exponent->numberArray[(int)(j)] = GetNumberFromNumberCharacterForBase(n[(int)(i + j)], base);
            }

            i = i + count;

            if(i == (double)nLength){
              success = true;
            }else{
              success = false;
              fprintf(stderr, "There cannot be any characters past the exponent of the number.");
            }
          }else{
            success = false;
            fprintf(stderr, "There must be at least one digit after the decimal point.");
          }
        }else{
          success = false;
          fprintf(stderr, "There must be at least one digit after the exponent symbol.");
        }
      }else{
        success = false;
        fprintf(stderr, "There must be at least one digit after the exponent symbol.");
      }
    }else{
      success = false;
      fprintf(stderr, "Expected exponent symbol.");
    }
  }

  return success;
}


double GetNumberFromNumberCharacterForBase(char c, double base){
  char *numberTable;
  size_t numberTableLength;
  double i;
  double position;

  numberTable = GetDigitCharacterTable(&numberTableLength);
  position = 0.0;

  for(i = 0.0; i < base; i = i + 1.0){
    if(numberTable[(int)(i)] == c){
      position = i;
    }
  }

  return position;
}


bool CharacterIsNumberCharacterInBase(char c, double base){
  char *numberTable;
  size_t numberTableLength;
  double i;
  bool found;

  numberTable = GetDigitCharacterTable(&numberTableLength);
  found = false;

  for(i = 0.0; i < base; i = i + 1.0){
    if(numberTable[(int)(i)] == c){
      found = true;
    }
  }

  return found;
}


double *StringToNumberArray(size_t *returnArrayLength, char *str, size_t strLength){
  NumberArrayReference *numberArrayReference;
  double *numbers;
  size_t numbersLength;

  numberArrayReference = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);

  StringToNumberArrayWithCheck(str, strLength, numberArrayReference);

  numbers = numberArrayReference->numberArray;
  numbersLength = numberArrayReference->numberArrayLength;

  free(numberArrayReference);

  *returnArrayLength = numbersLength;
  return numbers;
}


bool StringToNumberArrayWithCheck(char *str, size_t strLength, NumberArrayReference *numberArrayReference){
  StringReference **numberStrings;
  size_t numberStringsLength;
  double *numbers;
  size_t numbersLength;
  double i;
  char *numberString, *trimmedNumberString;
  size_t numberStringLength, trimmedNumberStringLength;
  bool success;
  NumberReference *numberReference;

  numberStrings = SplitByString(&numberStringsLength, str, strLength, strparam(","));

  numbers = (double*)calloc(sizeof(double) * ((double)numberStringsLength), 1);
  numbersLength = (double)numberStringsLength;
  success = true;
  numberReference = (NumberReference *)calloc(sizeof(NumberReference), 1);

  for(i = 0.0; i < (double)numberStringsLength; i = i + 1.0){
    numberString = numberStrings[(int)(i)]->string;
    numberStringLength = numberStrings[(int)(i)]->stringLength;
    trimmedNumberString = Trim(&trimmedNumberStringLength, numberString, numberStringLength);
    success = CreateNumberFromDecimalStringWithCheck(trimmedNumberString, trimmedNumberStringLength, numberReference);
    numbers[(int)(i)] = numberReference->numberValue;

    FreeStringReference(numberStrings[(int)(i)]);
    free(trimmedNumberString);
  }

  free(numberStrings);
  free(numberReference);

  numberArrayReference->numberArray = numbers;
  numberArrayReference->numberArrayLength = numbersLength;

  return success;
}


void QuickSortNumbers(double *list, size_t listLength){
  QuickSortNumbersBounds(list, listLength, 0.0, (double)listLength - 1.0);
}


static inline double QuickSortNumbersPartition(double *A, double lo, double hi){
  double pivot, lowPos, j;

  pivot = A[(int)(hi)];
  lowPos = lo;
  for(j = lo; j <= hi - 1.0; j = j + 1.0){
    if(A[(int)(j)] < pivot){
      aSwapElementsOfNumberArray(A, lowPos, j);
      lowPos = lowPos + 1.0;
    }
  }
  aSwapElementsOfNumberArray(A, lowPos, hi);

  return lowPos;
}


void QuickSortNumbersBounds(double *A, size_t ALength, double lo, double hi){
  double p;

  if(lo < hi){
    p = QuickSortNumbersPartition(A, lo, hi);
    QuickSortNumbersBounds(A, ALength, lo, p - 1.0);
    QuickSortNumbersBounds(A, ALength, p + 1.0, hi);
  }
}


double *QuickSortNumbersWithIndexes(size_t *returnArrayLength, double *A, size_t ALength){
  double *indexes;
  size_t indexesLength;
  double i;

  indexes = (double*)calloc(sizeof(double) * ((double)ALength), 1);
  indexesLength = (double)ALength;

  for(i = 0.0; i < (double)ALength; i = i + 1.0){
    indexes[(int)(i)] = i;
  }

  QuickSortNumbersBoundsWithIndexes(A, ALength, indexes, indexesLength, 0.0, (double)ALength - 1.0);

  *returnArrayLength = indexesLength;
  return indexes;
}


static inline double QuickSortNumbersPartitionWithIndexes(double *A, double *indexes, double lo, double hi){
  double pivot, i, j;

  pivot = A[(int)(hi)];
  i = lo - 1.0;
  for(j = lo; j <= hi - 1.0; j = j + 1.0){
    if(A[(int)(j)] < pivot){
      i = i + 1.0;
      aSwapElementsOfNumberArray(A, i, j);
      aSwapElementsOfNumberArray(indexes, i, j);
    }
  }
  aSwapElementsOfNumberArray(A, i + 1.0, hi);
  aSwapElementsOfNumberArray(indexes, i + 1.0, hi);

  return i + 1.0;
}


void QuickSortNumbersBoundsWithIndexes(double *A, size_t ALength, double *indexes, size_t indexesLength, double lo, double hi){
  double p;

  if(lo < hi){
    p = QuickSortNumbersPartitionWithIndexes(A, indexes, lo, hi);
    QuickSortNumbersBoundsWithIndexes(A, ALength, indexes, indexesLength, lo, p - 1.0);
    QuickSortNumbersBoundsWithIndexes(A, ALength, indexes, indexesLength, p + 1.0, hi);
  }
}


#define DivisibleBy(a, b) ((fmod(a, b)) == 0.0)


void aFillNumberArray(double *a, size_t aLength, double value){
  double i;

  for(i = 0.0; i < (double)aLength; i = i + 1.0){
    a[(int)(i)] = value;
  }
}


void aFillString(char *a, size_t aLength, char value){
  double i;

  for(i = 0.0; i < (double)aLength; i = i + 1.0){
    a[(int)(i)] = value;
  }
}


void aFillBooleanArray(bool *a, size_t aLength, bool value){
  double i;

  for(i = 0.0; i < (double)aLength; i = i + 1.0){
    a[(int)(i)] = value;
  }
}


bool aFillNumberArrayRange(double *a, size_t aLength, double value, double from, double to){
  double i, length;
  bool success;

  if(from >= 0.0 && from <= (double)aLength && to >= 0.0 && to <= (double)aLength && from <= to){
    length = to - from;
    for(i = 0.0; i < length; i = i + 1.0){
      a[(int)(from + i)] = value;
    }

    success = true;
  }else{
    success = false;
  }

  return success;
}


bool aFillBooleanArrayRange(bool *a, size_t aLength, bool value, double from, double to){
  double i, length;
  bool success;

  if(from >= 0.0 && from <= (double)aLength && to >= 0.0 && to <= (double)aLength && from <= to){
    length = to - from;
    for(i = 0.0; i < length; i = i + 1.0){
      a[(int)(from + i)] = value;
    }

    success = true;
  }else{
    success = false;
  }

  return success;
}


bool aFillStringRange(char *a, size_t aLength, char value, double from, double to){
  double i, length;
  bool success;

  if(from >= 0.0 && from <= (double)aLength && to >= 0.0 && to <= (double)aLength && from <= to){
    length = to - from;
    for(i = 0.0; i < length; i = i + 1.0){
      a[(int)(from + i)] = value;
    }

    success = true;
  }else{
    success = false;
  }

  return success;
}


double *aCopyNumberArray(size_t *returnArrayLength, double *a, size_t aLength){
  double i;
  double *n;
  size_t nLength;

  n = (double*)calloc(sizeof(double) * ((double)aLength), 1);
  nLength = (double)aLength;

  for(i = 0.0; i < (double)aLength; i = i + 1.0){
    n[(int)(i)] = a[(int)(i)];
  }

  *returnArrayLength = nLength;
  return n;
}


bool *aCopyBooleanArray(size_t *returnArrayLength, bool *a, size_t aLength){
  double i;
  bool *n;
  size_t nLength;

  n = (bool*)calloc(sizeof(bool) * ((double)aLength), 1);
  nLength = (double)aLength;

  for(i = 0.0; i < (double)aLength; i = i + 1.0){
    n[(int)(i)] = a[(int)(i)];
  }

  *returnArrayLength = nLength;
  return n;
}


char *aCopyString(size_t *returnArrayLength, char *a, size_t aLength){
  double i;
  char *n;
  size_t nLength;

  n = (char*)calloc(sizeof(char) * ((double)aLength), 1);
  nLength = (double)aLength;

  for(i = 0.0; i < (double)aLength; i = i + 1.0){
    n[(int)(i)] = a[(int)(i)];
  }

  *returnArrayLength = nLength;
  return n;
}


bool aCopyNumberArrayRange(double *a, size_t aLength, double from, double to, NumberArrayReference *copyReference){
  double i, length;
  double *n;
  size_t nLength;
  bool success;

  if(from >= 0.0 && from <= (double)aLength && to >= 0.0 && to <= (double)aLength && from <= to){
    length = to - from;
    n = (double*)calloc(sizeof(double) * (length), 1);
    nLength = length;

    for(i = 0.0; i < length; i = i + 1.0){
      n[(int)(i)] = a[(int)(from + i)];
    }

    copyReference->numberArray = n;
    copyReference->numberArrayLength = nLength;
    success = true;
  }else{
    success = false;
  }

  return success;
}


bool aCopyBooleanArrayRange(bool *a, size_t aLength, double from, double to, BooleanArrayReference *copyReference){
  double i, length;
  bool *n;
  size_t nLength;
  bool success;

  if(from >= 0.0 && from <= (double)aLength && to >= 0.0 && to <= (double)aLength && from <= to){
    length = to - from;
    n = (bool*)calloc(sizeof(bool) * (length), 1);
    nLength = length;

    for(i = 0.0; i < length; i = i + 1.0){
      n[(int)(i)] = a[(int)(from + i)];
    }

    copyReference->booleanArray = n;
    copyReference->booleanArrayLength = nLength;
    success = true;
  }else{
    success = false;
  }

  return success;
}


bool aCopyStringRange(char *a, size_t aLength, double from, double to, StringReference *copyReference){
  double i, length;
  char *n;
  size_t nLength;
  bool success;

  if(from >= 0.0 && from <= (double)aLength && to >= 0.0 && to <= (double)aLength && from <= to){
    length = to - from;
    n = (char*)calloc(sizeof(char) * (length), 1);
    nLength = length;

    for(i = 0.0; i < length; i = i + 1.0){
      n[(int)(i)] = a[(int)(from + i)];
    }

    copyReference->string = n;
    copyReference->stringLength = nLength;
    success = true;
  }else{
    success = false;
  }

  return success;
}


BooleanReference *CreateBooleanReference(bool value){
  BooleanReference *ref;

  ref = (BooleanReference *)calloc(sizeof(BooleanReference), 1);
  ref->booleanValue = value;

  return ref;
}


BooleanArrayReference *CreateBooleanArrayReferenceLengthValue(double length, bool value){
  BooleanArrayReference *ref;
  double i;

  ref = (BooleanArrayReference *)calloc(sizeof(BooleanArrayReference), 1);
  ref->booleanArray = (bool*)calloc(sizeof(bool) * (length), 1);
  ref->booleanArrayLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    ref->booleanArray[(int)(i)] = value;
  }

  return ref;
}


void freeBooleanArrayReference(BooleanArrayReference *booleanArrayReference){
  free(booleanArrayReference->booleanArray);
  free(booleanArrayReference);
}


CharacterReference *CreateCharacterReference(char value){
  CharacterReference *ref;

  ref = (CharacterReference *)calloc(sizeof(CharacterReference), 1);
  ref->characterValue = value;

  return ref;
}


NumberReference *CreateNumberReference(double value){
  NumberReference *ref;

  ref = (NumberReference *)calloc(sizeof(NumberReference), 1);
  ref->numberValue = value;

  return ref;
}


NumberArrayReference *CreateNumberArrayReference(double *value, size_t valueLength){
  NumberArrayReference *ref;

  ref = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);
  ref->numberArray = value;
  ref->numberArrayLength = valueLength;

  return ref;
}


NumberArrayReference *CreateNumberArrayReferenceLengthValue(double length, double value){
  NumberArrayReference *ref;
  double i;

  ref = (NumberArrayReference *)calloc(sizeof(NumberArrayReference), 1);
  ref->numberArray = (double*)calloc(sizeof(double) * (length), 1);
  ref->numberArrayLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    ref->numberArray[(int)(i)] = value;
  }

  return ref;
}


void freeNumberArrayReference(NumberArrayReference *numberArrayReference){
  free(numberArrayReference->numberArray);
  free(numberArrayReference);
}


StringReference *CreateStringReference(char *value, size_t valueLength){
  StringReference *ref;

  ref = (StringReference *)calloc(sizeof(StringReference), 1);
  ref->string = value;
  ref->stringLength = valueLength;

  return ref;
}


StringReference *CreateStringReferenceLengthValue(double length, char value){
  StringReference *ref;
  double i;

  ref = (StringReference *)calloc(sizeof(StringReference), 1);
  ref->string = (char*)calloc(sizeof(char) * (length), 1);
  ref->stringLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    ref->string[(int)(i)] = value;
  }

  return ref;
}


void FreeStringReference(StringReference *stringReference){
  free(stringReference->string);
  free(stringReference);
}


StringArrayReference *CreateStringArrayReference(StringReference **strings, size_t stringsLength){
  StringArrayReference *ref;

  ref = (StringArrayReference *)calloc(sizeof(StringArrayReference), 1);
  ref->stringArray = strings;
  ref->stringArrayLength = stringsLength;

  return ref;
}


StringArrayReference *CreateStringArrayReferenceLengthValue(double length, char *value, size_t valueLength){
  StringArrayReference *ref;
  double i;

  ref = (StringArrayReference *)calloc(sizeof(StringArrayReference), 1);
  ref->stringArray = (StringReference**)calloc(sizeof(StringReference) * (length), 1);
  ref->stringArrayLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    ref->stringArray[(int)(i)] = CreateStringReference(value, valueLength);
  }

  return ref;
}


void freeStringArrayReference(StringArrayReference *stringArrayReference){
  double i;

  for(i = 0.0; i < (double)stringArrayReference->stringArrayLength; i = i + 1.0){
    free(stringArrayReference->stringArray[(int)(i)]);
  }
  free(stringArrayReference->stringArray);
  free(stringArrayReference);
}


static char pixel_array_base16[] = 
"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe891412108153069c4ffffffffffffffffffffffffffffffffffffffff9409000000000000000049e7ffffffffffffffffffffffffffffffffff61000000000000000000000017ddffffffffffffffffffffffffffffff840000000573d3f5e5a62b00000028f0ffffffffffffffffffffffffffda04000008bcfffffffffff44200000073ffffffffffffffffffffffffff5700000088ffffffffffffffe812000008e3ffffffffffffffffffffffea02000015f9ffffffffffffffff8100000080ffffffffffffffffffffff9c00000072ffffffffffffffffffe40100002fffffffffffffffffffffff51000000b8ffffffffffffffffffff2a000000e2ffffffffffffffffffff21000001f0ffffffffffffffffffff65000000b3fffffffffffffffffff602000018ffffffffffffffffffffff8b0000008affffffffffffffffffd200000036ffffffffffffffffffffffa900000063ffffffffffffffffffc00000004effffffffffffffffffffffc100000052ffffffffffffffffffb500000057ffffffffffffffffffffffc900000046ffffffffffffffffffa90000005fffffffffffffffffffffffd20000003affffffffffffffffffa900000060ffffffffffffffffffffffd30000003affffffffffffffffffb400000057ffffffffffffffffffffffca00000046ffffffffffffffffffc00000004effffffffffffffffffffffc100000052ffffffffffffffffffd100000037ffffffffffffffffffffffa900000063fffffffffffffffffff602000019ffffffffffffffffffffff8b00000089ffffffffffffffffffff21000001f1ffffffffffffffffffff66000000b3ffffffffffffffffffff50000000b8ffffffffffffffffffff2a000000e1ffffffffffffffffffff9c00000073ffffffffffffffffffe40100002fffffffffffffffffffffffea02000015f9ffffffffffffffff8200000080ffffffffffffffffffffffff5700000088ffffffffffffffe812000008e2ffffffffffffffffffffffffda04000008bcfffffffffff44300000073ffffffffffffffffffffffffffff830000000674d3f6e6a72b00000028f0ffffffffffffffffffffffffffffff60000000000000000000000016ddfffffffffffffffffffffffffffffffffe9309000000000000000048e6ffffffffffffffffffffffffffffffffffffffe88f3f1f07132e68c3fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff9d7b28e69441f02000000afffffffffffffffffffffffffffffffffffff6300000000000000000000afffffffffffffffffffffffffffffffffffff6300000000000000000000afffffffffffffffffffffffffffffffffffff6a274c7095b9de64000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000afffffffffffffffffffffffffffffffffffffffffffffffffff67000000affffffffffffffffffffffffffffffffffffff7000000000000000000000000000000003bfffffffffffffffffffffffff7000000000000000000000000000000003bfffffffffffffffffffffffff7000000000000000000000000000000003bffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffd48b56271005142a5ea0f6ffffffffffffffffffffffffffffffffdb7c20000000000000000000001392feffffffffffffffffffffffffffff1f00000000000000000000000000004cf9ffffffffffffffffffffffffff1f0000003784c7e7f9e8b1480000000056ffffffffffffffffffffffffff1f015accffffffffffffffff9701000000b0ffffffffffffffffffffffff58caffffffffffffffffffffff770000003cfffffffffffffffffffffffffffffffffffffffffffffffffff107000002edffffffffffffffffffffffffffffffffffffffffffffffffff3a000000ccffffffffffffffffffffffffffffffffffffffffffffffffff4c000000baffffffffffffffffffffffffffffffffffffffffffffffffff32000000cbffffffffffffffffffffffffffffffffffffffffffffffffec05000002edffffffffffffffffffffffffffffffffffffffffffffffff8d00000039ffffffffffffffffffffffffffffffffffffffffffffffffeb140000009affffffffffffffffffffffffffffffffffffffffffffffff520000002afbffffffffffffffffffffffffffffffffffffffffffffff8c00000003c7ffffffffffffffffffffffffffffffffffffffffffffffb30300000085ffffffffffffffffffffffffffffffffffffffffffffffc50a0000005dfeffffffffffffffffffffffffffffffffffffffffffffd2110000004efbffffffffffffffffffffffffffffffffffffffffffffdb1800000042f8ffffffffffffffffffffffffffffffffffffffffffffe21f00000039f3ffffffffffffffffffffffffffffffffffffffffffffe92600000030efffffffffffffffffffffffffffffffffffffffffffffee2e00000029eafffffffffffffffffffffffffffffffffffffffffffff33700000022e5fffffffffffffffffffffffffffffffffffffffffffff7410000001cdffffffffffffffffffffffffffffffffffffffffffffffb4c00000017d9fffffffffffffffffffffffffffffffffffffffffffffd5900000012d2ffffffffffffffffffffffffffffffffffffffffffffff680000000ecbffffffffffffffffffffffffffffffffffffffffffffffef0000000000000000000000000000000000008bffffffffffffffffffffef0000000000000000000000000000000000008bffffffffffffffffffffef0000000000000000000000000000000000008bffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffe2af8058392817060a1a3f74c8ffffffffffffffffffffffffffffffffeb0000000000000000000000000036cfffffffffffffffffffffffffffffeb000000000000000000000000000004a7ffffffffffffffffffffffffffeb00000f5a9dd0edfbf0ca841900000003c2ffffffffffffffffffffffffec3da8f9fffffffffffffffff0410000002bffffffffffffffffffffffffffffffffffffffffffffffffffee12000000cbffffffffffffffffffffffffffffffffffffffffffffffffff6900000090ffffffffffffffffffffffffffffffffffffffffffffffffff9600000078ffffffffffffffffffffffffffffffffffffffffffffffffff9a0000007effffffffffffffffffffffffffffffffffffffffffffffffff73000000a5fffffffffffffffffffffffffffffffffffffffffffffffff51b000009edfffffffffffffffffffffffffffffffffffffffffffffff7540000007efffffffffffffffffffffffffffffffffffffffffff3d3912400000055fcffffffffffffffffffffffffffffffffff1700000000000000001692feffffffffffffffffffffffffffffffffffff17000000000000002db8feffffffffffffffffffffffffffffffffffffff170000000000000000002bc3fffffffffffffffffffffffffffffffffffffffffffdf0cf922e00000003a5fffffffffffffffffffffffffffffffffffffffffffffffffd8700000007d1ffffffffffffffffffffffffffffffffffffffffffffffffff780000004ffffffffffffffffffffffffffffffffffffffffffffffffffff308000006f6ffffffffffffffffffffffffffffffffffffffffffffffffff3c000000d0ffffffffffffffffffffffffffffffffffffffffffffffffff4d000000c6ffffffffffffffffffffffffffffffffffffffffffffffffff35000000ddffffffffffffffffffffffffffffffffffffffffffffffffea0300000bf9ffffffffffffffffffffffffffffffffffffffffffffffff6200000054ffffffffffffffffffffff47bafefffffffffffffffffff56b00000002cbffffffffffffffffffffff0b001e71a9d7edfbf6e4ba771a000000007cffffffffffffffffffffffff0b0000000000000000000000000000017dffffffffffffffffffffffffff0b000000000000000000000000003cc8ffffffffffffffffffffffffffffe9b989593827160608162a5689dbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffbd0100000000f3fffffffffffffffffffffffffffffffffffffffffffff3200000000000f3ffffffffffffffffffffffffffffffffffffffffffff69000000000000f3ffffffffffffffffffffffffffffffffffffffffffbf01000b0e000000f3fffffffffffffffffffffffffffffffffffffffff42100008e1f000000f3ffffffffffffffffffffffffffffffffffffffff6a000035fc1f000000f3ffffffffffffffffffffffffffffffffffffffc0010004d1ff1f000000f3fffffffffffffffffffffffffffffffffffff42200007affff1f000000f3ffffffffffffffffffffffffffffffffffff6c000026f7ffff1f000000f3ffffffffffffffffffffffffffffffffffc1010001c1ffffff1f000000f3fffffffffffffffffffffffffffffffff523000066ffffffff1f000000f3ffffffffffffffffffffffffffffffff6d000019f0ffffffff1f000000f3ffffffffffffffffffffffffffffffc2010000aeffffffffff1f000000f3fffffffffffffffffffffffffffff524000052ffffffffffff1f000000f3ffffffffffffffffffffffffffff6e00000fe6ffffffffffff1f000000f3ffffffffffffffffffffffffffc30200009affffffffffffff1f000000f3fffffffffffffffffffffffff62400003ffeffffffffffffff1f000000f3ffffffffffffffffffffffff70000008daffffffffffffffff1f000000f3fffffffffffffffffffffff602000086ffffffffffffffffff1f000000f3fffffffffffffffffffffff3000000000000000000000000000000000000000000cbfffffffffffffff3000000000000000000000000000000000000000000cbfffffffffffffff3000000000000000000000000000000000000000000cbffffffffffffffffffffffffffffffffffffffffff1f000000f3ffffffffffffffffffffffffffffffffffffffffffffffffff1f000000f3ffffffffffffffffffffffffffffffffffffffffffffffffff1f000000f3ffffffffffffffffffffffffffffffffffffffffffffffffff1f000000f3ffffffffffffffffffffffffffffffffffffffffffffffffff1f000000f3ffffffffffffffffffffffffffffffffffffffffffffffffff1f000000f3ffffffffffffffffffffffffffffffffffffffffffffffffff1f000000f3ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff4f00000000000000000000000000002fffffffffffffffffffffffffffff4f00000000000000000000000000002fffffffffffffffffffffffffffff4f00000000000000000000000000002fffffffffffffffffffffffffffff4f00000fffffffffffffffffffffffffffffffffffffffffffffffffffff4f00000fffffffffffffffffffffffffffffffffffffffffffffffffffff4f00000fffffffffffffffffffffffffffffffffffffffffffffffffffff4f00000fffffffffffffffffffffffffffffffffffffffffffffffffffff4f00000fffffffffffffffffffffffffffffffffffffffffffffffffffff4f00000fffffffffffffffffffffffffffffffffffffffffffffffffffff4f00000fffffffffffffffffffffffffffffffffffffffffffffffffffff4f000008672f120514275997efffffffffffffffffffffffffffffffffff4f00000000000000000000000b73f6ffffffffffffffffffffffffffffff4f000000000000000000000000002bdeffffffffffffffffffffffffffff60538cbad2e7faf0d599370000000025ebffffffffffffffffffffffffffffffffffffffffffffffffa0090000005bffffffffffffffffffffffffffffffffffffffffffffffffffb100000001d2ffffffffffffffffffffffffffffffffffffffffffffffffff560000007effffffffffffffffffffffffffffffffffffffffffffffffffb80000003dffffffffffffffffffffffffffffffffffffffffffffffffffec00000022fffffffffffffffffffffffffffffffffffffffffffffffffffd00000011ffffffffffffffffffffffffffffffffffffffffffffffffffec00000022ffffffffffffffffffffffffffffffffffffffffffffffffffb80000003cffffffffffffffffffffffffffffffffffffffffffffffffff580000007dffffffffffffffffffffffffffffffffffffffffffffffffb301000000cfffffffffffffffffffffff4cb1fdffffffffffffffffffa40a00000058ffffffffffffffffffffffff17001a6ea9d7eefbf2d69b380000000024e8ffffffffffffffffffffffff1700000000000000000000000000002de0ffffffffffffffffffffffffff17000000000000000000000000127ef9ffffffffffffffffffffffffffffebba8a59372615050a1a3569a6f7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffca753915050d233866a3e0ffffffffffffffffffffffffffffffffffd13f0000000000000000000000f7ffffffffffffffffffffffffffffff9d07000000000000000000000000f7ffffffffffffffffffffffffffff9700000000469fdbf3f5da9e490100f7ffffffffffffffffffffffffffca0300000eb3ffffffffffffffffd84df8fffffffffffffffffffffffffa2d000007c8ffffffffffffffffffffffffffffffffffffffffffffffff9100000081ffffffffffffffffffffffffffffffffffffffffffffffffff28000010f6ffffffffffffffffffffffffffffffffffffffffffffffffc20000006affffffffffffffffffffffffffffffffffffffffffffffffff79000000b2ffffffffffffffffffffffffffffffffffffffffffffffffff43000000ebffeb903d1a0616306fc0ffffffffffffffffffffffffffffff0f000015ffa211000000000000000041dcfffffffffffffffffffffffff30000003087000000000000000000000013c6ffffffffffffffffffffffe30000000f00000055beeef7d8881000000017e6ffffffffffffffffffffd30000000000019dffffffffffffe12200000056ffffffffffffffffffffd100000000006effffffffffffffffce04000002dbffffffffffffffffffdd0000000006eaffffffffffffffffff550000008bffffffffffffffffffe90000000043ffffffffffffffffffffa90000004dfffffffffffffffffff80200000074ffffffffffffffffffffdb0000002cffffffffffffffffffff2200000088ffffffffffffffffffffef00000019ffffffffffffffffffff4d00000088ffffffffffffffffffffee0000001affffffffffffffffffff7e00000074ffffffffffffffffffffdb0000002dffffffffffffffffffffcd00000042ffffffffffffffffffffa900000052ffffffffffffffffffffff21000005e9ffffffffffffffffff5400000093ffffffffffffffffffffff8f0000006dffffffffffffffffcd04000007e6fffffffffffffffffffffff9220000019effffffffffffe1230000006cffffffffffffffffffffffffffc00600000056beeff8d888110000002af3ffffffffffffffffffffffffffffa603000000000000000000000026ddffffffffffffffffffffffffffffffffc8280000000000000000025deffffffffffffffffffffffffffffffffffffffab25a2a1106193b7ed7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff47000000000000000000000000000000000000f7ffffffffffffffffffff47000000000000000000000000000000000003faffffffffffffffffffff4700000000000000000000000000000000004afffffffffffffffffffffffffffffffffffffffffffffffffc1a000000adffffffffffffffffffffffffffffffffffffffffffffffffb300000015faffffffffffffffffffffffffffffffffffffffffffffffff5100000073ffffffffffffffffffffffffffffffffffffffffffffffffea05000000d6ffffffffffffffffffffffffffffffffffffffffffffffff8d00000039ffffffffffffffffffffffffffffffffffffffffffffffffff2c0000009dffffffffffffffffffffffffffffffffffffffffffffffffc90000000cf3ffffffffffffffffffffffffffffffffffffffffffffffff6700000063fffffffffffffffffffffffffffffffffffffffffffffffff60f000000c6ffffffffffffffffffffffffffffffffffffffffffffffffa300000029ffffffffffffffffffffffffffffffffffffffffffffffffff410000008cffffffffffffffffffffffffffffffffffffffffffffffffdf01000005e9ffffffffffffffffffffffffffffffffffffffffffffffff7d00000052fffffffffffffffffffffffffffffffffffffffffffffffffd1e000000b5ffffffffffffffffffffffffffffffffffffffffffffffffb90000001bfcffffffffffffffffffffffffffffffffffffffffffffffff570000007bffffffffffffffffffffffffffffffffffffffffffffffffee07000001ddffffffffffffffffffffffffffffffffffffffffffffffff9300000042ffffffffffffffffffffffffffffffffffffffffffffffffff31000000a5ffffffffffffffffffffffffffffffffffffffffffffffffd000000010f7ffffffffffffffffffffffffffffffffffffffffffffffff6d0000006bfffffffffffffffffffffffffffffffffffffffffffffffff913000000ceffffffffffffffffffffffffffffffffffffffffffffffffa900000031ffffffffffffffffffffffffffffffffffffffffffffffffff4700000094ffffffffffffffffffffffffffffffffffffffffffffffffe302000008eeffffffffffffffffffffffffffffffffffffffffffffffff840000005afffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff9a8602c13050c1d4882dfffffffffffffffffffffffffffffffffffffa918000000000000000000025eeeffffffffffffffffffffffffffffff780000000000000000000000000023e5ffffffffffffffffffffffffff9f0000000037a8e4faf1c66d0500000033fdfffffffffffffffffffffff81600000065fdffffffffffffc40a0000009fffffffffffffffffffffffb600000021faffffffffffffffff8d00000047ffffffffffffffffffffff820000007bffffffffffffffffffeb01000014ffffffffffffffffffffff6d000000a2ffffffffffffffffffff15000001fdffffffffffffffffffff76000000a2ffffffffffffffffffff14000007ffffffffffffffffffffffa10000007bffffffffffffffffffec01000033ffffffffffffffffffffffec08000022fbffffffffffffffff8e00000087ffffffffffffffffffffffff7d00000068fdffffffffffffc70b00001ef2fffffffffffffffffffffffffb5500000039aae5fbf2c87006000013d0fffffffffffffffffffffffffffffe93160000000000000000000153e3ffffffffffffffffffffffffffffffffffbd2e000000000000000780f0ffffffffffffffffffffffffffffffffce3500000000000000000000000e87fcffffffffffffffffffffffffffb3060000004fb2e6faf0cd82150000004ffaffffffffffffffffffffffda0b000004a9ffffffffffffffe93600000076ffffffffffffffffffffff5600000084ffffffffffffffffffe80e000005e2fffffffffffffffffff606000008f4ffffffffffffffffffff6f0000008dffffffffffffffffffcb00000039ffffffffffffffffffffffac0000005cffffffffffffffffffbc0000004affffffffffffffffffffffbe0000004dffffffffffffffffffcc00000039ffffffffffffffffffffffac0000005effffffffffffffffffea00000008f4ffffffffffffffffffff6e0000007cffffffffffffffffffff2f00000085ffffffffffffffffffe70d000000c1ffffffffffffffffffff9300000004a9ffffffffffffffe83400000028fcfffffffffffffffffffffa2d0000000050b2e7fbf2cd821400000002b8ffffffffffffffffffffffffe523000000000000000000000000000299fffffffffffffffffffffffffffff16605000000000000000000002cc5ffffffffffffffffffffffffffffffffffe88e542512040b1b3d72c1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff8a259251008203f8be2ffffffffffffffffffffffffffffffffffffffa91d0000000000000000047ffaffffffffffffffffffffffffffffffff7b00000000000000000000000040f8ffffffffffffffffffffffffffff94000000004db9ecf7da8b1300000057ffffffffffffffffffffffffffdc050000008fffffffffffffe527000000acffffffffffffffffffffffff630000005fffffffffffffffffd406000025fbfffffffffffffffffffffb0c000002e0ffffffffffffffffff5f000000b2ffffffffffffffffffffc600000036ffffffffffffffffffffb50000005fffffffffffffffffffffa000000068ffffffffffffffffffffe700000011feffffffffffffffffff8d0000007cfffffffffffffffffffffb00000000dfffffffffffffffffff8c0000007cfffffffffffffffffffffb00000000b4ffffffffffffffffff9e00000069ffffffffffffffffffffe7000000008dffffffffffffffffffbe00000038ffffffffffffffffffffb6000000007bfffffffffffffffffff606000003e2ffffffffffffffffff62000000006fffffffffffffffffffff4f00000064ffffffffffffffffd8080000000062ffffffffffffffffffffc50000000096ffffffffffffe82b000000000064ffffffffffffffffffffff6c0000000051bbeff8dc8e1500001000000074fffffffffffffffffffffff94f0000000000000000000000288c00000084fffffffffffffffffffffffffd810b000000000000000052ea830000009fffffffffffffffffffffffffffffea8d471d090d2864c1ffff5b000000d4ffffffffffffffffffffffffffffffffffffffffffffffffff2100000dfdffffffffffffffffffffffffffffffffffffffffffffffffd900000052ffffffffffffffffffffffffffffffffffffffffffffffffff75000000b8ffffffffffffffffffffffffffffffffffffffffffffffffe30d000023fefffffffffffffffffffffffffffffffffffffffffffffff945000000b7ffffffffffffffffffffffffff7fa2fdffffffffffffffe8480000005effffffffffffffffffffffffffff63002080c4ecfae7c0740e00000034f4ffffffffffffffffffffffffffff6300000000000000000000000043f0ffffffffffffffffffffffffffffff6300000000000000000000118efdfffffffffffffffffffffffffffffffff4bb7f462b15040b25569ff4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";


char *DigitDataBase16(size_t *returnArrayLength){
    *returnArrayLength = sizeof(pixel_array_base16) / sizeof(*pixel_array_base16);
    return pixel_array_base16;
}


void DrawDigitCharacter(RGBABitmapImage *image, double topx, double topy, double digit){
  double x, y;
  char *allCharData, *colorChars;
  size_t allCharDataLength, colorCharsLength;
  NumberReference *colorReference;
  RGBA *color;

  colorReference = (NumberReference *)calloc(sizeof(NumberReference), 1);
  color = (RGBA *)calloc(sizeof(RGBA), 1);

  colorChars = (char*)calloc(sizeof(char) * (2.0), 1);
  colorCharsLength = 2.0;

  allCharData = DigitDataBase16(&allCharDataLength);

  for(y = 0.0; y < 37.0; y = y + 1.0){
    for(x = 0.0; x < 30.0; x = x + 1.0){
      colorChars[0] = allCharData[(int)(digit*30.0*37.0*2.0 + y*2.0*30.0 + x*2.0 + 0.0)];
      colorChars[1] = allCharData[(int)(digit*30.0*37.0*2.0 + y*2.0*30.0 + x*2.0 + 1.0)];

      ToUpperCase(colorChars, colorCharsLength);
      CreateNumberFromStringWithCheck(colorChars, colorCharsLength, 16.0, colorReference);
      color->r = colorReference->numberValue/255.0;
      color->g = colorReference->numberValue/255.0;
      color->b = colorReference->numberValue/255.0;
      color->a = 1.0;
      SetPixel(image, topx + x, topy + y, color);
    }
  }
}


static char pixel_array[] = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011000000110000000000000000000000110000001100000011000000110000001100000011000000110000000000000000000000000000000000000000000000000000000000000000000000000000110110001101100011011000110110000000000000000000000000001100110011001101111111101100110011001101111111101100110011001100000000000000000000000000000000000011000011111101111111111011000111110000111111000011111000110111111111101111110000110000000000000000000011100001101100011011011011101100000110000011000001100000110111011011011000110110000111000000000000000001111111001100011111100110001101100001110000011100001101100110011001100110011011000011100000000000000000000000000000000000000000000000000000000000000000000000000000110000011100000110000011100000000000000000000001100000001100000001100000011000000110000001100000011000000110000001100000110000011000000000000000000000000110000011000001100000011000000110000001100000011000000110000001100000001100000001100000000000000000000000000000000001001100101011010001111001111111100111100010110101001100100000000000000000000000000000000000000000001100000011000000110001111111111111111000110000001100000011000000000000000000000000000000000000000110000011000001110000011100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111111111111111100000000000000000000000000000000000000000000000000000000000000000001110000011100000000000000000000000000000000000000000000000000000000000000000000000000000001100000011000001100000011000001100000011000001100000011000001100000011000001100000011000000000000000000000000111100011001101100001111000111110011111101101111110011111000111100001101100110001111000000000000000000011111100001100000011000000110000001100000011000000110000001100000011110000111000001100000000000000000001111111100000011000000110000011000001100000110000011000001100000110000001110011101111110000000000000000001111110111001111100000011000000111000000111111011100000110000001100000011100111011111100000000000000000001100000011000000110000001100000011000011111111001100110011011000111100001110000011000000000000000000000111111011100111110000001100000011100000011111110000001100000011000000110000001111111111000000000000000001111110111001111100001111000011111000110111111100000011000000110000001111100111011111100000000000000000000011000000110000001100000011000001100000110000011000001100000011000000110000001111111100000000000000000111111011100111110000111100001111100111011111101110011111000011110000111110011101111110000000000000000001111110111001111100000011000000110000001111111011100111110000111100001111100111011111100000000000000000000000000001110000011100000000000000000000011100000111000000000000000000000000000000000000000000000000000000110000011000001110000011100000000000000000000011100000111000000000000000000000000000000000000000000001100000001100000001100000001100000001100000001100000110000011000001100000110000011000000000000000000000000000000000000011111111111111110000000011111111111111110000000000000000000000000000000000000000000000000000011000001100000110000011000001100000110000000110000000110000000110000000110000000110000000000000000000011000000000000000000000011000000110000011000001100000110000001100001111000011011111100000000000000000111111000000011011110011110110111100101110111011110000110111111000000000000000000000000000000000000000001100001111000011110000111100001111111111110000111100001111000011011001100011110000011000000000000000000001111111111000111100001111000011111000110111111111100011110000111100001111100011011111110000000000000000011111101110011100000011000000110000001100000011000000110000001100000011111001110111111000000000000000000011111101110011111000111100001111000011110000111100001111000011111000110111001100111111000000000000000011111111000000110000001100000011000000110011111100000011000000110000001100000011111111110000000000000000000000110000001100000011000000110000001100000011001111110000001100000011000000111111111100000000000000000111111011100111110000111100001111110011000000110000001100000011000000111110011101111110000000000000000011000011110000111100001111000011110000111111111111000011110000111100001111000011110000110000000000000000011111100001100000011000000110000001100000011000000110000001100000011000000110000111111000000000000000000011111001110111011000110110000001100000011000000110000001100000011000000110000001100000000000000000000011000011011000110011001100011011000011110000011100001111000110110011001101100011110000110000000000000000111111110000001100000011000000110000001100000011000000110000001100000011000000110000001100000000000000001100001111000011110000111100001111000011110000111101101111111111111111111110011111000011000000000000000011100011111000111111001111110011111110111101101111011111110011111100111111000111110001110000000000000000011111101110011111000011110000111100001111000011110000111100001111000011111001110111111000000000000000000000001100000011000000110000001100000011011111111110001111000011110000111110001101111111000000000000000011111100011101101111101111011011110000111100001111000011110000111100001101100110001111000000000000000000110000110110001100110011000110110000111101111111111000111100001111000011111000110111111100000000000000000111111011100111110000001100000011100000011111100000011100000011000000111110011101111110000000000000000000011000000110000001100000011000000110000001100000011000000110000001100000011000111111110000000000000000011111101110011111000011110000111100001111000011110000111100001111000011110000111100001100000000000000000001100000111100001111000110011001100110110000111100001111000011110000111100001111000011000000000000000011000011111001111111111111111111110110111101101111000011110000111100001111000011110000110000000000000000110000110110011001100110001111000011110000011000001111000011110001100110011001101100001100000000000000000001100000011000000110000001100000011000000110000011110000111100011001100110011011000011000000000000000011111111000000110000001100000110000011000111111000110000011000001100000011000000111111110000000000000000001111000000110000001100000011000000110000001100000011000000110000001100000011000011110000000000110000001100000001100000011000000011000000110000000110000001100000001100000011000000011000000110000000000000000000111100001100000011000000110000001100000011000000110000001100000011000000110000001111000000000000000000000000000000000000000000000000000000000000000000000000001100001101100110001111000001100011111111111111110000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110000001110000001100000011100000000000000000111111101100001111000011111111101100000011000011011111100000000000000000000000000000000000000000000000000111111111000011110000111100001111000011011111110000001100000011000000110000001100000011000000000000000001111110110000110000001100000011000000111100001101111110000000000000000000000000000000000000000000000000111111101100001111000011110000111100001111111110110000001100000011000000110000001100000000000000000000001111111000000011000000110111111111000011110000110111111000000000000000000000000000000000000000000000000000001100000011000000110000001100000011000011111100001100000011000000110011001100011110000111111011000011110000001100000011111110110000111100001111000011011111100000000000000000000000000000000000000000000000001100001111000011110000111100001111000011110000110111111100000011000000110000001100000011000000000000000000011000000110000001100000011000000110000001100000011000000000000000000000011000000000000001110000110110001100000011000000110000001100000011000000110000001100000000000000000000001100000000000000000000000000000110001100110011000111110000111100011011001100110110001100000011000000110000001100000011000000000000000001111110000110000001100000011000000110000001100000011000000110000001100000011000000111100000000000000000110110111101101111011011110110111101101111011011011111110000000000000000000000000000000000000000000000000110001101100011011000110110001101100011011000110011111100000000000000000000000000000000000000000000000000111110011000110110001101100011011000110110001100111110000000000000000000000000000000000000001100000011000000110111111111000011110000111100001111000011011111110000000000000000000000000000000011000000110000001100000011111110110000111100001111000011110000111111111000000000000000000000000000000000000000000000000000000011000000110000001100000011000000110000011101111111000000000000000000000000000000000000000000000000011111111100000011000000011111100000001100000011111111100000000000000000000000000000000000000000000000000011100001101100000011000000110000001100000011000011111100001100000011000000110000000000000000000000000001111110011000110110001101100011011000110110001101100011000000000000000000000000000000000000000000000000000110000011110000111100011001100110011011000011110000110000000000000000000000000000000000000000000000001100001111100111111111111101101111000011110000111100001100000000000000000000000000000000000000000000000011000011011001100011110000011000001111000110011011000011000000000000000000000000000000000000001100000110000001100000110000011000001111000110011001100110110000110000000000000000000000000000000000000000000000001111111100000110000011000001100000110000011000001111111100000000000000000000000000000000000000000000000011110000000110000001100000011000000111000000111100011100000110000001100000011000111100000001100000011000000110000001100000011000000110000001100000011000000110000001100000011000000110000001100000000000000000000000111100011000000110000001100000111000111100000011100000011000000110000001100000001111";


char *GetPixelFontData(size_t *returnArrayLength){
    *returnArrayLength = sizeof(pixel_array)/sizeof(*pixel_array);
    return pixel_array;
}


void DrawAsciiCharacter(RGBABitmapImage *image, double topx, double topy, char a, RGBA *color){
  double index, x, y, pixel, basis, ybasis;
  char *allCharData;
  size_t allCharDataLength;

  index = a;
  index = index - 32.0;
  allCharData = GetPixelFontData(&allCharDataLength);

  basis = index * 8.0 * 13.0;

  for(y = 0.0; y < 13.0; y = y + 1.0){
    ybasis = basis + y*8.0;
    for(x = 0.0; x < 8.0; x = x + 1.0){
      pixel = allCharData[(int)(ybasis + x)];
      if(pixel == '1'){
        DrawPixel(image, topx + 8.0 - 1.0 - x, topy + 13.0 - 1.0 - y, color);
      }
    }
  }
}


double GetTextWidth(size_t textLength){
  double charWidth, spacing, width;

  charWidth = 8.0;
  spacing = 2.0;

  if((double)textLength == 0.0){
    width = 0.0;
  } else {
    width = (double)textLength*charWidth + ((double)textLength - 1.0)*spacing;
  }

  return width;
}


ByteArray *ConvertToPNG(RGBABitmapImage *image){
  return ConvertToPNGWithOptions(image, 6.0, false, 0.0, 0.001);
}


ByteArray *ConvertToPNGGrayscale(RGBABitmapImage *image){
  return ConvertToPNGWithOptions(image, 0.0, false, 0.0, 0.001);
}


PHYS *PysicsHeader(double pixelsPerMeter){
  PHYS *phys;

  phys = (PHYS *)calloc(sizeof(PHYS), 1);

  phys->pixelsPerMeter = pixelsPerMeter;

  return phys;
}


ByteArray *ConvertToPNGWithOptions(RGBABitmapImage *image, double colorType, bool setPhys, double pixelsPerMeter, double compressionLevel){
  PNGImage *png;
  ByteArray *pngData, *colorData;

  png = (PNGImage *)calloc(sizeof(PNGImage), 1);

  png->signature = PNGSignature(&png->signatureLength);

  png->ihdr = PNGHeader(image, colorType);

  png->physPresent = setPhys;
  png->phys = PysicsHeader(pixelsPerMeter);

  if(colorType == 6.0){
    colorData = GetPNGColorData(image);
  }else{
    colorData = GetPNGColorDataGreyscale(image);
  }
  png->zlibStruct = ZLibCompressStaticHuffman(colorData, compressionLevel);

  pngData = PNGSerializeChunks(png);

  return pngData;
}


ByteArray *PNGSerializeChunks(PNGImage *png){
  double length, i, chunkLength;
  ByteArray *data;
  NumberReference *position;

  length = (double)png->signatureLength + 12.0 + PNGHeaderLength() + 12.0 + PNGIDATLength(png) + 12.0;
  if(png->physPresent){
    length = length + 4.0 + 4.0 + 1.0 + 12.0;
  }
  data = CreateByteArray(length);
  position = CreateNumberReference(0.0);

  /* Signature */
  for(i = 0.0; i < (double)png->signatureLength; i = i + 1.0){
    WriteByte(data, png->signature[(int)(i)], position);
  }

  /* Header */
  chunkLength = PNGHeaderLength();
  Write4BytesBE(data, chunkLength, position);
  WriteStringBytes(data, strparam("IHDR"), position);
  Write4BytesBE(data, png->ihdr->Width, position);
  Write4BytesBE(data, png->ihdr->Height, position);
  WriteByte(data, png->ihdr->BitDepth, position);
  WriteByte(data, png->ihdr->ColourType, position);
  WriteByte(data, png->ihdr->CompressionMethod, position);
  WriteByte(data, png->ihdr->FilterMethod, position);
  WriteByte(data, png->ihdr->InterlaceMethod, position);
  Write4BytesBE(data, CRC32OfInterval(data, position->numberValue - chunkLength - 4.0, chunkLength + 4.0), position);

  /* pHYs */
  if(png->physPresent){
    chunkLength = 4.0 + 4.0 + 1.0;
    Write4BytesBE(data, chunkLength, position);
    WriteStringBytes(data, strparam("pHYs"), position);

    Write4BytesBE(data, png->phys->pixelsPerMeter, position);
    Write4BytesBE(data, png->phys->pixelsPerMeter, position);
    WriteByte(data, 1.0, position);
    /* 1 = pixels per meter */
    Write4BytesBE(data, CRC32OfInterval(data, position->numberValue - chunkLength - 4.0, chunkLength + 4.0), position);
  }

  /* IDAT */
  chunkLength = PNGIDATLength(png);
  Write4BytesBE(data, chunkLength, position);
  WriteStringBytes(data, strparam("IDAT"), position);
  WriteByte(data, png->zlibStruct->CMF, position);
  WriteByte(data, png->zlibStruct->FLG, position);
  for(i = 0.0; i < ByteArrayLength(png->zlibStruct->CompressedDataBlocks); i = i + 1.0){
    WriteByte(data, GetByte(png->zlibStruct->CompressedDataBlocks, i), position);
  }
  Write4BytesBE(data, png->zlibStruct->Adler32CheckValue, position);
  Write4BytesBE(data, CRC32OfInterval(data, position->numberValue - chunkLength - 4.0, chunkLength + 4.0), position);

  /* IEND */
  chunkLength = 0.0;
  Write4BytesBE(data, chunkLength, position);
  WriteStringBytes(data, strparam("IEND"), position);
  Write4BytesBE(data, CRC32OfInterval(data, position->numberValue - 4.0, 4.0), position);

  return data;
}


double PNGIDATLength(PNGImage *png){
  return 2.0 + ByteArrayLength(png->zlibStruct->CompressedDataBlocks) + 4.0;
}


double PNGHeaderLength(){
  return 4.0 + 4.0 + 1.0 + 1.0 + 1.0 + 1.0 + 1.0;
}


ByteArray *GetPNGColorData(RGBABitmapImage *image){
  ByteArray *colordata;
  size_t pixel;
  double length, x, y, next;
	uint32_t rgba;

  length = 4.0*ImageWidth(image)*ImageHeight(image) + ImageHeight(image);

  colordata = CreateByteArray(length);

  next = 0.0;

  for(y = 0.0; y < ImageHeight(image); y = y + 1.0){
    SetByte(colordata, next, 0.0);
    next = next + 1.0;
    for(x = 0.0; x < ImageWidth(image); x = x + 1.0){
			pixel = x + y * ImageWidth(image);
			rgba = image->pixels[pixel];

      SetByte(colordata, next, (rgba >> 24) & 0xFF);
      next = next + 1.0;
      SetByte(colordata, next, (rgba >> 16) & 0xFF);
      next = next + 1.0;
      SetByte(colordata, next, (rgba >> 8) & 0xFF);
      next = next + 1.0;
      SetByte(colordata, next, (rgba >> 0) & 0xFF);
      next = next + 1.0;
    }
  }

  return colordata;
}


ByteArray *GetPNGColorDataGreyscale(RGBABitmapImage *image){
  ByteArray *colordata;
  double length, x, y, next;
	uint32_t rgba;
	size_t pixel;

  length = ImageWidth(image)*ImageHeight(image) + ImageHeight(image);

  colordata = CreateByteArray(length);

  next = 0.0;

  for(y = 0.0; y < ImageHeight(image); y = y + 1.0){
    SetByte(colordata, next, 0.0);
    next = next + 1.0;
    for(x = 0.0; x < ImageWidth(image); x = x + 1.0){
			pixel = x + y * ImageWidth(image);
			rgba = image->pixels[pixel];

      SetByte(colordata, next, rgba & 0xFF);
      next = next + 1.0;
    }
  }

  return colordata;
}


IHDR *PNGHeader(RGBABitmapImage *image, double colortype){
  IHDR *ihdr;

  ihdr = (IHDR *)calloc(sizeof(IHDR), 1);
  ihdr->Width = ImageWidth(image);
  ihdr->Height = ImageHeight(image);
  /* Truecolour with alpha */
  ihdr->BitDepth = 8.0;
  ihdr->ColourType = colortype;
  ihdr->FilterMethod = 0.0;
  /* None */
  ihdr->CompressionMethod = 0.0;
  /* zlib */
  ihdr->InterlaceMethod = 0.0;
  /* no interlace */
  return ihdr;
}


double *PNGSignature(size_t *returnArrayLength){
  double *s;
  size_t sLength;

  s = (double*)calloc(sizeof(double) * (8.0), 1);
  sLength = 8.0;
  s[0] = 137.0;
  s[1] = 80.0;
  s[2] = 78.0;
  s[3] = 71.0;
  s[4] = 13.0;
  s[5] = 10.0;
  s[6] = 26.0;
  s[7] = 10.0;

  *returnArrayLength = sLength;
  return s;
}


double *PNGReadDataChunks(size_t *returnArrayLength, Chunk **cs, size_t csLength){
  double i, j, length, zlibpos;
  Chunk *c;
  double *zlibData;
  size_t zlibDataLength;

  length = 0.0;
  for(i = 0.0; i < (double)csLength; i = i + 1.0){
    c = cs[(int)(i)];
    if(strncmp(c->type, "IDAT", c->typeLength) == 0){
      length = length + c->length;
    }
  }

  zlibData = (double*)calloc(sizeof(double) * (length), 1);
  zlibDataLength = length;
  zlibpos = 0.0;

  for(i = 0.0; i < (double)csLength; i = i + 1.0){
    c = cs[(int)(i)];
    if(strncmp(c->type,"IDAT", c->typeLength) == 0){
      for(j = 0.0; j < c->length; j = j + 1.0){
        zlibData[(int)(zlibpos)] = GetByte(c->data, j);
        zlibpos = zlibpos + 1.0;
      }
    }
  }

  *returnArrayLength = zlibDataLength;
  return zlibData;
}


bool PNGReadHeader(RGBABitmapImage *image, Chunk **cs, size_t csLength, StringReference *errorMessages){
  double i;
  IHDR *ihdr;
  Chunk *c;
  NumberReference *position;
  RGBABitmapImage *n;
  bool success;

  position = CreateNumberReference(0.0);
  success = false;

  for(i = 0.0; i < (double)csLength; i = i + 1.0){
    c = cs[(int)(i)];
    if(strncmp(c->type, "IHDR", c->typeLength) == 0){
      ihdr = (IHDR *)calloc(sizeof(IHDR), 1);

      ihdr->Width = Read4bytesBE(c->data, position);
      ihdr->Height = Read4bytesBE(c->data, position);
      ihdr->BitDepth = ReadByte(c->data, position);
      ihdr->ColourType = ReadByte(c->data, position);
      ihdr->CompressionMethod = ReadByte(c->data, position);
      ihdr->FilterMethod = ReadByte(c->data, position);
      ihdr->InterlaceMethod = ReadByte(c->data, position);

      n = CreateImage(ihdr->Width, ihdr->Height, GetTransparent());
      image->pixels = n->pixels;
      image->xLength = n->xLength;
      image->yLength = n->yLength;

      if(ihdr->ColourType == 6.0){
        if(ihdr->BitDepth == 8.0){
          if(ihdr->CompressionMethod == 0.0){
            if(ihdr->FilterMethod == 0.0){
              if(ihdr->InterlaceMethod == 0.0){
                success = true;
              }else{
                success = false;
                errorMessages->string = AppendString(&errorMessages->stringLength, errorMessages->string, errorMessages->stringLength, strparam("Interlace method not supported."));
              }
            }else{
              success = false;
              errorMessages->string = AppendString(&errorMessages->stringLength, errorMessages->string, errorMessages->stringLength, strparam("Filter method not supported."));
            }
          }else{
            success = false;
            errorMessages->string = AppendString(&errorMessages->stringLength, errorMessages->string, errorMessages->stringLength, strparam("Compression type not supported."));
          }
        }else{
          success = false;
          errorMessages->string = AppendString(&errorMessages->stringLength, errorMessages->string, errorMessages->stringLength, strparam("Bit depth not supported."));
        }
      }else{
        success = false;
        errorMessages->string = AppendString(&errorMessages->stringLength, errorMessages->string, errorMessages->stringLength, strparam("Color type not supported."));
      }
    }
  }

  return success;
}


Chunk **PNGReadChunks(size_t *returnArrayLength, ByteArray *data, NumberReference *position){
  bool done;
  double prepos;
  double chunks;
  Chunk *c;
  Chunk **cs;
  size_t csLength;
  double i;
  done = false;
  prepos = position->numberValue;
  for(chunks = 0.0;  !done ; chunks = chunks + 1.0){
    c = PNGReadChunk(data, position);
    if(strncmp(c->type, "IEND", c->typeLength) == 0){
      done = true;
    }
  }
  position->numberValue = prepos;
  cs = (Chunk**)calloc(sizeof(Chunk) * (chunks), 1);
  csLength = chunks;
  for(i = 0.0; i < chunks; i = i + 1.0){
    cs[(int)(i)] = PNGReadChunk(data, position);
  }

  *returnArrayLength = csLength;
  return cs;
}


Chunk *PNGReadChunk(ByteArray *data, NumberReference *position){
  Chunk *c;

  c = (Chunk *)calloc(sizeof(Chunk), 1);

  c->length = Read4bytesBE(data, position);
  c->type = (char*)calloc(sizeof(char) * (4.0), 1);
  c->typeLength = 4.0;
  c->type[0] = ReadByte(data, position);
  c->type[1] = ReadByte(data, position);
  c->type[2] = ReadByte(data, position);
  c->type[3] = ReadByte(data, position);
  c->data = ReadXbytes(data, position, c->length);
  c->crc = Read4bytesBE(data, position);

  return c;
}


void WriteStringToStingStream(char *stream, NumberReference *index, char *src, size_t srcLength){
  for(size_t i = 0; i < srcLength; i++)
    stream[(size_t) index->numberValue + i] = src[i];

  index->numberValue = index->numberValue + (double)srcLength;
}


void WriteCharacterToStingStream(char *stream, NumberReference *index, char src){
  stream[(int)(index->numberValue)] = src;
  index->numberValue = index->numberValue + 1.0;
}


void WriteBooleanToStingStream(char *stream, NumberReference *index, bool src){
  if(src){
    WriteStringToStingStream(stream, index, strparam("true"));
  }else{
    WriteStringToStingStream(stream, index, strparam("false"));
  }
}


static char * Substring(size_t *returnArrayLength, char *string, double from, double to){
  char *n;
  size_t nLength;
  double i, length;

  length = to - from;

  n = (char*)calloc(sizeof(char) * (length), 1);
  nLength = length;

  for(i = from; i < to; i = i + 1.0){
    n[(int)(i - from)] = string[(int)(i)];
  }

  *returnArrayLength = nLength;
  return n;
}


bool SubstringWithCheck(char *string, size_t stringLength, double from, double to, StringReference *stringReference){
  bool success;

  if(from >= 0.0 && from <= (double)stringLength && to >= 0.0 && to <= (double)stringLength && from <= to){
    stringReference->string = Substring(&stringReference->stringLength, string, from, to);
    success = true;
  }else{
    success = false;
  }

  return success;
}


char *AppendString(size_t *returnArrayLength, char *s1, size_t s1Length, char *s2, size_t s2Length){
  char *newString;
  size_t newStringLength;

  newString = ConcatenateString(&newStringLength, s1, s1Length, s2, s2Length);

  free(s1);

  *returnArrayLength = newStringLength;
  return newString;
}


char *ConcatenateString(size_t *returnArrayLength, char *s1, size_t s1Length, char *s2, size_t s2Length){
  char *newString;
  size_t newStringLength;
  double i;

  newString = (char*)calloc(sizeof(char) * ((double)s1Length + (double)s2Length), 1);
  newStringLength = (double)s1Length + (double)s2Length;

  for(i = 0.0; i < (double)s1Length; i = i + 1.0){
    newString[(int)(i)] = s1[(int)(i)];
  }

  for(i = 0.0; i < (double)s2Length; i = i + 1.0){
    newString[(int)((double)s1Length + i)] = s2[(int)(i)];
  }

  *returnArrayLength = newStringLength;
  return newString;
}


char *AppendCharacter(size_t *returnArrayLength, char *string, size_t stringLength, char c){
  char *newString;
  size_t newStringLength;

  newString = ConcatenateCharacter(&newStringLength, string, stringLength, c);

  free(string);

  *returnArrayLength = newStringLength;
  return newString;
}


char *ConcatenateCharacter(size_t *returnArrayLength, char *string, size_t stringLength, char c){
  char *newString;
  size_t newStringLength;
  double i;
  newString = (char*)calloc(sizeof(char) * ((double)stringLength + 1.0), 1);
  newStringLength = (double)stringLength + 1.0;

  for(i = 0.0; i < (double)stringLength; i = i + 1.0){
    newString[(int)(i)] = string[(int)(i)];
  }

  newString[(int)((double)stringLength)] = c;

  *returnArrayLength = newStringLength;
  return newString;
}


bool SubstringEqualsWithCheck(char *string, size_t stringLength, double from, char *substring, size_t substringLength, BooleanReference *equalsReference){
  bool success;

  if(from < (double)stringLength){
    success = true;
    equalsReference->booleanValue = SubstringEquals(string, stringLength, from, substring, substringLength);
  }else{
    success = false;
  }

  return success;
}


bool SubstringEquals(char *string, size_t stringLength, double from, char *substring, size_t substringLength){
  double i;
  bool equal;

  equal = true;
  if((double)stringLength - from >= (double)substringLength){
    for(i = 0.0; i < (double)substringLength && equal; i = i + 1.0){
      if(string[(int)(from + i)] != substring[(int)(i)]){
        equal = false;
      }
    }
  }else{
    equal = false;
  }

  return equal;
}


void ToUpperCase(char *string, size_t stringLength){
  double i;

  for(i = 0.0; i < (double)stringLength; i = i + 1.0){
    string[(int)(i)] = toupper(string[(int)(i)]);
  }
}


void ToLowerCase(char *string, size_t stringLength){
  double i;

  for(i = 0.0; i < (double)stringLength; i = i + 1.0){
    string[(int)(i)] = tolower(string[(int)(i)]);
  }
}


char *Trim(size_t *returnArrayLength, char *string, size_t stringLength){
  char *result;
  size_t resultLength;
  double i, lastWhitespaceLocationStart, lastWhitespaceLocationEnd;
  bool firstNonWhitespaceFound;

  /* Find whitepaces at the start. */
  lastWhitespaceLocationStart =  -1.0;
  firstNonWhitespaceFound = false;
  for(i = 0.0; i < (double)stringLength &&  !firstNonWhitespaceFound ; i = i + 1.0){
    if(isspace((int) string[(int)(i)])){
      lastWhitespaceLocationStart = i;
    }else{
      firstNonWhitespaceFound = true;
    }
  }

  /* Find whitepaces at the end. */
  lastWhitespaceLocationEnd = (double)stringLength;
  firstNonWhitespaceFound = false;
  for(i = (double)stringLength - 1.0; i >= 0.0 &&  !firstNonWhitespaceFound ; i = i - 1.0){
    if(isspace((int) string[(int)(i)])){
      lastWhitespaceLocationEnd = i;
    }else{
      firstNonWhitespaceFound = true;
    }
  }

  if(lastWhitespaceLocationStart < lastWhitespaceLocationEnd){
    result = Substring(&resultLength, string, lastWhitespaceLocationStart + 1.0, lastWhitespaceLocationEnd);
  }else{
    result = (char*)calloc(sizeof(char) * (0.0), 1);
    resultLength = 0.0;
  }

  *returnArrayLength = resultLength;
  return result;
}


StringReference **SplitByString(size_t *returnArrayLength, char *toSplit, size_t toSplitLength, char *splitBy, size_t splitByLength){
  StringReference **split;
  size_t splitLength;
  char *next;
  size_t nextLength;
  double i;
  char c;
  StringReference *n;

  split = (StringReference**)calloc(sizeof(StringReference) * (0.0), 1);
  splitLength = 0.0;

  next = (char*)calloc(sizeof(char) * (0.0), 1);
  nextLength = 0.0;
  for(i = 0.0; i < (double)toSplitLength; ){
    c = toSplit[(int)(i)];

    if(SubstringEquals(toSplit, toSplitLength, i, splitBy, splitByLength)){
      n = (StringReference *)calloc(sizeof(StringReference), 1);
      n->string = next;
      n->stringLength = nextLength;
      split = AddString(&splitLength, split, splitLength, n);
      next = (char*)calloc(sizeof(char) * (0.0), 1);
      nextLength = 0.0;
      i = i + (double)splitByLength;
    }else{
      next = AppendCharacter(&nextLength, next, nextLength, c);
      i = i + 1.0;
    }
  }

  n = (StringReference *)calloc(sizeof(StringReference), 1);
  n->string = next;
  n->stringLength = nextLength;
  split = AddString(&splitLength, split, splitLength, n);

  *returnArrayLength = splitLength;
  return split;
}


bool charCharacterIsBefore(char a, char b){
  double ad, bd;

  ad = a;
  bd = b;

  return ad < bd;
}


bool StringIsBefore(char *a, size_t aLength, char *b, size_t bLength){
  bool before, equal, done;
  double i;

  before = false;
  equal = true;
  done = false;

  if((double)aLength == 0.0 && (double)bLength > 0.0){
    before = true;
  }else{
    for(i = 0.0; i < (double)aLength && i < (double)bLength &&  !done ; i = i + 1.0){
      if(a[(int)(i)] != b[(int)(i)]){
        equal = false;
      }
      if(charCharacterIsBefore(a[(int)(i)], b[(int)(i)])){
        before = true;
      }
      if(charCharacterIsBefore(b[(int)(i)], a[(int)(i)])){
        done = true;
      }
    }

    if(equal){
      if((double)aLength < (double)bLength){
        before = true;
      }
    }
  }

  return before;
}


StringReference **AddString(size_t *returnArrayLength, StringReference **list, size_t listLength, StringReference *a){
  StringReference **newlist;
  size_t newlistLength;
  double i;

  newlist = (StringReference**)calloc(sizeof(StringReference) * ((double)listLength + 1.0), 1);
  newlistLength = (double)listLength + 1.0;

  for(i = 0.0; i < (double)listLength; i = i + 1.0){
    newlist[(int)(i)] = list[(int)(i)];
  }
  newlist[(int)((double)listLength)] = a;
		
  free(list);
		
  *returnArrayLength = newlistLength;
  return newlist;
}


void AddStringRef(StringArrayReference *list, StringReference *i){
  list->stringArray = AddString(&list->stringArrayLength, list->stringArray, list->stringArrayLength, i);
}


StringReference **RemoveString(size_t *returnArrayLength, StringReference **list, size_t listLength, double n){
  StringReference **newlist;
  size_t newlistLength;
  double i;

  newlist = (StringReference**)calloc(sizeof(StringReference) * ((double)listLength - 1.0), 1);
  newlistLength = (double)listLength - 1.0;

  if(n >= 0.0 && n < (double)listLength){
    for(i = 0.0; i < (double)listLength; i = i + 1.0){
      if(i < n){
        newlist[(int)(i)] = list[(int)(i)];
      }
      if(i > n){
        newlist[(int)(i - 1.0)] = list[(int)(i)];
      }
    }

    free(list);
  }else{
    free(newlist);
  }
		
  *returnArrayLength = newlistLength;
  return newlist;
}


StringReference *GetStringRef(StringArrayReference *list, double i){
  return list->stringArray[(int)(i)];
}


void RemoveStringRef(StringArrayReference *list, double i){
  list->stringArray = RemoveString(&list->stringArrayLength, list->stringArray, list->stringArrayLength, i);
}


DynamicArrayCharacters *CreateDynamicArrayCharacters(){
  DynamicArrayCharacters *da;

  da = (DynamicArrayCharacters *)calloc(sizeof(DynamicArrayCharacters), 1);
  da->array = (char*)calloc(sizeof(char) * (10.0), 1);
  da->arrayLength = 10.0;
  da->length = 0.0;

  return da;
}


DynamicArrayCharacters *CreateDynamicArrayCharactersWithInitialCapacity(double capacity){
  DynamicArrayCharacters *da;

  da = (DynamicArrayCharacters *)calloc(sizeof(DynamicArrayCharacters), 1);
  da->array = (char*)calloc(sizeof(char) * (capacity), 1);
  da->arrayLength = capacity;
  da->length = 0.0;

  return da;
}


void DynamicArrayAddCharacter(DynamicArrayCharacters *da, char value){
  if(da->length == (double)da->arrayLength){
    DynamicArrayCharactersIncreaseSize(da);
  }

  da->array[(int)(da->length)] = value;
  da->length = da->length + 1.0;
}

void DynamicArrayCharactersIncreaseSize(DynamicArrayCharacters *da){
  double newLength, i;
  char *newArray;
  size_t newArrayLength;

  newLength = round((double)da->arrayLength*3.0/2.0);
  newArray = (char*)calloc(sizeof(char) * (newLength), 1);
  newArrayLength = newLength;

  for(i = 0.0; i < (double)da->arrayLength; i = i + 1.0){
    newArray[(int)(i)] = da->array[(int)(i)];
  }

  free(da->array);

  da->array = newArray;
  da->arrayLength = newArrayLength;
}


bool DynamicArrayCharactersDecreaseSizeNecessary(DynamicArrayCharacters *da){
  bool needsDecrease;

  needsDecrease = false;

  if(da->length > 10.0){
    needsDecrease = da->length <= round((double)da->arrayLength*2.0/3.0);
  }

  return needsDecrease;
}


void DynamicArrayCharactersDecreaseSize(DynamicArrayCharacters *da){
  double newLength, i;
  char *newArray;
  size_t newArrayLength;

  newLength = round((double)da->arrayLength*2.0/3.0);
  newArray = (char*)calloc(sizeof(char) * (newLength), 1);
  newArrayLength = newLength;

  for(i = 0.0; i < newLength; i = i + 1.0){
    newArray[(int)(i)] = da->array[(int)(i)];
  }

  free(da->array);

  da->array = newArray;
  da->arrayLength = newArrayLength;
}


double DynamicArrayCharactersIndex(DynamicArrayCharacters *da, double index){
  return da->array[(int)(index)];
}


double DynamicArrayCharactersLength(DynamicArrayCharacters *da){
  return da->length;
}


void DynamicArrayInsertCharacter(DynamicArrayCharacters *da, double index, char value){
  double i;

  if(da->length == (double)da->arrayLength){
    DynamicArrayCharactersIncreaseSize(da);
  }

  for(i = da->length; i > index; i = i - 1.0){
    da->array[(int)(i)] = da->array[(int)(i - 1.0)];
  }

  da->array[(int)(index)] = value;

  da->length = da->length + 1.0;
}


bool DynamicArrayCharacterSet(DynamicArrayCharacters *da, double index, char value){
  bool success;

  if(index < da->length){
    da->array[(int)(index)] = value;
    success = true;
  }else{
    success = false;
  }

  return success;
}


void DynamicArrayRemoveCharacter(DynamicArrayCharacters *da, double index){
  double i;

  for(i = index; i < da->length - 1.0; i = i + 1.0){
    da->array[(int)(i)] = da->array[(int)(i + 1.0)];
  }

  da->length = da->length - 1.0;

  if(DynamicArrayCharactersDecreaseSizeNecessary(da)){
    DynamicArrayCharactersDecreaseSize(da);
  }
}


void FreeDynamicArrayCharacters(DynamicArrayCharacters *da){
  free(da->array);
  free(da);
}


char *DynamicArrayCharactersToArray(size_t *returnArrayLength, DynamicArrayCharacters *da){
  char *array;
  size_t arrayLength;
  double i;

  array = (char*)calloc(sizeof(char) * (da->length), 1);
  arrayLength = da->length;

  for(i = 0.0; i < da->length; i = i + 1.0){
    array[(int)(i)] = da->array[(int)(i)];
  }

  *returnArrayLength = arrayLength;
  return array;
}


DynamicArrayCharacters *ArrayToDynamicArrayCharactersWithOptimalSize(char *array, size_t arrayLength){
  DynamicArrayCharacters *da;
  double i;
  double c, n, newCapacity;

  c = (double)arrayLength;
  n = (log(c) - 1.0)/log(3.0/2.0);
  newCapacity = floor(n) + 1.0;

  da = CreateDynamicArrayCharactersWithInitialCapacity(newCapacity);

  for(i = 0.0; i < (double)arrayLength; i = i + 1.0){
    da->array[(int)(i)] = array[(int)(i)];
  }

  return da;
}


DynamicArrayCharacters *ArrayToDynamicArrayCharacters(char *array, size_t arrayLength){
  DynamicArrayCharacters *da;

  da = (DynamicArrayCharacters *)calloc(sizeof(DynamicArrayCharacters), 1);
  da->array = aCopyString(&da->arrayLength, array, arrayLength);
  da->length = (double)arrayLength;

  return da;
}


bool DynamicArrayCharactersEqual(DynamicArrayCharacters *a, DynamicArrayCharacters *b){
  bool equal;
  double i;

  equal = true;
  if(a->length == b->length){
    for(i = 0.0; i < a->length && equal; i = i + 1.0){
      if(a->array[(int)(i)] != b->array[(int)(i)]){
        equal = false;
      }
    }
  }else{
    equal = false;
  }

  return equal;
}


LinkedListCharacters *DynamicArrayCharactersToLinkedList(DynamicArrayCharacters *da){
  LinkedListCharacters *ll;
  double i;

  ll = CreateLinkedListCharacter();

  for(i = 0.0; i < da->length; i = i + 1.0){
    LinkedListAddCharacter(ll, da->array[(int)(i)]);
  }

  return ll;
}


DynamicArrayCharacters *LinkedListToDynamicArrayCharacters(LinkedListCharacters *ll){
  DynamicArrayCharacters *da;
  double i;
  LinkedListNodeCharacters *node;

  node = ll->first;

  da = (DynamicArrayCharacters *)calloc(sizeof(DynamicArrayCharacters), 1);
  da->length = LinkedListCharactersLength(ll);

  da->array = (char*)calloc(sizeof(char) * (da->length), 1);
  da->arrayLength = da->length;

  for(i = 0.0; i < da->length; i = i + 1.0){
    da->array[(int)(i)] = node->value;
    node = node->next;
  }

  return da;
}


bool *AddBoolean(size_t *returnArrayLength, bool *list, size_t listLength, bool a){
  bool *newlist;
  size_t newlistLength;
  double i;

  newlist = (bool*)calloc(sizeof(bool) * ((double)listLength + 1.0), 1);
  newlistLength = (double)listLength + 1.0;
  for(i = 0.0; i < (double)listLength; i = i + 1.0){
    newlist[(int)(i)] = list[(int)(i)];
  }
  newlist[(int)((double)listLength)] = a;
		
  free(list);
		
  *returnArrayLength = newlistLength;
  return newlist;
}


void AddBooleanRef(BooleanArrayReference *list, bool i){
  list->booleanArray = AddBoolean(&list->booleanArrayLength, list->booleanArray, list->booleanArrayLength, i);
}


bool *RemoveBoolean(size_t *returnArrayLength, bool *list, size_t listLength, double n){
  bool *newlist;
  size_t newlistLength;
  double i;

  newlist = (bool*)calloc(sizeof(bool) * ((double)listLength - 1.0), 1);
  newlistLength = (double)listLength - 1.0;

  if(n >= 0.0 && n < (double)listLength){
    for(i = 0.0; i < (double)listLength; i = i + 1.0){
      if(i < n){
        newlist[(int)(i)] = list[(int)(i)];
      }
      if(i > n){
        newlist[(int)(i - 1.0)] = list[(int)(i)];
      }
    }

    free(list);
  }else{
    free(newlist);
  }
		
  *returnArrayLength = newlistLength;
  return newlist;
}


bool GetBooleanRef(BooleanArrayReference *list, double i){
  return list->booleanArray[(int)(i)];
}


void RemoveDecimalRef(BooleanArrayReference *list, double i){
  list->booleanArray = RemoveBoolean(&list->booleanArrayLength, list->booleanArray, list->booleanArrayLength, i);
}


LinkedListStrings *CreateLinkedListString(){
  LinkedListStrings *ll;

  ll = (LinkedListStrings *)calloc(sizeof(LinkedListStrings), 1);
  ll->first = (LinkedListNodeStrings *)calloc(sizeof(LinkedListNodeStrings), 1);
  ll->last = ll->first;
  ll->last->end = true;

  return ll;
}


void LinkedListAddString(LinkedListStrings *ll, char *value, size_t valueLength){
  ll->last->end = false;
  ll->last->value = value;
  ll->last->valueLength = valueLength;
  ll->last->next = (LinkedListNodeStrings *)calloc(sizeof(LinkedListNodeStrings), 1);
  ll->last->next->end = true;
  ll->last = ll->last->next;
}


StringReference **LinkedListStringsToArray(size_t *returnArrayLength, LinkedListStrings *ll){
  StringReference **array;
  size_t arrayLength;
  double length, i;
  LinkedListNodeStrings *node;

  node = ll->first;

  length = LinkedListStringsLength(ll);

  array = (StringReference**)calloc(sizeof(StringReference) * (length), 1);
  arrayLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    array[(int)(i)] = (StringReference *)calloc(sizeof(StringReference), 1);
    array[(int)(i)]->string = node->value;
    array[(int)(i)]->stringLength = node->valueLength;
    node = node->next;
  }

  *returnArrayLength = arrayLength;
  return array;
}


double LinkedListStringsLength(LinkedListStrings *ll){
  double l;
  LinkedListNodeStrings *node;

  l = 0.0;
  node = ll->first;
  for(;  !node->end ; ){
    node = node->next;
    l = l + 1.0;
  }

  return l;
}


void freeLinkedListString(LinkedListStrings *ll){
  LinkedListNodeStrings *node, *prev;

  node = ll->first;

  for(;  !node->end ; ){
    prev = node;
    node = node->next;
    free(prev);
  }

  free(node);
}


LinkedListNumbers *CreateLinkedListNumbers(){
  LinkedListNumbers *ll;

  ll = (LinkedListNumbers *)calloc(sizeof(LinkedListNumbers), 1);
  ll->first = (LinkedListNodeNumbers *)calloc(sizeof(LinkedListNodeNumbers), 1);
  ll->last = ll->first;
  ll->last->end = true;

  return ll;
}


LinkedListNumbers **CreateLinkedListNumbersArray(size_t *returnArrayLength, double length){
  LinkedListNumbers **lls;
  size_t llsLength;
  double i;

  lls = (LinkedListNumbers**)calloc(sizeof(LinkedListNumbers) * (length), 1);
  llsLength = length;
  for(i = 0.0; i < (double)llsLength; i = i + 1.0){
    lls[(int)(i)] = CreateLinkedListNumbers();
  }

  *returnArrayLength = llsLength;
  return lls;
}


void LinkedListAddNumber(LinkedListNumbers *ll, double value){
  ll->last->end = false;
  ll->last->value = value;
  ll->last->next = (LinkedListNodeNumbers *)calloc(sizeof(LinkedListNodeNumbers), 1);
  ll->last->next->end = true;
  ll->last = ll->last->next;
}


double LinkedListNumbersLength(LinkedListNumbers *ll){
  double l;
  LinkedListNodeNumbers *node;

  l = 0.0;
  node = ll->first;
  for(;  !node->end ; ){
    node = node->next;
    l = l + 1.0;
  }

  return l;
}


double LinkedListNumbersIndex(LinkedListNumbers *ll, double index){
  double i;
  LinkedListNodeNumbers *node;

  node = ll->first;
  for(i = 0.0; i < index; i = i + 1.0){
    node = node->next;
  }

  return node->value;
}


void LinkedListInsertNumber(LinkedListNumbers *ll, double index, double value){
  double i;
  LinkedListNodeNumbers *node, *tmp;

  if(index == 0.0){
    tmp = ll->first;
    ll->first = (LinkedListNodeNumbers *)calloc(sizeof(LinkedListNodeNumbers), 1);
    ll->first->next = tmp;
    ll->first->value = value;
    ll->first->end = false;
  }else{
    node = ll->first;
    for(i = 0.0; i < index - 1.0; i = i + 1.0){
      node = node->next;
    }

    tmp = node->next;
    node->next = (LinkedListNodeNumbers *)calloc(sizeof(LinkedListNodeNumbers), 1);
    node->next->next = tmp;
    node->next->value = value;
    node->next->end = false;
  }
}


void LinkedListSet(LinkedListNumbers *ll, double index, double value){
  double i;
  LinkedListNodeNumbers *node;

  node = ll->first;
  for(i = 0.0; i < index; i = i + 1.0){
    node = node->next;
  }

  node->next->value = value;
}


void LinkedListRemoveNumber(LinkedListNumbers *ll, double index){
  double i;
  LinkedListNodeNumbers *node, *prev;

  node = ll->first;
  prev = ll->first;

  for(i = 0.0; i < index; i = i + 1.0){
    prev = node;
    node = node->next;
  }

  if(index == 0.0){
    ll->first = prev->next;
  }
  if( !prev->next->end ){
    prev->next = prev->next->next;
  }
}


void freeLinkedListNumbers(LinkedListNumbers *ll){
  LinkedListNodeNumbers *node, *prev;

  node = ll->first;

  for(;  !node->end ; ){
    prev = node;
    node = node->next;
    free(prev);
  }

  free(node);
}


void freeLinkedListNumbersArray(LinkedListNumbers **lls, size_t llsLength){
  double i;

  for(i = 0.0; i < (double)llsLength; i = i + 1.0){
    freeLinkedListNumbers(lls[(int)(i)]);
  }
  free(lls);
}


double *LinkedListNumbersToArray(size_t *returnArrayLength, LinkedListNumbers *ll){
  double *array;
  size_t arrayLength;
  double length, i;
  LinkedListNodeNumbers *node;

  node = ll->first;

  length = LinkedListNumbersLength(ll);

  array = (double*)calloc(sizeof(double) * (length), 1);
  arrayLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    array[(int)(i)] = node->value;
    node = node->next;
  }

  *returnArrayLength = arrayLength;
  return array;
}


LinkedListNumbers *ArrayToLinkedListNumbers(double *array, size_t arrayLength){
  LinkedListNumbers *ll;
  double i;

  ll = CreateLinkedListNumbers();

  for(i = 0.0; i < (double)arrayLength; i = i + 1.0){
    LinkedListAddNumber(ll, array[(int)(i)]);
  }

  return ll;
}


bool LinkedListNumbersEqual(LinkedListNumbers *a, LinkedListNumbers *b){
  bool equal, done;
  LinkedListNodeNumbers *an, *bn;

  an = a->first;
  bn = b->first;

  equal = true;
  done = false;
  for(; equal &&  !done ; ){
    if(an->end == bn->end){
      if(an->end){
        done = true;
      }else if(an->value == bn->value){
        an = an->next;
        bn = bn->next;
      }else{
        equal = false;
      }
    }else{
      equal = false;
    }
  }

  return equal;
}


LinkedListCharacters *CreateLinkedListCharacter(){
  LinkedListCharacters *ll;

  ll = (LinkedListCharacters *)calloc(sizeof(LinkedListCharacters), 1);
  ll->first = (LinkedListNodeCharacters *)calloc(sizeof(LinkedListNodeCharacters), 1);
  ll->last = ll->first;
  ll->last->end = true;

  return ll;
}


void LinkedListAddCharacter(LinkedListCharacters *ll, char value){
  ll->last->end = false;
  ll->last->value = value;
  ll->last->next = (LinkedListNodeCharacters *)calloc(sizeof(LinkedListNodeCharacters), 1);
  ll->last->next->end = true;
  ll->last = ll->last->next;
}


char *LinkedListCharactersToArray(size_t *returnArrayLength, LinkedListCharacters *ll){
  char *array;
  size_t arrayLength;
  double length, i;
  LinkedListNodeCharacters *node;

  node = ll->first;

  length = LinkedListCharactersLength(ll);

  array = (char*)calloc(sizeof(char) * (length), 1);
  arrayLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    array[(int)(i)] = node->value;
    node = node->next;
  }

  *returnArrayLength = arrayLength;
  return array;
}


double LinkedListCharactersLength(LinkedListCharacters *ll){
  double l;
  LinkedListNodeCharacters *node;

  l = 0.0;
  node = ll->first;
  for(;  !node->end ; ){
    node = node->next;
    l = l + 1.0;
  }

  return l;
}


void freeLinkedListCharacter(LinkedListCharacters *ll){
  LinkedListNodeCharacters *node, *prev;

  node = ll->first;

  for(;  !node->end ; ){
    prev = node;
    node = node->next;
    free(prev);
  }

  free(node);
}


void LinkedListCharactersAddString(LinkedListCharacters *ll, char *str, size_t strLength){
  double i;

  for(i = 0.0; i < (double)strLength; i = i + 1.0){
    LinkedListAddCharacter(ll, str[(int)(i)]);
  }
}


DynamicArrayNumbers *CreateDynamicArrayNumbers(){
  DynamicArrayNumbers *da;

  da = (DynamicArrayNumbers *)calloc(sizeof(DynamicArrayNumbers), 1);
  da->array = (double*)calloc(sizeof(double) * (10.0), 1);
  da->arrayLength = 10.0;
  da->length = 0.0;

  return da;
}


DynamicArrayNumbers *CreateDynamicArrayNumbersWithInitialCapacity(double capacity){
  DynamicArrayNumbers *da;

  da = (DynamicArrayNumbers *)calloc(sizeof(DynamicArrayNumbers), 1);
  da->array = (double*)calloc(sizeof(double) * (capacity), 1);
  da->arrayLength = capacity;
  da->length = 0.0;

  return da;
}


void DynamicArrayAddNumber(DynamicArrayNumbers *da, double value){
  if(da->length == (double)da->arrayLength){
    DynamicArrayNumbersIncreaseSize(da);
  }

  da->array[(int)(da->length)] = value;
  da->length = da->length + 1.0;
}


void DynamicArrayNumbersIncreaseSize(DynamicArrayNumbers *da){
  double newLength, i;
  double *newArray;
  size_t newArrayLength;

  newLength = round((double)da->arrayLength*3.0/2.0);
  newArray = (double*)calloc(sizeof(double) * (newLength), 1);
  newArrayLength = newLength;

  for(i = 0.0; i < (double)da->arrayLength; i = i + 1.0){
    newArray[(int)(i)] = da->array[(int)(i)];
  }

  free(da->array);

  da->array = newArray;
  da->arrayLength = newArrayLength;
}


bool DynamicArrayNumbersDecreaseSizeNecessary(DynamicArrayNumbers *da){
  bool needsDecrease;

  needsDecrease = false;

  if(da->length > 10.0){
    needsDecrease = da->length <= round((double)da->arrayLength*2.0/3.0);
  }

  return needsDecrease;
}


void DynamicArrayNumbersDecreaseSize(DynamicArrayNumbers *da){
  double newLength, i;
  double *newArray;
  size_t newArrayLength;

  newLength = round((double)da->arrayLength*2.0/3.0);
  newArray = (double*)calloc(sizeof(double) * (newLength), 1);
  newArrayLength = newLength;

  for(i = 0.0; i < newLength; i = i + 1.0){
    newArray[(int)(i)] = da->array[(int)(i)];
  }

  free(da->array);

  da->array = newArray;
  da->arrayLength = newArrayLength;
}


double DynamicArrayNumbersIndex(DynamicArrayNumbers *da, double index){
  return da->array[(int)(index)];
}


double DynamicArrayNumbersLength(DynamicArrayNumbers *da){
  return da->length;
}


void DynamicArrayInsertNumber(DynamicArrayNumbers *da, double index, double value){
  double i;

  if(da->length == (double)da->arrayLength){
    DynamicArrayNumbersIncreaseSize(da);
  }

  for(i = da->length; i > index; i = i - 1.0){
    da->array[(int)(i)] = da->array[(int)(i - 1.0)];
  }

  da->array[(int)(index)] = value;
  da->length = da->length + 1.0;
}


bool DynamicArrayNumberSet(DynamicArrayNumbers *da, double index, double value){
  bool success;

  if(index < da->length){
    da->array[(int)(index)] = value;
    success = true;
  }else{
    success = false;
  }

  return success;
}


void DynamicArrayRemoveNumber(DynamicArrayNumbers *da, double index){
  double i;

  for(i = index; i < da->length - 1.0; i = i + 1.0){
    da->array[(int)(i)] = da->array[(int)(i + 1.0)];
  }

  da->length = da->length - 1.0;

  if(DynamicArrayNumbersDecreaseSizeNecessary(da)){
    DynamicArrayNumbersDecreaseSize(da);
  }
}


void freeDynamicArrayNumbers(DynamicArrayNumbers *da){
  free(da->array);
  free(da);
}


double *DynamicArrayNumbersToArray(size_t *returnArrayLength, DynamicArrayNumbers *da){
  double *array;
  size_t arrayLength;
  double i;

  array = (double*)calloc(sizeof(double) * (da->length), 1);
  arrayLength = da->length;

  for(i = 0.0; i < da->length; i = i + 1.0){
    array[(int)(i)] = da->array[(int)(i)];
  }

  *returnArrayLength = arrayLength;
  return array;
}


DynamicArrayNumbers *ArrayToDynamicArrayNumbersWithOptimalSize(double *array, size_t arrayLength){
  DynamicArrayNumbers *da;
  double i;
  double c, n, newCapacity;

  c = (double)arrayLength;
  n = (log(c) - 1.0)/log(3.0/2.0);
  newCapacity = floor(n) + 1.0;

  da = CreateDynamicArrayNumbersWithInitialCapacity(newCapacity);

  for(i = 0.0; i < (double)arrayLength; i = i + 1.0){
    da->array[(int)(i)] = array[(int)(i)];
  }

  return da;
}


DynamicArrayNumbers *ArrayToDynamicArrayNumbers(double *array, size_t arrayLength){
  DynamicArrayNumbers *da;

  da = (DynamicArrayNumbers *)calloc(sizeof(DynamicArrayNumbers), 1);
  da->array = aCopyNumberArray(&da->arrayLength, array, arrayLength);
  da->length = (double)arrayLength;

  return da;
}


bool DynamicArrayNumbersEqual(DynamicArrayNumbers *a, DynamicArrayNumbers *b){
  bool equal;
  double i;

  equal = true;
  if(a->length == b->length){
    for(i = 0.0; i < a->length && equal; i = i + 1.0){
      if(a->array[(int)(i)] != b->array[(int)(i)]){
        equal = false;
      }
    }
  }else{
    equal = false;
  }

  return equal;
}


LinkedListNumbers *DynamicArrayNumbersToLinkedList(DynamicArrayNumbers *da){
  LinkedListNumbers *ll;
  double i;

  ll = CreateLinkedListNumbers();

  for(i = 0.0; i < da->length; i = i + 1.0){
    LinkedListAddNumber(ll, da->array[(int)(i)]);
  }

  return ll;
}


DynamicArrayNumbers *LinkedListToDynamicArrayNumbers(LinkedListNumbers *ll){
  DynamicArrayNumbers *da;
  double i;
  LinkedListNodeNumbers *node;

  node = ll->first;

  da = (DynamicArrayNumbers *)calloc(sizeof(DynamicArrayNumbers), 1);
  da->length = LinkedListNumbersLength(ll);

  da->array = (double*)calloc(sizeof(double) * (da->length), 1);
  da->arrayLength = da->length;

  for(i = 0.0; i < da->length; i = i + 1.0){
    da->array[(int)(i)] = node->value;
    node = node->next;
  }

  return da;
}


double DynamicArrayNumbersIndexOf(DynamicArrayNumbers *arr, double n, BooleanReference *foundReference){
  bool found;
  double i;

  found = false;
  for(i = 0.0; i < arr->length &&  !found ; i = i + 1.0){
    if(arr->array[(int)(i)] == n){
      found = true;
    }
  }
  if( !found ){
    i =  -1.0;
  }else{
    i = i - 1.0;
  }

  foundReference->booleanValue = found;

  return i;
}


bool DynamicArrayNumbersIsInArray(DynamicArrayNumbers *arr, double n){
  bool found;
  double i;

  found = false;
  for(i = 0.0; i < arr->length &&  !found ; i = i + 1.0){
    if(arr->array[(int)(i)] == n){
      found = true;
    }
  }

  return found;
}


char *AddCharacter(size_t *returnArrayLength, char *list, size_t listLength, char a){
  char *newlist;
  size_t newlistLength;
  double i;

  newlist = (char*)calloc(sizeof(char) * ((double)listLength + 1.0), 1);
  newlistLength = (double)listLength + 1.0;
  for(i = 0.0; i < (double)listLength; i = i + 1.0){
    newlist[(int)(i)] = list[(int)(i)];
  }
  newlist[(int)((double)listLength)] = a;
		
  free(list);
		
  *returnArrayLength = newlistLength;
  return newlist;
}


void AddCharacterRef(StringReference *list, char i){
  list->string = AddCharacter(&list->stringLength, list->string, list->stringLength, i);
}


char *RemoveCharacter(size_t *returnArrayLength, char *list, size_t listLength, double n){
  char *newlist;
  size_t newlistLength;
  double i;

  newlist = (char*)calloc(sizeof(char) * ((double)listLength - 1.0), 1);
  newlistLength = (double)listLength - 1.0;

  if(n >= 0.0 && n < (double)listLength){
    for(i = 0.0; i < (double)listLength; i = i + 1.0){
      if(i < n){
        newlist[(int)(i)] = list[(int)(i)];
      }
      if(i > n){
        newlist[(int)(i - 1.0)] = list[(int)(i)];
      }
    }

    free(list);
  }else{
    free(newlist);
  }

  *returnArrayLength = newlistLength;
  return newlist;
}


char GetCharacterRef(StringReference *list, double i){
  return list->string[(int)(i)];
}


void RemoveCharacterRef(StringReference *list, double i){
  list->string = RemoveCharacter(&list->stringLength, list->string, list->stringLength, i);
}


ByteArray *ReadXbytes(ByteArray *data, NumberReference *position, double length){
  ByteArray *r;
  double i;

  r = CreateByteArray(length);

  for(i = 0.0; i < length; i = i + 1.0){
    r->bytes[(int)(i)] = ReadByte(data, position);
  }

  return r;
}


double Read4bytesBE(ByteArray *data, NumberReference *position){
  double r;

  r = 0.0;
  r = r + pow(2.0, 24.0)*ReadByte(data, position);
  r = r + pow(2.0, 16.0)*ReadByte(data, position);
  r = r + pow(2.0, 8.0)*ReadByte(data, position);
  r = r + ReadByte(data, position);

  return r;
}


double Read2bytesBE(ByteArray *data, NumberReference *position){
  double r;

  r = 0.0;
  r = r + pow(2.0, 8.0)*ReadByte(data, position);
  r = r + ReadByte(data, position);

  return r;
}


double ReadByte(ByteArray *data, NumberReference *position){
  double next;

  next = data->bytes[(int)(position->numberValue)];
  position->numberValue = position->numberValue + 1.0;

  return next;
}


double Read4bytesLE(ByteArray *data, NumberReference *position){
  double r;

  r = 0.0;
  r = r + ReadByte(data, position);
  r = r + pow(2.0, 8.0)*ReadByte(data, position);
  r = r + pow(2.0, 16.0)*ReadByte(data, position);
  r = r + pow(2.0, 24.0)*ReadByte(data, position);

  return r;
}


void WriteByte(ByteArray *data, double b, NumberReference *position){
  data->bytes[(int)(position->numberValue)] = b;
  position->numberValue = position->numberValue + 1.0;
}


void Write2BytesLE(ByteArray *data, double b, NumberReference *position){
  data->bytes[(int)(position->numberValue)] = BytesRound(fmod(b, pow(2.0, 8.0)));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 8.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
}


void Write4BytesLE(ByteArray *data, double b, NumberReference *position){
  data->bytes[(int)(position->numberValue)] = BytesRound(fmod(b, pow(2.0, 8.0)));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 8.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 16.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 24.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
}


void Write2BytesBE(ByteArray *data, double b, NumberReference *position){
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 8.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = BytesRound(fmod(b, pow(2.0, 8.0)));
  position->numberValue = position->numberValue + 1.0;
}


void Write4BytesBE(ByteArray *data, double b, NumberReference *position){
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 24.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 16.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = fmod(floor(b/pow(2.0, 8.0)), pow(2.0, 8.0));
  position->numberValue = position->numberValue + 1.0;
  data->bytes[(int)(position->numberValue)] = BytesRound(fmod(b, pow(2.0, 8.0)));
  position->numberValue = position->numberValue + 1.0;
}


void WriteStringBytes(ByteArray *data, char *cs, size_t csLength, NumberReference *position){
  double i, v;

  for(i = 0.0; i < (double)csLength; i = i + 1.0){
    v = cs[(int)(i)];
    WriteByte(data, v, position);
  }
}


double BytesRound(double x){
  return floor(x + 0.5);
}


double *ByteArrayToNumberArray(size_t *returnArrayLength, ByteArray *src){
  double *dst;
  size_t dstLength;
  double i, len;

  len = (double)src->bytesLength;
  dst = (double*)calloc(sizeof(double) * (len), 1);
  dstLength = len;
  for(i = 0.0; i < len; i = i + 1.0){
    dst[(int)(i)] = src->bytes[(int)(i)];
  }

  *returnArrayLength = dstLength;
  return dst;
}


ByteArray *NumberArrayToByteArray(double *src, size_t srcLength){
  ByteArray *dst;
  double i, len;

  len = (double)srcLength;
  dst = (ByteArray *)calloc(sizeof(ByteArray), 1);
  dst->bytes = (uint8_t*)calloc(sizeof(uint8_t) * (len), 1);
  dst->bytesLength = len;
  for(i = 0.0; i < len; i = i + 1.0){
    dst->bytes[(int)(i)] = src[(int)(i)];
  }

  return dst;
}


bool ByteArraysEqual(ByteArray *a, ByteArray *b){
  bool equal;
  double i;

  equal = true;
  if((double)a->bytesLength == (double)b->bytesLength){
    for(i = 0.0; i < (double)a->bytesLength && equal; i = i + 1.0){
      if(a->bytes[(int)(i)] != b->bytes[(int)(i)]){
        equal = false;
      }
    }
  }else{
    equal = false;
  }

  return equal;
}


ByteArray *CopyByteArray(ByteArray *a){
  double i;
  ByteArray *n;

  n = (ByteArray *)calloc(sizeof(ByteArray), 1);
  n->bytes = (uint8_t*)calloc(sizeof(uint8_t) * ((double)a->bytesLength), 1);
  n->bytesLength = (double)a->bytesLength;

  for(i = 0.0; i < (double)a->bytesLength; i = i + 1.0){
    n->bytes[(int)(i)] = a->bytes[(int)(i)];
  }

  return n;
}


double ByteArrayLength(ByteArray *response){
  return (double)response->bytesLength;
}


ByteArray *CreateAndFillByteArray(double length, double value){
  ByteArray *bytes;
  double i;

  bytes = (ByteArray *)calloc(sizeof(ByteArray), 1);
  bytes->bytes = (uint8_t*)calloc(sizeof(uint8_t) * (length), 1);
  bytes->bytesLength = length;

  for(i = 0.0; i < length; i = i + 1.0){
    bytes->bytes[(int)(i)] = value;
  }

  return bytes;
}


ByteArray *CreateByteArray(double length){
  ByteArray *bytes;

  bytes = (ByteArray *)calloc(sizeof(ByteArray), 1);
  bytes->bytes = (uint8_t*)calloc(sizeof(uint8_t) * (length), 1);
  bytes->bytesLength = length;

  return bytes;
}


void SetByte(ByteArray *array, double index, double value){
  array->bytes[(int)(index)] = value;
}


double GetByte(ByteArray *array, double index){
  return array->bytes[(int)(index)];
}


void AssertByteArraysEqual(ByteArray *a, ByteArray *b, NumberReference *failures){
  double i;

  if((double)a->bytesLength == (double)b->bytesLength){
    for(i = 0.0; i < (double)a->bytesLength; i = i + 1.0){
      if(GetByte(a, i) != GetByte(b, i)){
        failures->numberValue = failures->numberValue + 1.0;
      }
    }
  }else{
    failures->numberValue = failures->numberValue + 1.0;
  }
}


void freeByteArray(ByteArray *byteArray){
  free(byteArray->bytes);
  free(byteArray);
}


bool CopyByteArrayRange(ByteArray *a, double from, double to, ByteArray *b){
  double i, length;
  uint8_t *n;
  size_t nLength;
  bool success;

  if(from >= 0.0 && from <= (double)a->bytesLength && to >= 0.0 && to <= (double)a->bytesLength && from <= to){
    length = to - from;
    n = (uint8_t*)calloc(sizeof(uint8_t) * (length), 1);
    nLength = length;

    for(i = 0.0; i < length; i = i + 1.0){
      n[(int)(i)] = a->bytes[(int)(from + i)];
    }

    b->bytes = n;
    b->bytesLength = nLength;
    success = true;
  }else{
    success = false;
  }

  return success;
}


static void FreeBase64ByteCombinations(StringReference **comb){
  for(double i = 0.0; i < 256.0; i = i + 1.0)
    FreeStringReference(comb[(int)(i)]);

  free(comb);
}


char *BytesToTextBase16(size_t *returnArrayLength, double *bytes, size_t bytesLength){
  StringReference **comb;
  size_t combLength;
  double i, v;
  char *r;
  size_t rLength;

  r = (char*)calloc(sizeof(char) * ((double)bytesLength*2.0), 1);
  rLength = (double)bytesLength*2.0;
  comb = GenerateBase16ByteCombinations(&combLength);

  for(i = 0.0; i < (double)bytesLength; i = i + 1.0){
    v = bytes[(int)(i)];
    r[(int)(2.0*i + 0.0)] = comb[(int)(v)]->string[0];
    r[(int)(2.0*i + 1.0)] = comb[(int)(v)]->string[1];
  }

  FreeBase64ByteCombinations(comb);

  *returnArrayLength = rLength;
  return r;
}


double *TextToBytesBase16(size_t *returnArrayLength, char *string, size_t stringLength){
  StringReference **comb;
  size_t combLength;
  double i, n, c, d;
  double *r;
  size_t rLength;
  double *reverseBase16;
  char *codes;
  size_t codesLength;

  codes = "0123456789ABCDEF";
  codesLength = strlen(codes);

  reverseBase16 = (double*)calloc(sizeof(double) * (128.0), 1);
  for(i = 0.0; i < (double)codesLength; i = i + 1.0){
    n = codes[(int)(i)];
    reverseBase16[(int)(n)] = i;
  }

  if(fmod((double)stringLength, 2.0) == 0.0){
    r = (double*)calloc(sizeof(double) * ((double)stringLength/2.0), 1);
    rLength = (double)stringLength/2.0;
    comb = GenerateBase16ByteCombinations(&combLength);

    for(i = 0.0; i < (double)stringLength/2.0; i = i + 1.0){
      c = string[(int)(2.0*i + 0.0)];
      d = string[(int)(2.0*i + 1.0)];

      n = reverseBase16[(int)(c)]*16.0 + reverseBase16[(int)(d)];

      r[(int)(i)] = n;
    }

    FreeBase64ByteCombinations(comb);
  }else{
    r = (double*)calloc(sizeof(double) * (0.0), 1);
    rLength = 0.0;
  }

  free(reverseBase16);

  *returnArrayLength = rLength;
  return r;
}


StringReference **GenerateBase16ByteCombinations(size_t *returnArrayLength){
  StringReference **comb;
  size_t combLength;
  char *codes;
  double i, j;
  comb = (StringReference**)calloc(sizeof(StringReference) * (256.0), 1);
  combLength = 256.0;

  codes = "0123456789ABCDEF";

  for(i = 0.0; i < 256.0; i = i + 1.0){
    comb[(int)(i)] = (StringReference *)calloc(sizeof(StringReference), 1);
    comb[(int)(i)]->string = (char*)calloc(sizeof(char) * (2.0), 1);
    comb[(int)(i)]->stringLength = 2.0;
  }

  for(j = 0.0; j < 16.0; j = j + 1.0){
    for(i = 0.0; i < 16.0; i = i + 1.0){
      comb[(int)(j*16.0 + i)]->string[1] = codes[(int)(i)];
      comb[(int)(j*16.0 + i)]->string[0] = codes[(int)(j)];
    }
  }

  *returnArrayLength = combLength;
  return comb;
}


double *MakeCRC32Table(size_t *returnArrayLength){
  double c, n, k;
  double *crcTable;
  size_t crcTableLength;

  crcTable = (double*)calloc(sizeof(double) * (256.0), 1);
  crcTableLength = 256.0;

  for(n = 0.0; n < 256.0; n = n + 1.0){
    c = n;
    for(k = 0.0; k < 8.0; k = k + 1.0){
      if( !DivisibleBy(c, 2.0) ){
        c = Xor4Byte(3988292384.0, floor(c/2.0));
      }else{
        c = floor(c/2.0);
      }
    }
    crcTable[(int)(n)] = c;
  }

  *returnArrayLength = crcTableLength;
  return crcTable;
}


static inline double UpdateCRC32(double crc, ByteArray *buf, double *crc_table){
  double n, index;

  for(n = 0.0; n < ByteArrayLength(buf); n = n + 1.0){
    index = And4Byte(Xor4Byte(crc, GetByte(buf, n)), pow(2.0, 8.0) - 1.0);
    crc = Xor4Byte(crc_table[(int)(index)], floor(crc/pow(2.0, 8.0)));
  }

  return crc;
}


double CalculateCRC32(ByteArray *buf){
  double *crcTable;
  size_t crcTableLength;
  double b32max, value;

  crcTable = MakeCRC32Table(&crcTableLength);

  b32max = pow(2.0, 32.0) - 1.0;
  value = UpdateCRC32(b32max, buf, crcTable);

  return Xor4Byte(value, b32max);
}


double CRC32OfInterval(ByteArray *data, double from, double length){
  ByteArray *crcBase;
  double i, crc;

  crcBase = CreateByteArray(length);

  for(i = 0.0; i < length; i = i + 1.0){
    SetByte(crcBase, i, GetByte(data, from + i));
  }

  crc = CalculateCRC32(crcBase);

  free(crcBase);

  return crc;
}


static double ComputeAdler32(ByteArray *data){
  double a, b, m, i;

  a = 1.0;
  b = 0.0;
  m = 65521.0;

  for(i = 0.0; i < ByteArrayLength(data); i = i + 1.0){
    a = fmod(a + GetByte(data, i), m);
    b = fmod(b + a, m);
  }

  return b*pow(2.0, 16.0) + a;
}


ZLIBStruct *ZLibCompressNoCompression(ByteArray *data){
  ZLIBStruct *zlibStruct;

  zlibStruct = (ZLIBStruct *)calloc(sizeof(ZLIBStruct), 1);

  zlibStruct->CMF = 120.0;
  zlibStruct->FLG = 1.0;
  zlibStruct->CompressedDataBlocks = DeflateDataNoCompression(data);
  zlibStruct->Adler32CheckValue = ComputeAdler32(data);

  return zlibStruct;
}


ZLIBStruct *ZLibCompressStaticHuffman(ByteArray *data, double level){
  ZLIBStruct *zlibStruct;

  zlibStruct = (ZLIBStruct *)calloc(sizeof(ZLIBStruct), 1);

  zlibStruct->CMF = 120.0;
  zlibStruct->FLG = 1.0;
  zlibStruct->CompressedDataBlocks = DeflateDataStaticHuffman(data, level);
  zlibStruct->Adler32CheckValue = ComputeAdler32(data);

  return zlibStruct;
}


double And4Byte(double n1, double n2){
    if((double)n1 >= 0.0 && (double)n1 <= (double)0xFFFFFFFFUL && (double)n2 >= 0.0 && (double)n2 <= (double)0xFFFFFFFFUL){
      return (unsigned long)n1 & (unsigned long)n2;
    }else{
      return 0.0;
    }
}


double Or4Byte(double n1, double n2){
    if((double)n1 >= 0.0 && (double)n1 <= (double)0xFFFFFFFFUL && (double)n2 >= 0.0 && (double)n2 <= (double)0xFFFFFFFFUL){
      return (unsigned long)n1 | (unsigned long)n2;
    }else{
      return 0.0;
    }
}


double Xor4Byte(double n1, double n2){
    if((double)n1 >= 0.0 && (double)n1 <= (double)0xFFFFFFFFUL && (double)n2 >= 0.0 && (double)n2 <= (double)0xFFFFFFFFUL){
      return (unsigned long)n1 ^ (unsigned long)n2;
    }else{
      return 0.0;
    }
}


double Not2Byte(double b){
    if((double)b >= 0.0 && (double)b <= (double)0xFFFFUL){
      return ~(unsigned long)b & 0xFFFFUL;
    }else{
      return 0.0;
    }
}


double ShiftLeft4Byte(double b, double amount){
    if((double)b >= 0.0 && (double)b <= (double)0xFFFFFFFF && (double)amount >= 0.0 && (double)amount < (double)32){
      return (unsigned long)b << (unsigned long)amount;
    }else{
      return 0.0;
    }
}


double ShiftRight4Byte(double b, double amount){
    if((double)b >= 0.0 && (double)b <= (double)0xFFFFFFFF && (double)amount >= 0.0 && (double)amount < (double)32){
      return (unsigned long)b >> (unsigned long)amount;
    }else{
      return 0.0;
    }
}


ByteArray *Pack(ByteArray *data, double level){
  ByteArray *output, *trimmed;
  double read, i, j, write, start, len, beststart, bestlen, maxforward, size, plain, maxlength;
  bool done;

  /* 5-65536 */
  size = level - 1.0;
  maxlength = 255.0;

  output = CreateByteArray(ceil(ByteArrayLength(data)*1.1));

  write = 0.0;
  read = 0.0;
  plain = 0.0;
  for(; read < ByteArrayLength(data); ){
    /* Search */
    bestlen = 0.0;
    beststart = 0.0;
    /* Lower limit at most 0. */
    start = fmax(0.0, read - size);
    for(i = read - 1.0; i >= start; i = i - 1.0){
      /* Maximum forward is limited by 1) size, 2) data length or 3) distance to down to search start (i) */
      maxforward = fmin(maxlength, fmin(ByteArrayLength(data) - read, read - i));
      len = 0.0;
      done = false;
      for(j = 0.0; j < maxforward &&  !done ; j = j + 1.0){
        if(GetByte(data, i + j) == GetByte(data, read + j)){
          len = len + 1.0;
        }else{
          done = true;
        }
      }

      if(len > bestlen){
        bestlen = len;
        beststart = read - i;
      }
    }

    /* Write data. */
    if(read + 1.0 == ByteArrayLength(data)){
      plain = plain + 1.0;
      read = read + 1.0;
    }

    if(bestlen >= 4.0 || plain == 254.0 || read == ByteArrayLength(data)){
      if(plain > 0.0){
        SetByte(output, write, plain);
        write = write + 1.0;
        for(i = 0.0; i < plain; i = i + 1.0){
          SetByte(output, write, GetByte(data, read - plain + i));
          write = write + 1.0;
        }
      }
      plain = 0.0;
    }

    if(bestlen >= 4.0){
      if(beststart > 255.0){
        SetByte(output, write + 0.0, 0.0);
        SetByte(output, write + 1.0, floor(beststart/256.0));
        SetByte(output, write + 2.0, fmod(beststart, 256.0));
        SetByte(output, write + 3.0, bestlen);
        write = write + 4.0;
      }else{
        SetByte(output, write + 0.0, 255.0);
        SetByte(output, write + 1.0, beststart);
        SetByte(output, write + 2.0, bestlen);
        write = write + 3.0;
      }

      read = read + bestlen;
    }else{
      plain = plain + 1.0;
      read = read + 1.0;
    }
  }

  trimmed = CreateByteArray(write);

  for(i = 0.0; i < write; i = i + 1.0){
    SetByte(trimmed, i, GetByte(output, i));
  }

  return trimmed;
}


ByteArray *Unpack(ByteArray *data){
  ByteArray *output;
  double read, write, i, dist, length, pos, start, unpackedLength;

  /* Compute length */
  unpackedLength = 0.0;
  read = 0.0;
  for(; read < ByteArrayLength(data); ){
    length = GetByte(data, read);

    if(length == 0.0 || length == 255.0){
      if(length == 0.0){
        unpackedLength = unpackedLength + GetByte(data, read + 3.0);
        read = read + 4.0;
      }else{
        unpackedLength = unpackedLength + GetByte(data, read + 2.0);
        read = read + 3.0;
      }
    }else{
      read = read + 1.0;
      unpackedLength = unpackedLength + length;
      read = read + length;
    }
  }

  /* calloc result memory. */
  output = CreateByteArray(unpackedLength);

  /* Unpack */
  write = 0.0;
  read = 0.0;
  for(; read < ByteArrayLength(data); ){
    length = GetByte(data, read);

    if(length == 0.0 || length == 255.0){
      if(length == 0.0){
        dist = GetByte(data, read + 1.0)*256.0;
        dist = dist + GetByte(data, read + 2.0);
        length = GetByte(data, read + 3.0);
        read = read + 4.0;
      }else{
        dist = GetByte(data, read + 1.0);
        length = GetByte(data, read + 2.0);
        read = read + 3.0;
      }

      start = write - dist;

      for(i = 0.0; i < length; i = i + 1.0){
        pos = start + i;
        SetByte(output, write, GetByte(output, pos));
        write = write + 1.0;
      }
    }else{
      read = read + 1.0;
      for(i = 0.0; i < length; i = i + 1.0){
        SetByte(output, write, GetByte(data, read));
        write = write + 1.0;
        read = read + 1.0;
      }
    }
  }

  return output;
}


static void GetDeflateStaticHuffmanCode(double b, NumberReference *code, NumberReference *length, double *bitReverseLookupTable){
  double reversed;

  if(b >= 0.0 && b <= 143.0){
    code->numberValue = 48.0 + b;
    length->numberValue = 8.0;
  }else if(b >= 144.0 && b <= 255.0){
    code->numberValue = b - 144.0 + 400.0;
    length->numberValue = 9.0;
  }else if(b >= 256.0 && b <= 279.0){
    code->numberValue = b - 256.0 + 0.0;
    length->numberValue = 7.0;
  }else if(b >= 280.0 && b <= 287.0){
    code->numberValue = b - 280.0 + 192.0;
    length->numberValue = 8.0;
  }

  reversed = bitReverseLookupTable[(int)(code->numberValue)];
  code->numberValue = ShiftRight4Byte(reversed, 32.0 - length->numberValue);
}


static void GetDeflateDistanceCode(
        double distance
        , NumberReference *code
        , NumberReference *distanceAdditionReference
        , NumberReference *distanceAdditionLengthReference
        , double *bitReverseLookupTable){
  double reversed;

  if(distance >= 1.0 && distance <= 4.0){
    code->numberValue = distance - 1.0;
    distanceAdditionLengthReference->numberValue = 0.0;
  }else if(distance >= 5.0 && distance <= 8.0){
    code->numberValue = 4.0 + floor((distance - 5.0)/2.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 5.0, 2.0));
    distanceAdditionLengthReference->numberValue = 1.0;
  }else if(distance >= 9.0 && distance <= 16.0){
    code->numberValue = 6.0 + floor((distance - 9.0)/4.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 9.0, 4.0));
    distanceAdditionLengthReference->numberValue = 2.0;
  }else if(distance >= 17.0 && distance <= 32.0){
    code->numberValue = 8.0 + floor((distance - 17.0)/8.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 17.0, 8.0));
    distanceAdditionLengthReference->numberValue = 3.0;
  }else if(distance >= 33.0 && distance <= 64.0){
    code->numberValue = 10.0 + floor((distance - 33.0)/16.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 33.0, 16.0));
    distanceAdditionLengthReference->numberValue = 4.0;
  }else if(distance >= 65.0 && distance <= 128.0){
    code->numberValue = 12.0 + floor((distance - 65.0)/32.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 65.0, 32.0));
    distanceAdditionLengthReference->numberValue = 5.0;
  }else if(distance >= 129.0 && distance <= 256.0){
    code->numberValue = 14.0 + floor((distance - 129.0)/64.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 129.0, 64.0));
    distanceAdditionLengthReference->numberValue = 6.0;
  }else if(distance >= 257.0 && distance <= 512.0){
    code->numberValue = 16.0 + floor((distance - 257.0)/128.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 257.0, 128.0));
    distanceAdditionLengthReference->numberValue = 7.0;
  }else if(distance >= 513.0 && distance <= 1024.0){
    code->numberValue = 18.0 + floor((distance - 513.0)/256.0);
    distanceAdditionReference->numberValue = floor(fmod(distance - 513.0, 256.0));
    distanceAdditionLengthReference->numberValue = 8.0;
  }else if(distance >= 1025.0 && distance <= 2048.0){
    code->numberValue = 20.0 + floor((distance - 1025.0)/pow(2.0, 9.0));
    distanceAdditionReference->numberValue = floor(fmod(distance - 1025.0, pow(2.0, 9.0)));
    distanceAdditionLengthReference->numberValue = 9.0;
  }else if(distance >= 2049.0 && distance <= 4096.0){
    code->numberValue = 22.0 + floor((distance - 2049.0)/pow(2.0, 10.0));
    distanceAdditionReference->numberValue = floor(fmod(distance - 2049.0, pow(2.0, 10.0)));
    distanceAdditionLengthReference->numberValue = 10.0;
  }else if(distance >= 4097.0 && distance <= 8192.0){
    code->numberValue = 24.0 + floor((distance - 4097.0)/pow(2.0, 11.0));
    distanceAdditionReference->numberValue = floor(fmod(distance - 4097.0, pow(2.0, 11.0)));
    distanceAdditionLengthReference->numberValue = 11.0;
  }else if(distance >= 8193.0 && distance <= 16384.0){
    code->numberValue = 26.0 + floor((distance - 8193.0)/pow(2.0, 12.0));
    distanceAdditionReference->numberValue = floor(fmod(distance - 8193.0, pow(2.0, 12.0)));
    distanceAdditionLengthReference->numberValue = 12.0;
  }else if(distance >= 16385.0 && distance <= 32768.0){
    code->numberValue = 28.0 + floor((distance - 16385.0)/pow(2.0, 13.0));
    distanceAdditionReference->numberValue = floor(fmod(distance - 16385.0, pow(2.0, 13.0)));
    distanceAdditionLengthReference->numberValue = 13.0;
  }

  reversed = bitReverseLookupTable[(int)(code->numberValue)];
  code->numberValue = ShiftRight4Byte(reversed, 32.0 - 5.0);
}


ByteArray *DeflateDataStaticHuffman(ByteArray *data, double level){
  ByteArray *bytes;
  NumberReference *currentBit;
  double i;
  ByteArray *copy;
  NumberReference *code, *length, *compressedCode, *lengthAdditionLength, *distanceCode;
  NumberReference *distanceReference, *lengthReference, *lengthAddition;
  NumberReference *distanceAdditionReference, *distanceAdditionLengthReference;
  double *bitReverseLookupTable;
  size_t bitReverseLookupTableLength;
  BooleanReference *match;

  code = CreateNumberReference(0.0);
  length = CreateNumberReference(0.0);
  compressedCode = CreateNumberReference(0.0);
  lengthAdditionLength = CreateNumberReference(0.0);
  distanceCode = CreateNumberReference(0.0);
  distanceReference = CreateNumberReference(0.0);
  lengthReference = CreateNumberReference(0.0);
  lengthAddition = CreateNumberReference(0.0);
  distanceAdditionReference = CreateNumberReference(0.0);
  distanceAdditionLengthReference = CreateNumberReference(0.0);
  match = (BooleanReference *)calloc(sizeof(BooleanReference), 1);

  bytes = CreateAndFillByteArray(fmax(ByteArrayLength(data)*1.0, 100.0), 0.0);
  currentBit = CreateNumberReference(0.0);

  bitReverseLookupTable = GenerateBitReverseLookupTable(&bitReverseLookupTableLength, 9.0);

  /* Final block */
  AppendBitsToBytesRight(bytes, currentBit, 1.0, 1.0);
  /* Fixed code */
  AppendBitsToBytesRight(bytes, currentBit, 1.0, 2.0);

  for(i = 0.0; i < ByteArrayLength(data); ){
    FindMatch(data, i, distanceReference, lengthReference, match, level);

    if(match->booleanValue){
      GetDeflateLengthCode(lengthReference->numberValue, compressedCode, lengthAddition, lengthAdditionLength);
      GetDeflateDistanceCode(distanceReference->numberValue, distanceCode, distanceAdditionReference, distanceAdditionLengthReference, bitReverseLookupTable);
    }

    if( !match->booleanValue ){
      GetDeflateStaticHuffmanCode(GetByte(data, i), code, length, bitReverseLookupTable);
      AppendBitsToBytesRight(bytes, currentBit, code->numberValue, length->numberValue);
      i = i + 1.0;
    }else{
      GetDeflateStaticHuffmanCode(compressedCode->numberValue, code, length, bitReverseLookupTable);
      AppendBitsToBytesRight(bytes, currentBit, code->numberValue, length->numberValue);
      AppendBitsToBytesRight(bytes, currentBit, lengthAddition->numberValue, lengthAdditionLength->numberValue);
      AppendBitsToBytesRight(bytes, currentBit, distanceCode->numberValue, 5.0);
      AppendBitsToBytesRight(bytes, currentBit, distanceAdditionReference->numberValue, distanceAdditionLengthReference->numberValue);
      i = i + lengthReference->numberValue;
    }
  }

  /* Stop symbol */
  GetDeflateStaticHuffmanCode(256.0, code, length, bitReverseLookupTable);
  AppendBitsToBytesRight(bytes, currentBit, code->numberValue, length->numberValue);

  copy = (ByteArray *)calloc(sizeof(ByteArray), 1);
  CopyByteArrayRange(bytes, 0.0, ceil(currentBit->numberValue/8.0), copy);
  free(bytes);

  return copy;
}


void FindMatch(ByteArray *data, double pos, NumberReference *distanceReference, NumberReference *lengthReference, BooleanReference *match, double level){
  double i, j;
  double deflateMinMength, deflateMaxLength, deflateMaxDistance;
  double longest, maxLength, distanceForMax;
  double startDistance, matchLength;
  bool done;

  deflateMinMength = 3.0;
  deflateMaxLength = 258.0;

  longest = fmin(pos - 1.0, deflateMaxLength);
  longest = fmin(ByteArrayLength(data) - pos, longest);

  deflateMaxDistance = floor(32768.0/10.0*level);

  startDistance = fmin(pos, deflateMaxDistance);

  if(longest >= deflateMinMength){
    maxLength = 0.0;
    distanceForMax = 0.0;

    for(i = pos - 1.0; i >= pos - startDistance && maxLength != longest; i = i - 1.0){
      matchLength = 0.0;
      done = false;
      for(j = 0.0; j < longest &&  !done ; j = j + 1.0){
        if(GetByte(data, i + j) == GetByte(data, pos + j)){
          matchLength = matchLength + 1.0;
        }else{
          done = true;
        }
      }

      if(matchLength >= deflateMinMength && matchLength > maxLength){
        maxLength = matchLength;
        distanceForMax = pos - i;
      }
    }

    if(maxLength >= deflateMinMength){
      match->booleanValue = true;
      lengthReference->numberValue = maxLength;
      distanceReference->numberValue = distanceForMax;
    }else{
      match->booleanValue = false;
    }
  }else{
    match->booleanValue = false;
  }
}


double *GenerateBitReverseLookupTable(size_t *returnArrayLength, double bits){
  double *table;
  size_t tableLength;
  double i;

  table = (double*)calloc(sizeof(double) * (pow(2.0, bits)), 1);
  tableLength = pow(2.0, bits);

  for(i = 0.0; i < (double)tableLength; i = i + 1.0){
    table[(int)(i)] = ReverseBits(i, 32.0);
  }

  *returnArrayLength = tableLength;
  return table;
}


double ReverseBits(double x, double bits){
  double b, bit, i;

  b = 0.0;

  for(i = 0.0; i < bits; i = i + 1.0){
    b = ShiftLeft4Byte(b, 1.0);
    bit = And4Byte(x, 1.0);
    b = Or4Byte(b, bit);
    x = ShiftRight4Byte(x, 1.0);
  }

  return b;
}


ByteArray *DeflateDataNoCompression(ByteArray *data){
  ByteArray *deflated;
  NumberReference *position;
  double block, i, blocks, blocklength, maxblocksize;

  maxblocksize = pow(2.0, 16.0) - 1.0;
  blocks = ceil(ByteArrayLength(data)/maxblocksize);

  position = CreateNumberReference(0.0);

  deflated = CreateByteArray((1.0 + 4.0)*blocks + ByteArrayLength(data));

  for(block = 0.0; block < blocks; block = block + 1.0){
    if(block + 1.0 == blocks){
      WriteByte(deflated, 1.0, position);
    }else{
      WriteByte(deflated, 0.0, position);
    }
    blocklength = fmin(ByteArrayLength(data) - block*maxblocksize, maxblocksize);
    Write2BytesLE(deflated, blocklength, position);
    Write2BytesLE(deflated, Not2Byte(blocklength), position);

    for(i = 0.0; i < blocklength; i = i + 1.0){
      WriteByte(deflated, GetByte(data, block*maxblocksize + i), position);
    }
  }

  return deflated;
}


void GetDeflateLengthCode(double length, NumberReference *code, NumberReference *lengthAddition, NumberReference *lengthAdditionLength){
  if(length >= 3.0 && length <= 10.0){
    code->numberValue = 257.0 + length - 3.0;
    lengthAdditionLength->numberValue = 0.0;
  }else if(length >= 11.0 && length <= 18.0){
    code->numberValue = 265.0 + floor((length - 11.0)/2.0);
    lengthAddition->numberValue = floor(fmod(length - 11.0, 2.0));
    lengthAdditionLength->numberValue = 1.0;
  }else if(length >= 19.0 && length <= 34.0){
    code->numberValue = 269.0 + floor((length - 19.0)/4.0);
    lengthAddition->numberValue = floor(fmod(length - 19.0, 4.0));
    lengthAdditionLength->numberValue = 2.0;
  }else if(length >= 35.0 && length <= 66.0){
    code->numberValue = 273.0 + floor((length - 35.0)/8.0);
    lengthAddition->numberValue = floor(fmod(length - 35.0, 8.0));
    lengthAdditionLength->numberValue = 3.0;
  }else if(length >= 67.0 && length <= 130.0){
    code->numberValue = 277.0 + floor((length - 67.0)/16.0);
    lengthAddition->numberValue = floor(fmod(length - 67.0, 16.0));
    lengthAdditionLength->numberValue = 4.0;
  }else if(length >= 131.0 && length <= 257.0){
    code->numberValue = 281.0 + floor((length - 131.0)/32.0);
    lengthAddition->numberValue = floor(fmod(length - 131.0, 32.0));
    lengthAdditionLength->numberValue = 5.0;
  }else if(length == 258.0){
    code->numberValue = 285.0;
    lengthAdditionLength->numberValue = 0.0;
  }
}


void AppendBitsToBytesLeft(ByteArray *bytes, NumberReference *nextbit, double data, double length){
  double bytePos, bitPos, segment, part, remove;

  for(; length > 0.0; ){
    bytePos = Truncate(nextbit->numberValue/8.0);
    bitPos = fmod(nextbit->numberValue, 8.0);

    if(length < 8.0 - bitPos){
      part = ShiftLeft4Byte(data, 8.0 - bitPos - length);

      SetByte(bytes, bytePos, Or4Byte(GetByte(bytes, bytePos), part));

      nextbit->numberValue = nextbit->numberValue + length;

      length = 0.0;
    }else{
      segment = 8.0 - bitPos;

      part = ShiftRight4Byte(data, length - segment);
      SetByte(bytes, bytePos, Or4Byte(GetByte(bytes, bytePos), part));
      nextbit->numberValue = nextbit->numberValue + segment;

      remove = ShiftLeft4Byte(part, length - segment);
      data = Xor4Byte(data, remove);

      length = length - segment;
    }
  }
}


void AppendBitsToBytesRight(ByteArray *bytes, NumberReference *nextbit, double data, double length){
    double bytePos, bitPos, segment, part;
    double mask;

    for(; length > 0.0; ) {
        bytePos = Truncate(nextbit->numberValue/8.0);
        bitPos = fmod(nextbit->numberValue, 8.0);

        if(length < 8.0 - bitPos){
            part = ShiftLeft4Byte(data, bitPos);

            SetByte(bytes, bytePos, Or4Byte(GetByte(bytes, bytePos), part));

            nextbit->numberValue = nextbit->numberValue + length;
            length = 0.0;
        } else {
            segment = 8.0 - bitPos;

            mask = 1.0;
            mask = ShiftLeft4Byte(mask, segment);
            mask = mask - 1.0;

            part = And4Byte(mask, data);
            part = ShiftLeft4Byte(part, bitPos);
            SetByte(bytes, bytePos, Or4Byte(GetByte(bytes, bytePos), part));
            nextbit->numberValue = nextbit->numberValue + segment;

            data = ShiftRight4Byte(data, segment);
            length = length - segment;
        }
    }
}



