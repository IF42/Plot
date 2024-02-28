#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define strparam(str) (str), strlen(str)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct RGBABitmapImage
{
  uint32_t *pixels;
  size_t yLength;
  size_t xLength;
}RGBABitmapImage;


typedef struct RGBABitmapImageReference
{
  RGBABitmapImage *image;
}RGBABitmapImageReference;

typedef struct Rectangle
{
  double x1;
  double x2;
  double y1;
  double y2;
}Rectangle;


typedef struct RGBA
{
  double r;
  double g;
  double b;
  double a;
}RGBA;


typedef enum {
    Plot_LineType_Solid
    , Plot_LineType_Dashed
    , Plot_LineType_Dotted
    , Plot_LineType_Dotdash
    , Plot_LineType_Longdash
    , Plot_LineType_Twodash
}Plot_LineType;


typedef enum {
    Plot_PointType_Crosses
    , Plot_PointType_Circles
    , Plot_PointType_Dots
    , Plot_PointType_Triangles
    , Plot_PointType_FilledTriangles
    , Plot_PointType_Pixels
    , Plot_PointType_DotlineToXAxis
}Plot_PointType;


typedef struct ScatterPlotSeries
{
  bool linearInterpolation;
  Plot_PointType pointType;
  Plot_LineType lineType; 
  double lineThickness;
  double *xs;
  size_t xsLength;
  double *ys;
  size_t ysLength;
  RGBA *color;
}ScatterPlotSeries;


typedef struct ScatterPlotSettings
{
  ScatterPlotSeries *scatterPlotSeries;
  size_t scatterPlotSeriesLength;
  bool autoBoundaries;
  double xMax;
  double xMin;
  double yMax;
  double yMin;
  bool autoPadding;
  double xPadding;
  double yPadding;
  char *xLabel;
  char *yLabel;
  char *title;
  bool showGrid;
  RGBA *gridColor;
  bool xAxisAuto;
  bool xAxisTop;
  bool xAxisBottom;
  bool yAxisAuto;
  bool yAxisLeft;
  bool yAxisRight;
  double width;
  double height;
}ScatterPlotSettings;


typedef struct BarPlotSeries
{
  double *ys;
  size_t ysLength;
  RGBA *color;
}BarPlotSeries;


typedef struct StringReference
{
  char *string;
  size_t stringLength;
}StringReference;


typedef struct BarPlotSettings
{
  double width;
  double height;
  bool autoBoundaries;
  double yMax;
  double yMin;
  bool autoPadding;
  double xPadding;
  double yPadding;
  char *title;
  bool showGrid;
  RGBA *gridColor;
  BarPlotSeries *barPlotSeries;
  size_t barPlotSeriesLength;
  char *yLabel;
  bool autoColor;
  bool grayscaleAutoColor;
  bool autoSpacing;
  double groupSeparation;
  double barSeparation;
  bool autoLabels;
  char ** xLabels;
  size_t xLabelsLength;
  bool barBorder;
}BarPlotSettings;


typedef struct BooleanArrayReference
{
  bool *booleanArray;
  size_t booleanArrayLength;
}BooleanArrayReference;


typedef struct BooleanReference
{
  bool booleanValue;
}BooleanReference;


typedef struct CharacterReference
{
  char characterValue;
}CharacterReference;


typedef struct NumberArrayReference
{
  double *numberArray;
  size_t numberArrayLength;
}NumberArrayReference;


typedef struct NumberReference
{
  double numberValue;
}NumberReference;


typedef struct StringArrayReference
{
  StringReference **stringArray;
  size_t stringArrayLength;
}StringArrayReference;


typedef struct ByteArray
{
  uint8_t *bytes;
  size_t bytesLength;
}ByteArray;


void WriteToFile(ByteArray *data, char *filename);


typedef struct Chunk
{
  double length;
  char *type;
  size_t typeLength;
  ByteArray *data;
  double crc;
}Chunk;


typedef struct IHDR
{
  double Width;
  double Height;
  double BitDepth;
  double ColourType;
  double CompressionMethod;
  double FilterMethod;
  double InterlaceMethod;
}IHDR;


typedef struct PHYS
{
  double pixelsPerMeter;
}PHYS;


typedef struct ZLIBStruct
{
  double CMF;
  double CM;
  double CINFO;
  double FLG;
  double FCHECK;
  double FDICT;
  double FLEVEL;
  ByteArray *CompressedDataBlocks;
  double Adler32CheckValue;
}ZLIBStruct;


typedef struct PNGImage
{
  double *signature;
  size_t signatureLength;
  IHDR *ihdr;
  ZLIBStruct *zlibStruct;
  bool physPresent;
  PHYS *phys;
}PNGImage;


typedef struct DynamicArrayCharacters
{
  char *array;
  size_t arrayLength;
  double length;
}DynamicArrayCharacters;


typedef struct LinkedListNodeStrings
{
  bool end;
  char *value;
  size_t valueLength;
  struct LinkedListNodeStrings *next;
}LinkedListNodeStrings;


typedef struct LinkedListStrings
{
  struct LinkedListNodeStrings *first;
  struct LinkedListNodeStrings *last;
}LinkedListStrings;


typedef struct LinkedListNodeNumbers
{
  struct LinkedListNodeNumbers *next;
  bool end;
  double value;
}LinkedListNodeNumbers;


typedef struct LinkedListNumbers
{
  struct LinkedListNodeNumbers *first;
  struct LinkedListNodeNumbers *last;
}LinkedListNumbers;


typedef struct LinkedListCharacters
{
  struct LinkedListNodeCharacters *first;
  struct LinkedListNodeCharacters *last;
}LinkedListCharacters;


typedef struct LinkedListNodeCharacters
{
  bool end;
  char value;
  struct LinkedListNodeCharacters *next;
}LinkedListNodeCharacters;


typedef struct DynamicArrayNumbers
{
  double *array;
  size_t arrayLength;
  double length;
}DynamicArrayNumbers;


typedef struct ByteArrayReference
{
  ByteArray *bytes;
}ByteArrayReference;


bool Loess(double *xs, size_t xsLength, double *ys, size_t ysLength, double bandwidth, double robustnessIters, double accuracy, NumberArrayReference *resultXs);
bool Lowess(double *xs, size_t xsLength, double *ys, size_t ysLength, double *weights, size_t weightsLength, double bandwidth, double robustnessIters, double accuracy, NumberArrayReference *resultXs);
void AssignNumberArray(double *as, size_t asLength, double *bs, size_t bsLength);
double FindNextNonZeroElement(double *array, size_t arrayLength, double offset);
double Tricube(double x);

bool CropLineWithinBoundary(NumberReference *x1Ref, NumberReference *y1Ref, NumberReference *x2Ref, NumberReference *y2Ref, double xMin, double xMax, double yMin, double yMax);

RGBA **Get8HighContrastColors(size_t *returnArrayLength);

void DrawFilledRectangleWithBorder(RGBABitmapImage *image, double x, double y, double w, double h, RGBA *borderColor, RGBA *fillColor);
RGBABitmapImageReference *CreateRGBABitmapImageReference();

bool RectanglesOverlap(Rectangle *r1, Rectangle *r2);
Rectangle *CreateRectangle(double x1, double y1, double x2, double y2);
void CopyRectangleValues(Rectangle *rd, Rectangle *rs);

void DrawXLabelsForPriority(double p, double xMin, double oy, double xMax, double xPixelMin, double xPixelMax, NumberReference *nextRectangle, RGBA *gridLabelColor, RGBABitmapImage *canvas, double *xGridPositions, size_t xGridPositionsLength, StringArrayReference *xLabels, NumberArrayReference *xLabelPriorities, Rectangle **occupied, bool textOnBottom);
void DrawYLabelsForPriority(double p, double yMin, double ox, double yMax, double yPixelMin, double yPixelMax, NumberReference *nextRectangle, RGBA *gridLabelColor, RGBABitmapImage *canvas, double *yGridPositions, size_t yGridPositionsLength, StringArrayReference *yLabels, NumberArrayReference *yLabelPriorities, Rectangle **occupied, bool textOnLeft);
double *ComputeGridLinePositions(size_t *returnArrayLength, double cMin, double cMax, StringArrayReference *labels, NumberArrayReference *priorities);
double MapYCoordinate(double y, double yMin, double yMax, double yPixelMin, double yPixelMax);
double MapXCoordinate(double x, double xMin, double xMax, double xPixelMin, double xPixelMax);
double MapXCoordinateAutoSettings(double x, RGBABitmapImage *image, double *xs, size_t xsLength);
double MapYCoordinateAutoSettings(double y, RGBABitmapImage *image, double *ys, size_t ysLength);
double MapXCoordinateBasedOnSettings(double x, ScatterPlotSettings *settings);
double MapYCoordinateBasedOnSettings(double y, ScatterPlotSettings *settings);
double GetDefaultPaddingPercentage();

void DrawText(RGBABitmapImage *canvas, double x, double y, char *text, RGBA *color);
void DrawTextUpwards(RGBABitmapImage *canvas, double x, double y, char *text, size_t textLength, RGBA *color);


#define GetDefaultScatterPlotSettings()              \
    (ScatterPlotSettings) {                          \
          .autoBoundaries = true                     \
          , .xMax = 0.0                              \
          , .xMin = 0.0                              \
          , .yMax = 0.0                              \
          , .yMin = 0.0                              \
          , .autoPadding = true                      \
          , .xPadding = 0.0                          \
          , .yPadding = 0.0                          \
          , .title = ""                              \
          , .xLabel = ""                             \
          , .yLabel = ""                             \
          , .scatterPlotSeries = NULL                \
          , .scatterPlotSeriesLength = 0.0           \
          , .showGrid = true                         \
          , .gridColor = GetGray(0.1)                \
          , .xAxisAuto = true                        \
          , .xAxisTop = false                        \
          , .xAxisBottom = false                     \
          , .yAxisAuto = true                        \
          , .yAxisLeft = false                       \
          , .yAxisRight = false}


#define GetDefaultScatterPlotSeriesSettings()       \
    (ScatterPlotSeries) {                           \
          .linearInterpolation = true               \
          , .pointType = Plot_PointType_Pixels      \
          , .lineType = Plot_LineType_Solid         \
          , .lineThickness = 1.0                    \
          , .xs = NULL                              \
          , .xsLength = 0.0                         \
          , .ys = NULL                              \
          , .ysLength = 0.0                         \
          , .color = GetBlack()}


bool DrawScatterPlot(RGBABitmapImageReference *canvasReference, double width, double height, double *xs, size_t xsLength, double *ys, size_t ysLength);
bool DrawScatterPlotFromSettings(RGBABitmapImageReference *canvasReference, ScatterPlotSettings *settings);
void ComputeBoundariesBasedOnSettings(ScatterPlotSettings *settings, Rectangle *boundaries);
bool ScatterPlotFromSettingsValid(ScatterPlotSettings *settings);


#define GetDefaultBarPlotSettings()             \
    (BarPlotSettings) {                         \
      .width = 800.0                            \
      , .height = 600.0                         \
      , .autoBoundaries = true                  \
      , .yMax = 0.0                             \
      , .yMin = 0.0                             \
      , .autoPadding = true                     \
      , .xPadding = 0.0                         \
      , .yPadding = 0.0                         \
      , .title = ""                             \
      , .yLabel = ""                            \
      , .barPlotSeries = NULL                   \
      , .barPlotSeriesLength = 0.0              \
      , .showGrid = true                        \
      , .gridColor = GetGray(0.1)               \
      , .autoColor = true                       \
      , .grayscaleAutoColor = false             \
      , .autoSpacing = true                     \
      , .groupSeparation = 0.0                  \
      , .barSeparation = 0.0                    \
      , .autoLabels = true                      \
      , .xLabels = NULL                         \
      , .xLabelsLength = 0.0                    \
      , .barBorder = false}


#define GetDefaultBarPlotSeriesSettings()       \
    (BarPlotSeries) {                           \
        .ys  = NULL                             \
        , .ysLength = 0.0                       \
        , .color = GetBlack()}


RGBABitmapImage *DrawBarPlotNoErrorCheck(double width, double height, double *ys, size_t ysLength);
bool DrawBarPlot(RGBABitmapImageReference *canvasReference, double width, double height, double *ys, size_t ysLength);
bool DrawBarPlotFromSettings(RGBABitmapImageReference *canvasReference, BarPlotSettings *settings);
bool BarPlotSettingsIsValid(BarPlotSettings *settings);

double GetMinimum(double *data, size_t dataLength);
double GetMaximum(double *data, size_t dataLength);

double RoundToDigits(double element, double digitsAfterPoint);

RGBA *GetBlack();
RGBA *GetWhite();
RGBA *GetTransparent();
RGBA *GetGray(double percentage);
RGBA *CreateRGBColor(double r, double g, double b);
RGBA *CreateRGBAColor(double r, double g, double b, double a);

RGBABitmapImage *CreateImage(double w, double h, RGBA *color);
void DeleteImage(RGBABitmapImage *image);
double ImageWidth(RGBABitmapImage *image);
double ImageHeight(RGBABitmapImage *image);
void SetPixel(RGBABitmapImage *image, double x, double y, RGBA *color);
void DrawPixel(RGBABitmapImage *image, double x, double y, RGBA *color);
double CombineAlpha(double as, double ad);
double AlphaBlend(double cs, double as, double cd, double ad, double ao);
void DrawHorizontalLine1px(RGBABitmapImage *image, double x, double y, double length, RGBA *color);
void DrawVerticalLine1px(RGBABitmapImage *image, double x, double y, double height, RGBA *color);
void DrawRectangle1px(RGBABitmapImage *image, double x, double y, double width, double height, RGBA *color);
void DrawImageOnImage(RGBABitmapImage *dst, RGBABitmapImage *src, double topx, double topy);
void DrawLine1px(RGBABitmapImage *image, double x0, double y0, double x1, double y1, RGBA *color);
void XiaolinWusLineAlgorithm(RGBABitmapImage *image, double x0, double y0, double x1, double y1, RGBA *color);
double OneMinusFractionalPart(double x);
double FractionalPart(double x);
RGBA *SetBrightness(RGBA *color, double newBrightness);
void DrawQuadraticBezierCurve(RGBABitmapImage *image, double x0, double y0, double cx, double cy, double x1, double y1, RGBA *color);
void QuadraticBezierPoint(double x0, double y0, double cx, double cy, double x1, double y1, double t, NumberReference *x, NumberReference *y);
void DrawCubicBezierCurve(RGBABitmapImage *image, double x0, double y0, double c0x, double c0y, double c1x, double c1y, double x1, double y1, RGBA *color);
void CubicBezierPoint(double x0, double y0, double c0x, double c0y, double c1x, double c1y, double x1, double y1, double t, NumberReference *x, NumberReference *y);
RGBABitmapImage *CopyImage(RGBABitmapImage *image);
RGBA *GetImagePixel(RGBABitmapImage *image, double x, double y);
RGBA GetImagePixelStruct(RGBABitmapImage *image, double x, double y);
void HorizontalFlip(RGBABitmapImage *img);
void DrawFilledRectangle(RGBABitmapImage *image, double x, double y, double w, double h, RGBA *color);
RGBABitmapImage *RotateAntiClockwise90Degrees(RGBABitmapImage *image);
void DrawCircle(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color);
void BresenhamsCircleDrawingAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color);
void DrawCircleMidpointAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color);
void DrawCircleBasicAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color);
void DrawFilledCircle(RGBABitmapImage *canvas, double x, double y, double r, RGBA *color);
void DrawFilledCircleMidpointAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color);
void DrawFilledCircleBasicAlgorithm(RGBABitmapImage *canvas, double xCenter, double yCenter, double radius, RGBA *color);
void DrawTriangle(RGBABitmapImage *canvas, double xCenter, double yCenter, double height, RGBA *color);
void DrawFilledTriangle(RGBABitmapImage *canvas, double xCenter, double yCenter, double height, RGBA *color);
void DrawLine(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, double thickness, RGBA *color);
void DrawLineBresenhamsAlgorithmThick(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, double thickness, RGBA *color);
void DrawLineBresenhamsAlgorithm(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, RGBA *color);
void DrawLineBresenhamsAlgorithmThickPatterned(RGBABitmapImage *canvas, double x1, double y1, double x2, double y2, double thickness, bool *pattern, size_t patternLength, NumberReference *offset, RGBA *color);

bool *GetLinePattern5(size_t *returnArrayLength);
bool *GetLinePattern4(size_t *returnArrayLength);
bool *GetLinePattern3(size_t *returnArrayLength);
bool *GetLinePattern2(size_t *returnArrayLength);
bool *GetLinePattern1(size_t *returnArrayLength);

RGBABitmapImage *Blur(RGBABitmapImage *src, double pixels);
RGBA *CreateBlurForPoint(RGBABitmapImage *src, double x, double y, double pixels);

char *CreateStringDecimalFromNumber(size_t *returnArrayLength, double decimal);
char *GetDigitCharacterTable(size_t *returnArrayLength);

bool CreateNumberFromDecimalStringWithCheck(char *string, size_t stringLength, NumberReference *decimalReference);
double CreateNumberFromDecimalString(char *string, size_t stringLength);
bool CreateNumberFromStringWithCheck(char *string, size_t stringLength, double base, NumberReference *numberReference);
double CreateNumberFromParts(double base, bool numberIsPositive, double *beforePoint, size_t beforePointLength, double *afterPoint, size_t afterPointLength, bool exponentIsPositive, double *exponent, size_t exponentLength);
bool ExtractPartsFromNumberString(char *n, size_t nLength, double base, BooleanReference *numberIsPositive, NumberArrayReference *beforePoint, NumberArrayReference *afterPoint, BooleanReference *exponentIsPositive, NumberArrayReference *exponent);
double GetNumberFromNumberCharacterForBase(char c, double base);
bool CharacterIsNumberCharacterInBase(char c, double base);
double *StringToNumberArray(size_t *returnArrayLength, char *str, size_t strLength);
bool StringToNumberArrayWithCheck(char *str, size_t strLength, NumberArrayReference *numberArrayReference);

void QuickSortNumbers(double *list, size_t listLength);
void QuickSortNumbersBounds(double *A, size_t ALength, double lo, double hi);
double *QuickSortNumbersWithIndexes(size_t *returnArrayLength, double *A, size_t ALength);
void QuickSortNumbersBoundsWithIndexes(double *A, size_t ALength, double *indexes, size_t indexesLength, double lo, double hi);

void aFillNumberArray(double *a, size_t aLength, double value);
void aFillString(char *a, size_t aLength, char value);
void aFillBooleanArray(bool *a, size_t aLength, bool value);
bool aFillNumberArrayRange(double *a, size_t aLength, double value, double from, double to);
bool aFillBooleanArrayRange(bool *a, size_t aLength, bool value, double from, double to);
bool aFillStringRange(char *a, size_t aLength, char value, double from, double to);
double *aCopyNumberArray(size_t *returnArrayLength, double *a, size_t aLength);
bool *aCopyBooleanArray(size_t *returnArrayLength, bool *a, size_t aLength);
char *aCopyString(size_t *returnArrayLength, char *a, size_t aLength);
bool aCopyNumberArrayRange(double *a, size_t aLength, double from, double to, NumberArrayReference *copyReference);
bool aCopyBooleanArrayRange(bool *a, size_t aLength, double from, double to, BooleanArrayReference *copyReference);
bool aCopyStringRange(char *a, size_t aLength, double from, double to, StringReference *copyReference);
void aSwapElementsOfStringArray(StringArrayReference *A, double ai, double bi);


BooleanReference *CreateBooleanReference(bool value);
BooleanArrayReference *CreateBooleanArrayReferenceLengthValue(double length, bool value);
void FreeBooleanArrayReference(BooleanArrayReference *booleanArrayReference);
CharacterReference *CreateCharacterReference(char value);
NumberReference *CreateNumberReference(double value);
NumberArrayReference *CreateNumberArrayReference(double *value, size_t valueLength);
NumberArrayReference *CreateNumberArrayReferenceLengthValue(double length, double value);
void FreeNumberArrayReference(NumberArrayReference *numberArrayReference);
StringReference *CreateStringReference(char *value, size_t valueLength);
StringReference *CreateStringReferenceLengthValue(double length, char value);
void FreeStringReference(StringReference *stringReference);
StringArrayReference *CreateStringArrayReference(StringReference **strings, size_t stringsLength);
StringArrayReference *CreateStringArrayReferenceLengthValue(double length, char *value, size_t valueLength);
void FreeStringArrayReference(StringArrayReference *stringArrayReference);

char *DigitDataBase16(size_t *returnArrayLength);
void DrawDigitCharacter(RGBABitmapImage *image, double topx, double topy, double digit);

char *GetPixelFontData(size_t *returnArrayLength);
void DrawAsciiCharacter(RGBABitmapImage *image, double topx, double topy, char a, RGBA *color);

double GetTextWidth(size_t textLength);
#define GetTextHeight 13.0

ByteArray *ConvertToPNG(RGBABitmapImage *image);
ByteArray *ConvertToPNGGrayscale(RGBABitmapImage *image);
PHYS *PysicsHeader(double pixelsPerMeter);
ByteArray *ConvertToPNGWithOptions(RGBABitmapImage *image, double colorType, bool setPhys, double pixelsPerMeter, double compressionLevel);
ByteArray *PNGSerializeChunks(PNGImage *png);
double PNGIDATLength(PNGImage *png);
double PNGHeaderLength();
ByteArray *GetPNGColorData(RGBABitmapImage *image);
ByteArray *GetPNGColorDataGreyscale(RGBABitmapImage *image);
IHDR *PNGHeader(RGBABitmapImage *image, double colortype);
double *PNGSignature(size_t *returnArrayLength);
double *PNGReadDataChunks(size_t *returnArrayLength, Chunk **cs, size_t csLength);
bool PNGReadHeader(RGBABitmapImage *image, Chunk **cs, size_t csLength, StringReference *errorMessages);
Chunk **PNGReadChunks(size_t *returnArrayLength, ByteArray *data, NumberReference *position);
Chunk *PNGReadChunk(ByteArray *data, NumberReference *position);

void WriteStringToStingStream(char *stream, NumberReference *index, char *src, size_t srcLength);
void WriteCharacterToStingStream(char *stream, NumberReference *index, char src);
void WriteBooleanToStingStream(char *stream, NumberReference *index, bool src);

bool SubstringWithCheck(char *string, size_t stringLength, double from, double to, StringReference *stringReference);
char *AppendString(size_t *returnArrayLength, char *s1, size_t s1Length, char *s2, size_t s2Length);
char *ConcatenateString(size_t *returnArrayLength, char *s1, size_t s1Length, char *s2, size_t s2Length);
char *AppendCharacter(size_t *returnArrayLength, char *string, size_t stringLength, char c);
char *ConcatenateCharacter(size_t *returnArrayLength, char *string, size_t stringLength, char c);
bool SubstringEqualsWithCheck(char *string, size_t stringLength, double from, char *substring, size_t substringLength, BooleanReference *equalsReference);
bool SubstringEquals(char *string, size_t stringLength, double from, char *substring, size_t substringLength);
void ToUpperCase(char *string, size_t stringLength);
void ToLowerCase(char *string, size_t stringLength);
char *Trim(size_t *returnArrayLength, char *string, size_t stringLength);
bool EndsWith(char *string, size_t stringLength, char *end, size_t endLength);
StringReference **SplitByString(size_t *returnArrayLength, char *toSplit, size_t toSplitLength, char *splitBy, size_t splitByLength);
bool StringIsBefore(char *a, size_t aLength, char *b, size_t bLength);


StringReference **AddString(size_t *returnArrayLength, StringReference **list, size_t listLength, StringReference *a);
void AddStringRef(StringArrayReference *list, StringReference *i);
StringReference **RemoveString(size_t *returnArrayLength, StringReference **list, size_t listLength, double n);
StringReference *GetStringRef(StringArrayReference *list, double i);
void RemoveStringRef(StringArrayReference *list, double i);


DynamicArrayCharacters *CreateDynamicArrayCharacters();
DynamicArrayCharacters *CreateDynamicArrayCharactersWithInitialCapacity(double capacity);
void DynamicArrayAddCharacter(DynamicArrayCharacters *da, char value);
void DynamicArrayCharactersIncreaseSize(DynamicArrayCharacters *da);
bool DynamicArrayCharactersDecreaseSizeNecessary(DynamicArrayCharacters *da);
void DynamicArrayCharactersDecreaseSize(DynamicArrayCharacters *da);
double DynamicArrayCharactersIndex(DynamicArrayCharacters *da, double index);
double DynamicArrayCharactersLength(DynamicArrayCharacters *da);
void DynamicArrayInsertCharacter(DynamicArrayCharacters *da, double index, char value);
bool DynamicArrayCharacterSet(DynamicArrayCharacters *da, double index, char value);
void DynamicArrayRemoveCharacter(DynamicArrayCharacters *da, double index);
void FreeDynamicArrayCharacters(DynamicArrayCharacters *da);
char *DynamicArrayCharactersToArray(size_t *returnArrayLength, DynamicArrayCharacters *da);
DynamicArrayCharacters *ArrayToDynamicArrayCharactersWithOptimalSize(char *array, size_t arrayLength);
DynamicArrayCharacters *ArrayToDynamicArrayCharacters(char *array, size_t arrayLength);
bool DynamicArrayCharactersEqual(DynamicArrayCharacters *a, DynamicArrayCharacters *b);
LinkedListCharacters *DynamicArrayCharactersToLinkedList(DynamicArrayCharacters *da);
DynamicArrayCharacters *LinkedListToDynamicArrayCharacters(LinkedListCharacters *ll);

bool *AddBoolean(size_t *returnArrayLength, bool *list, size_t listLength, bool a);
void AddBooleanRef(BooleanArrayReference *list, bool i);
bool *RemoveBoolean(size_t *returnArrayLength, bool *list, size_t listLength, double n);
bool GetBooleanRef(BooleanArrayReference *list, double i);
void RemoveDecimalRef(BooleanArrayReference *list, double i);


LinkedListStrings *CreateLinkedListString();
void LinkedListAddString(LinkedListStrings *ll, char *value, size_t valueLength);
StringReference **LinkedListStringsToArray(size_t *returnArrayLength, LinkedListStrings *ll);
double LinkedListStringsLength(LinkedListStrings *ll);
void FreeLinkedListString(LinkedListStrings *ll);


LinkedListNumbers *CreateLinkedListNumbers();
LinkedListNumbers **CreateLinkedListNumbersArray(size_t *returnArrayLength, double length);
void LinkedListAddNumber(LinkedListNumbers *ll, double value);
double LinkedListNumbersLength(LinkedListNumbers *ll);
double LinkedListNumbersIndex(LinkedListNumbers *ll, double index);
void LinkedListInsertNumber(LinkedListNumbers *ll, double index, double value);
void LinkedListSet(LinkedListNumbers *ll, double index, double value);
void LinkedListRemoveNumber(LinkedListNumbers *ll, double index);
void FreeLinkedListNumbers(LinkedListNumbers *ll);
void FreeLinkedListNumbersArray(LinkedListNumbers **lls, size_t llsLength);
double *LinkedListNumbersToArray(size_t *returnArrayLength, LinkedListNumbers *ll);
LinkedListNumbers *ArrayToLinkedListNumbers(double *array, size_t arrayLength);
bool LinkedListNumbersEqual(LinkedListNumbers *a, LinkedListNumbers *b);

LinkedListCharacters *CreateLinkedListCharacter();
void LinkedListAddCharacter(LinkedListCharacters *ll, char value);
char *LinkedListCharactersToArray(size_t *returnArrayLength, LinkedListCharacters *ll);
double LinkedListCharactersLength(LinkedListCharacters *ll);
void FreeLinkedListCharacter(LinkedListCharacters *ll);
void LinkedListCharactersAddString(LinkedListCharacters *ll, char *str, size_t strLength);

DynamicArrayNumbers *CreateDynamicArrayNumbers();
DynamicArrayNumbers *CreateDynamicArrayNumbersWithInitialCapacity(double capacity);
void DynamicArrayAddNumber(DynamicArrayNumbers *da, double value);
void DynamicArrayNumbersIncreaseSize(DynamicArrayNumbers *da);
bool DynamicArrayNumbersDecreaseSizeNecessary(DynamicArrayNumbers *da);
void DynamicArrayNumbersDecreaseSize(DynamicArrayNumbers *da);
double DynamicArrayNumbersIndex(DynamicArrayNumbers *da, double index);
double DynamicArrayNumbersLength(DynamicArrayNumbers *da);
void DynamicArrayInsertNumber(DynamicArrayNumbers *da, double index, double value);
bool DynamicArrayNumberSet(DynamicArrayNumbers *da, double index, double value);
void DynamicArrayRemoveNumber(DynamicArrayNumbers *da, double index);
void FreeDynamicArrayNumbers(DynamicArrayNumbers *da);
double *DynamicArrayNumbersToArray(size_t *returnArrayLength, DynamicArrayNumbers *da);
DynamicArrayNumbers *ArrayToDynamicArrayNumbersWithOptimalSize(double *array, size_t arrayLength);
DynamicArrayNumbers *ArrayToDynamicArrayNumbers(double *array, size_t arrayLength);
bool DynamicArrayNumbersEqual(DynamicArrayNumbers *a, DynamicArrayNumbers *b);
LinkedListNumbers *DynamicArrayNumbersToLinkedList(DynamicArrayNumbers *da);
DynamicArrayNumbers *LinkedListToDynamicArrayNumbers(LinkedListNumbers *ll);
double DynamicArrayNumbersIndexOf(DynamicArrayNumbers *arr, double n, BooleanReference *foundReference);
bool DynamicArrayNumbersIsInArray(DynamicArrayNumbers *arr, double n);

char *AddCharacter(size_t *returnArrayLength, char *list, size_t listLength, char a);
void AddCharacterRef(StringReference *list, char i);
char *RemoveCharacter(size_t *returnArrayLength, char *list, size_t listLength, double n);
char GetCharacterRef(StringReference *list, double i);
void RemoveCharacterRef(StringReference *list, double i);

ByteArray *ReadXbytes(ByteArray *data, NumberReference *position, double length);
double Read4bytesBE(ByteArray *data, NumberReference *position);
double Read2bytesBE(ByteArray *data, NumberReference *position);
double ReadByte(ByteArray *data, NumberReference *position);
double Read4bytesLE(ByteArray *data, NumberReference *position);
void WriteByte(ByteArray *data, double b, NumberReference *position);
void Write2BytesLE(ByteArray *data, double b, NumberReference *position);
void Write4BytesLE(ByteArray *data, double b, NumberReference *position);
void Write2BytesBE(ByteArray *data, double b, NumberReference *position);
void Write4BytesBE(ByteArray *data, double b, NumberReference *position);
void WriteStringBytes(ByteArray *data, char *cs, size_t csLength, NumberReference *position);
double BytesRound(double x);
double *ByteArrayToNumberArray(size_t *returnArrayLength, ByteArray *src);
ByteArray *NumberArrayToByteArray(double *src, size_t srcLength);
bool ByteArraysEqual(ByteArray *a, ByteArray *b);
ByteArray *CopyByteArray(ByteArray *a);
double ByteArrayLength(ByteArray *response);
ByteArray *CreateAndFillByteArray(double length, double value);
ByteArray *CreateByteArray(double length);
void SetByte(ByteArray *array, double index, double value);
double GetByte(ByteArray *array, double index);
void AssertByteArraysEqual(ByteArray *a, ByteArray *b, NumberReference *failures);
void FreeByteArray(ByteArray *byteArray);
bool CopyByteArrayRange(ByteArray *a, double from, double to, ByteArray *b);

char *BytesToTextBase16(size_t *returnArrayLength, double *bytes, size_t bytesLength);
double *TextToBytesBase16(size_t *returnArrayLength, char *string, size_t stringLength);
StringReference **GenerateBase16ByteCombinations(size_t *returnArrayLength);

double *MakeCRC32Table(size_t *returnArrayLength);
double CalculateCRC32(ByteArray *buf);
double CRC32OfInterval(ByteArray *data, double from, double length);

ZLIBStruct *ZLibCompressNoCompression(ByteArray *data);
ZLIBStruct *ZLibCompressStaticHuffman(ByteArray *data, double level);


double And4Byte(double n1, double n2);
double Or4Byte(double n1, double n2);
double Xor4Byte(double n1, double n2);
double Not2Byte(double b);
double ShiftLeft4Byte(double b, double amount);
double ShiftLeftByte(double b, double amount);
double ShiftRight4Byte(double b, double amount);


ByteArray *Pack(ByteArray *data, double level);
ByteArray *Unpack(ByteArray *data);

ByteArray *DeflateDataStaticHuffman(ByteArray *data, double level);
void FindMatch(ByteArray *data, double pos, NumberReference *distanceReference, NumberReference *lengthReference, BooleanReference *match, double level);
double *GenerateBitReverseLookupTable(size_t *returnArrayLength, double bits);
double ReverseBits(double x, double bits);
ByteArray *DeflateDataNoCompression(ByteArray *data);
void GetDeflateLengthCode(double length, NumberReference *code, NumberReference *lengthAddition, NumberReference *lengthAdditionLength);
void AppendBitsToBytesLeft(ByteArray *bytes, NumberReference *nextbit, double data, double length);
void AppendBitsToBytesRight(ByteArray *bytes, NumberReference *nextbit, double data, double length);





