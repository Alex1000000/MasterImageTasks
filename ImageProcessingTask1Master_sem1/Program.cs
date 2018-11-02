using ImageReadCS;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace ImageProcessingTask1Master_sem1
{

    class Program
    {static bool answ = false;
    static int reqXstat = 0;
    static int reqYstat = 0;
         static void InvertImage(ColorFloatImage image)
        {
            for (int y = 0; y < image.Height; y++)
                for (int x = 0; x < image.Width ; x++)
                {
                    float r = image[x, y].r;
                    float g = image[x, y].g;
                    float b = image[x, y].b;
                    ColorFloatPixel pixel = new ColorFloatPixel();
                    pixel.r = 255 - r;
                    pixel.g = 255 - g;
                    pixel.b = 255 - b;
                    image[x, y] = pixel;
                    //image[image.Width - 1 - x, y] = p;
                }
        }
         static void MirrowYImage(ColorFloatImage image)
         {
             for (int y = 0; y < image.Height; y++)
                 for (int x = 0; x < image.Width/2; x++)
                 {
                     ColorFloatPixel pixel = new ColorFloatPixel();
                     pixel.r = image[x, y].r;
                     pixel.g = image[x, y].g;
                     pixel.b = image[x, y].b;
                     image[x, y] = image[image.Width - 1 - x, y];
                     image[image.Width - 1 - x, y] = pixel;
 
                 }
         }
         static void MirrowXImage(ColorFloatImage image)
         {
             for (int y = 0; y < image.Height/2; y++)
                 for (int x = 0; x < image.Width; x++)
                 {
                     ColorFloatPixel pixel = new ColorFloatPixel();
                     pixel.r = image[x, y].r;
                     pixel.g = image[x, y].g;
                     pixel.b = image[x, y].b;
                     image[x, y] = image[ x, image.Height - 1 - y];
                     image[ x, image.Height - 1 - y] = pixel;

                 }
         }



         static void rotateImage(ColorFloatImage image, bool isCW, float angleDeg, string fileOutName )
         {
             float angle = angleDeg * (float)Math.PI / 180;
             int x = 0;
             int y = 0;
             int iCentreX = image.Width / 2;
             int iCentreY = image.Height / 2;
             //int deltaPlus=500;
             int bilinearRotationWidth = (int)Math.Round( image.Width * Math.Abs(Math.Cos(angle)) + image.Height * Math.Abs(Math.Sin(angle)));//image.Width + deltaPlus;
             int bilinearRotationHeight = (int)Math.Round(image.Width * Math.Abs(Math.Sin(angle)) + image.Height * Math.Abs(Math.Cos(angle)));//image.Height + deltaPlus;
             ColorFloatImage bilinearInterpolatationForRotation = new ColorFloatImage(bilinearRotationWidth, bilinearRotationHeight);
             int radius = 0;//1;
             ColorFloatImage imageCopyExtended = new ColorFloatImage (image.Width + 2 * radius, image.Height + 2 * radius);
             
             if (angleDeg % 90 == 0)
             {
                 
                 if (angleDeg % 360 == 0)
                 {
                     ImageIO.ImageToFile(image, fileOutName);
                     return;
                 }
                 else if (angleDeg % 180 == 0)
                 {
                     for (int i = 0; i < image.Width; i++)
                     {
                         for (int j = 0; j < image.Height; j++)
                         {
                             bilinearInterpolatationForRotation[image.Width - i - 1, image.Height - j - 1] = image[i, j];
                         }
                     }
                     ImageIO.ImageToFile(bilinearInterpolatationForRotation, fileOutName);
                     return;
                 }
                 else if (angleDeg % 270 == 0)
                 {
                     int isPlus270 = (int)angleDeg / 270;
                     if (((isPlus270 > 0) && (isPlus270 % 2 == 0)) || ((isPlus270 < 0) && (isPlus270 % 2 != 0)))//+90
                     {
                         for (int i = 0; i < image.Width; i++)
                         {
                             for (int j = 0; j < image.Height; j++)
                             {
                                 bilinearInterpolatationForRotation[image.Height - j - 1, i] = image[i, j];
                             }
                         }

                     }
                     else
                     {
                         for (int i = 0; i < image.Width; i++)
                         {
                             for (int j = 0; j < image.Height; j++)
                             {
                                 bilinearInterpolatationForRotation[j, image.Width - i - 1] = image[i, j];
                             }
                         }
                     }
                     ImageIO.ImageToFile(bilinearInterpolatationForRotation, fileOutName);
                     return;
                 }
                 
                 //+90
                 int isPlus90 = (int)angleDeg / 90;

                 if (((isPlus90 > 0) && (isPlus90 % 2 != 0)) || ((isPlus90 < 0) && (isPlus90 % 2 == 0)))//+90
                 {
                     for (int i = 0; i < image.Width; i++)
                     {
                         for (int j = 0; j < image.Height; j++)
                         {
                             bilinearInterpolatationForRotation[image.Height - j - 1, i] = image[i, j];
                         }
                     }

                 }
                 else
                 {
                     for (int i = 0; i < image.Width; i++)
                     {
                         for (int j = 0; j < image.Height; j++)
                         {
                             bilinearInterpolatationForRotation[j, image.Width - i - 1] = image[i, j];
                         }
                     }
                 }

                 ImageIO.ImageToFile(bilinearInterpolatationForRotation, fileOutName);
                 return;
             }
             
             //mirrow edges
             //copy center
             //копируем в середину
             for (int k = 0; k < image.Width; k++)
             {
                 for (int l = 0; l < image.Height; l++)
                 {
                     imageCopyExtended[k + radius, l + radius] = image[k, l];
                 }
             }
             //копируем края
             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[k, radius-l-1] = imageCopyExtended[k, l + radius];

                 }
             }

             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[k, image.Height + radius +image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[radius- k-1, l] = imageCopyExtended[k + radius, l];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[image.Width + radius+image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                 }
             }
             /////////////////////////////////////////////////////
             //копируем уголки
             for (int k = 0; k < radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[radius-k-1,radius -l-1] = imageCopyExtended[k + radius, l + radius];

                 }
             }

             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                 }
             }
             ImageIO.ImageToFile(imageCopyExtended, fileOutName);


             if ((isCW)||(!isCW))
             {
                 for (int i = 0; i < bilinearRotationHeight; i++)
                 {
                     for (int j = 0; j < bilinearRotationWidth; j++)
                     {

                         // convert raster to Cartesian
                         x = j - (int)bilinearRotationWidth/2;
                         y = (int)bilinearRotationHeight/2 - i;

                         // convert Cartesian to polar
                         float fDistance = (float)Math.Sqrt(x * x + y * y);
                         float fPolarAngle = 0.0f;
                         if (x == 0)
                         {
                             if (y == 0)
                             {
                                 // centre of image, no rotation needed
                                 //bilinearInterpolatationForRotation[j, i] = image[x + image.Width / 2, y + image.Width / 2];
                                 bilinearInterpolatationForRotation[j, i] =imageCopyExtended[x + image.Width / 2+radius, y + image.Width / 2+radius];
                                 continue;
                             }
                             else if (y < 0)
                             {
                                 fPolarAngle = (float)(1.5f * Math.PI);
                             }
                             else
                             {
                                 fPolarAngle = (float)(0.5 * Math.PI);
                             }
                         }
                         else
                         {
                             fPolarAngle = (float)Math.Atan2((double)y, (double)x);
                         }

                         // the crucial rotation part
                         // "reverse" rotate, so minus instead of plus
                         fPolarAngle += angle;

                         // convert polar to Cartesian
                         float fTrueX = (float)(fDistance * Math.Cos(fPolarAngle));
                         float fTrueY = (float)(fDistance * Math.Sin(fPolarAngle));

                         // convert Cartesian to raster
                         fTrueX = (float)(fTrueX + (double)iCentreX);
                         fTrueY = (float)((double)iCentreY - fTrueY);

                         int iFloorX = (int)(Math.Floor(fTrueX));
                         int iFloorY = (int)(Math.Floor(fTrueY));
                         int iCeilingX = (int)(Math.Ceiling(fTrueX));
                         int iCeilingY = (int)(Math.Ceiling(fTrueY));

                          //check bounds
                         if (iFloorX < 0 || iCeilingX < 0 || iFloorX >= image.Width ||
                             iCeilingX >= image.Width || iFloorY < 0 || iCeilingY < 0 ||
                             iFloorY >= image.Height || iCeilingY >= image.Height) continue;

                         //double fDeltaX = fTrueX - (double)iFloorX;
                         //double fDeltaY = fTrueY - (double)iFloorY;

                         //ColorFloatPixel clrTopLeft = image[iFloorX, iFloorY];
                         //ColorFloatPixel clrTopRight = image[iCeilingX, iFloorY];
                         //ColorFloatPixel clrBottomLeft = image[iFloorX, iCeilingY];
                         //ColorFloatPixel clrBottomRight = image[iCeilingX, iCeilingY];
                         
                         //check bounds
                         //if (iFloorX < -1 || iCeilingX < -1 || iFloorX >= (image.Width+1) ||
                         //    iCeilingX >= (image.Width+1) || iFloorY < -1 || iCeilingY < -1 ||
                         //    iFloorY >= (image.Height+1) || iCeilingY >= (image.Height+1)) continue;

                         //double fDeltaX = fTrueX - (double)iFloorX;
                         //double fDeltaY = fTrueY - (double)iFloorY;

                         //ColorFloatPixel clrTopLeft = imageCopyExtended[iFloorX+radius, iFloorY+radius];
                         //ColorFloatPixel clrTopRight = imageCopyExtended[iCeilingX + radius, iFloorY + radius];
                         //ColorFloatPixel clrBottomLeft = imageCopyExtended[iFloorX+radius, iCeilingY+radius];
                         //ColorFloatPixel clrBottomRight = imageCopyExtended[iCeilingX + radius, iCeilingY + radius];
                         // check bounds
                         //if (iFloorX < -2 && iCeilingX < -2 || iFloorX >= image.Width+2 ||
                         //   iCeilingX >= image.Width+2 || iFloorY < -2 || iCeilingY < -2 ||
                         //   iFloorY >= image.Height+2 || iCeilingY >= image.Height + 2) continue;
                         
                         
                         //if (iFloorX < 0) iFloorX = 0;
                         //if (iFloorY < 0) iFloorY = 0;
                         //if (iCeilingX < 0) iCeilingX = 0;
                         //if (iCeilingY < 0) iCeilingY = 0;
                         //if (iFloorX >= image.Width) iFloorX = image.Width-1;
                         //if (iFloorY >= image.Height) iFloorY = image.Height-1;
                         //if (iCeilingX >= image.Width) iCeilingX = image.Width-1;
                         //if (iCeilingY >= image.Height) iCeilingY = image.Height-1;



                         double fDeltaX = fTrueX - (double)iFloorX;
                         double fDeltaY = fTrueY - (double)iFloorY;
                         
                         ColorFloatPixel clrTopLeft = image[iFloorX , iFloorY];
                         ColorFloatPixel clrTopRight = image[iCeilingX , iFloorY];
                         ColorFloatPixel clrBottomLeft = image[iFloorX , iCeilingY];
                         ColorFloatPixel clrBottomRight = image[iCeilingX , iCeilingY];

                         // linearly interpolate horizontally between top neighbours
                         float fTopRed = (float)((1 - fDeltaX) * clrTopLeft.r + fDeltaX * clrTopRight.r);
                         float fTopGreen = (float)((1 - fDeltaX) * clrTopLeft.g + fDeltaX * clrTopRight.g);
                         float fTopBlue = (float)((1 - fDeltaX) * clrTopLeft.b + fDeltaX * clrTopRight.b);

                         // linearly interpolate horizontally between bottom neighbours
                         float fBottomRed = (float)((1 - fDeltaX) * clrBottomLeft.r + fDeltaX * clrBottomRight.r);
                         float fBottomGreen = (float)((1 - fDeltaX) * clrBottomLeft.g + fDeltaX * clrBottomRight.g);
                         float fBottomBlue = (float)((1 - fDeltaX) * clrBottomLeft.b + fDeltaX * clrBottomRight.b);

                         // linearly interpolate vertically between top and bottom interpolated results
                         float iRed = (int)(Math.Round((1 - fDeltaY) * fTopRed + fDeltaY * fBottomRed));
                         float iGreen = (int)(Math.Round((1 - fDeltaY) * fTopGreen + fDeltaY * fBottomGreen));
                         float iBlue = (int)(Math.Round((1 - fDeltaY) * fTopBlue + fDeltaY * fBottomBlue));

                         // make sure colour values are valid
                         if (iRed < 0) iRed = 0;
                         if (iRed > 255) iRed = 255;
                         if (iGreen < 0) iGreen = 0;
                         if (iGreen > 255) iGreen = 255;
                         if (iBlue < 0) iBlue = 0;
                         if (iBlue > 255) iBlue = 255;

                         ColorFloatPixel newPixelInterpolated = new ColorFloatPixel();
                         newPixelInterpolated.r = iRed;
                         newPixelInterpolated.g = iGreen;
                         newPixelInterpolated.b = iBlue;

                         bilinearInterpolatationForRotation[j, i] = newPixelInterpolated;



                     }
                 }
                 ImageIO.ImageToFile(bilinearInterpolatationForRotation, fileOutName);
             }
         }
         static void PrewittImage(GrayscaleFloatImage image,bool isHorizontal)
         {
             
             int radius = 1;
             GrayscaleFloatImage imageCopyExtended = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
             //mirrow edges
             //copy center
             //копируем в середину
             for (int k = 0; k < image.Width; k++)
             {
                 for (int l = 0; l < image.Height; l++)
                 {
                     imageCopyExtended[k + radius, l + radius] = image[k, l];
                 }
             }
             //копируем края
             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[k, radius - l - 1] = imageCopyExtended[k, l + radius];

                 }
             }

             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, l] = imageCopyExtended[k + radius, l];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                 }
             }
             /////////////////////////////////////////////////////
             //копируем уголки
             for (int k = 0; k < radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, radius - l - 1] = imageCopyExtended[k + radius, l + radius];

                 }
             }

             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                 }
             }
             

             
             for (int x = 1; x < imageCopyExtended.Width-1; x++)
             {
                 for (int y = 1; y < imageCopyExtended.Height-1; y++)
                 {

                     if (isHorizontal == true)
                     {
                         image[x-1 , y-1 ] = 128+ Math.Abs(imageCopyExtended[x - 1, y - 1] + imageCopyExtended[x - 1, y] + imageCopyExtended[x - 1, y + 1]) - 
                             (imageCopyExtended[x + 1, y - 1] + imageCopyExtended[x + 1, y] + imageCopyExtended[x + 1, y + 1]);
                  
                     }
                     else
                     {
                         image[x-1,y-1] =128+ (imageCopyExtended[x - 1, y - 1] + imageCopyExtended[x, y - 1] + imageCopyExtended[x + 1, y - 1])
                             - (imageCopyExtended[x - 1, y + 1] + imageCopyExtended[x, y + 1] + imageCopyExtended[x + 1, y + 1]);
                         
                     }

                 }
             }
             


         }
         static void SobelImage(GrayscaleFloatImage image, bool isHorizontal)
         {

             int radius = 1;
             GrayscaleFloatImage imageCopyExtended = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
             //mirrow edges
             //copy center
             //копируем в середину
             for (int k = 0; k < image.Width; k++)
             {
                 for (int l = 0; l < image.Height; l++)
                 {
                     imageCopyExtended[k + radius, l + radius] = image[k, l];
                 }
             }
             //копируем края
             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[k, radius - l - 1] = imageCopyExtended[k, l + radius];

                 }
             }

             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, l] = imageCopyExtended[k + radius, l];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                 }
             }
             /////////////////////////////////////////////////////
             //копируем уголки
             for (int k = 0; k < radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, radius - l - 1] = imageCopyExtended[k + radius, l + radius];

                 }
             }

             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                 }
             }

             for (int x = 1; x < imageCopyExtended.Width-1 ; x++)
             {
                 for (int y = 1; y < imageCopyExtended.Height-1 ; y++)
                 {

                     if (isHorizontal == true)
                     {
                         image[x-1 , y-1 ] = 128 + Math.Abs(imageCopyExtended[x - 1, y - 1] + 2*imageCopyExtended[x - 1, y] + imageCopyExtended[x - 1, y + 1]) -
                             (imageCopyExtended[x + 1, y - 1] + 2*imageCopyExtended[x + 1, y] + imageCopyExtended[x + 1, y + 1]);

                     }
                     else
                     {
                         image[x-1, y-1] = 128 + (imageCopyExtended[x - 1, y - 1] + 2*imageCopyExtended[x, y - 1] + imageCopyExtended[x + 1, y - 1])
                             - (imageCopyExtended[x - 1, y + 1] + 2*imageCopyExtended[x, y + 1] + imageCopyExtended[x + 1, y + 1]);

                     }

                 }
             }



         }

         static void RobertsImage(GrayscaleFloatImage image, bool isRigthDiagonal)
         {

             int radius = 1;
             GrayscaleFloatImage imageCopyExtended = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
             //mirrow edges
             //copy center
             //копируем в середину
             for (int k = 0; k < image.Width; k++)
             {
                 for (int l = 0; l < image.Height; l++)
                 {
                     imageCopyExtended[k + radius, l + radius] = image[k, l];
                 }
             }
             //копируем края
             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[k, radius - l - 1] = imageCopyExtended[k, l + radius];

                 }
             }

             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, l] = imageCopyExtended[k + radius, l];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                 }
             }
             /////////////////////////////////////////////////////
             //копируем уголки
             for (int k = 0; k < radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, radius - l - 1] = imageCopyExtended[k + radius, l + radius];

                 }
             }

             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                 }
             }


             for (int x = 0; x < image.Width ; x++)
             {
                 for (int y = 0; y < image.Height ; y++)
                 {

                     if (isRigthDiagonal == true)
                     {
                         image[x , y ] = 128 + imageCopyExtended[x + 1, y + 1] -imageCopyExtended[x , y ] ;

                     }
                     else
                     {
                         image[x, y] = 128 + imageCopyExtended[x , y + 1] - imageCopyExtended[x + 1, y ];

                     }

                 }
             }




         }


         static void MedianImage(ColorFloatImage image, int radius)
         {

             
             ColorFloatImage imageCopyExtended = new ColorFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
             //mirrow edges
             //copy center
             //копируем в середину
             for (int k = 0; k < image.Width; k++)
             {
                 for (int l = 0; l < image.Height; l++)
                 {
                     imageCopyExtended[k + radius, l + radius] = image[k, l];
                 }
             }
             //копируем края
             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[k, radius - l - 1] = imageCopyExtended[k, l + radius];

                 }
             }

             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, l] = imageCopyExtended[k + radius, l];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                 }
             }
             /////////////////////////////////////////////////////
             //копируем уголки
             for (int k = 0; k < radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, radius - l - 1] = imageCopyExtended[k + radius, l + radius];

                 }
             }

             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                 }
             }


             int index = 0;
             float[] SortArrR = new float[(2 * radius + 1) * (2 * radius + 1)];
             float[] SortArrG = new float[(2 * radius + 1) * (2 * radius + 1)];
             float[] SortArrB = new float[(2 * radius + 1) * (2 * radius + 1)];
             for (int x = radius; x < image.Width+radius; x++)
             {
                 for (int y = radius; y < image.Height+radius; y++)
                 {


                     index = 0;
                     for (int k = -radius; k < radius + 1; k++)
                     {
                         for (int n = -radius; n < radius + 1; n++)
                         {
                             SortArrR[index] = imageCopyExtended[x + k, y + n].r;
                             SortArrG[index] = imageCopyExtended[x + k, y + n].g;
                             SortArrB[index] = imageCopyExtended[x + k, y + n].b;
                             index++;
                         }

                     }
                     Array.Sort(SortArrR);
                     Array.Sort(SortArrG);
                     Array.Sort(SortArrB);
                     int midle = ((2 * radius + 1) * (2 * radius + 1) / 2);
                     ColorFloatPixel pixel = new ColorFloatPixel();
                     pixel.r = SortArrR[midle];
                     pixel.g = SortArrG[midle];
                     pixel.b = SortArrB[midle];
                     image[x - radius, y - radius] = pixel;
                     

                 }
             }




         }

         static double GaussFunction(float sigma, float i)
         {
             double result = (1 / (Math.Sqrt(2 * Math.PI) * sigma)) * Math.Exp((-i * i) / (2 * sigma * sigma));
             //////return (float)result;
             return result;
             ////            return Math.Exp (
             ////-(Math.Pow (i, 2)) / (2 * Math.Pow (sigma, 2)));
         }
         static double GaussFunctionDerivative(float sigma, float i)
         {
             //////double result = -(1 / (Math.Sqrt(2 * Math.PI) * Math.Pow( sigma,3)))* i * Math.Exp( (-i * i) / (2 * sigma * sigma));
             //////return (float)result;
             //double result = -(1 / (i * Math.Pow(sigma, 2))) *GaussFunction(sigma, i);
             //return result;
             return -GaussFunction(sigma, i) * (i) / Math.Pow(sigma, 2);
         }
         static void GaussImage(ColorFloatImage image, int radius,float sigma)
         {


             ColorFloatImage imageCopyExtended = new ColorFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
             //mirrow edges
             //copy center
             //копируем в середину
             for (int k = 0; k < image.Width; k++)
             {
                 for (int l = 0; l < image.Height; l++)
                 {
                     imageCopyExtended[k + radius, l + radius] = image[k, l];
                 }
             }
             //копируем края
             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[k, radius - l - 1] = imageCopyExtended[k, l + radius];

                 }
             }

             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, l] = imageCopyExtended[k + radius, l];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                 }
             }
             /////////////////////////////////////////////////////
             //копируем уголки
             for (int k = 0; k < radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, radius - l - 1] = imageCopyExtended[k + radius, l + radius];

                 }
             }

             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                 }
             }
             float[] Kernel = new float[radius * 2 + 1];
             float SumKernel = 0;
             //Console.WriteLine("GausKernel_without norm");
             for (int k = -radius; k < radius + 1; k++)
             {

                 Kernel[k + radius] = (float)GaussFunction(sigma, k);


                 SumKernel += Kernel[k + radius];
                 //Console.WriteLine(Kernel[k + radius] + "|");
             }
             for (int k = -radius; k < radius + 1; k++)
             {
                 Kernel[k + radius] /= SumKernel;

                 //Console.WriteLine(Kernel[k + radius] + "|");
             }
             string direction = "x";//no matter here(it make sence only for derivative)
             float sumR = 0;
             float sumG = 0;
             float sumB = 0;
             ////if (direction == "x")
             //{
             //    for (int j = radius; j < imageCopyExtended.Height - radius; j++)
             //    {
             //        for (int i = radius; i < imageCopyExtended.Width - radius; i++)
             //        {
             //            sumR = 0;
             //            sumG = 0;
             //            sumB = 0;
             //            for (int k = -radius; k < radius + 1; k++)
             //            {
             //                sumR += imageCopyExtended[i, j + k].r * Kernel[k + radius];
             //                sumG += imageCopyExtended[i, j + k].g * Kernel[k + radius];
             //                sumB += imageCopyExtended[i, j + k].b * Kernel[k + radius];
             //            }
             //            ColorFloatPixel pixel = new ColorFloatPixel();
             //            pixel.r = sumR;
             //            pixel.g = sumG;
             //            pixel.b = sumB;
             //            image[i - radius, j - radius] = pixel; 

             //        }
             //    }
             //}
             ////else //y
             //{
             //    for (int j = radius; j < imageCopyExtended.Height - radius; j++)
             //    {
             //        for (int i = radius; i < imageCopyExtended.Width - radius; i++)
             //        {
             //            sumR = 0;
             //            sumG = 0;
             //            sumB = 0;
             //            for (int k = -radius; k < radius + 1; k++)
             //            {
             //                sumR += imageCopyExtended[i, j + k].r * Kernel[k + radius];
             //                sumG += imageCopyExtended[i, j + k].g * Kernel[k + radius];
             //                sumB += imageCopyExtended[i, j + k].b * Kernel[k + radius];
             //            }

             //            ColorFloatPixel pixel = new ColorFloatPixel();
             //            pixel.r = sumR;
             //            pixel.g = sumG;
             //            pixel.b = sumB;
             //            image[i - radius, j - radius] = pixel; 
             //        }
             //    }
             //}
             float S = 0;
             float[,] GaussMatrix = new float[2 * radius + 1, 2 * radius + 1];
             for (int k = -radius; k < radius + 1; k++)
             {
                 for (int n = -radius; n < radius + 1; n++)
                 {
                     GaussMatrix[k + radius, n + radius] = (float)((1 / (2 * Math.PI * Math.Pow(radius, 2))) * Math.Exp(-((Math.Pow(k, 2) + Math.Pow(n, 2)) / (2 * Math.Pow(radius, 2)))));
                     //Console.Write(GaussMatrix[k+radius, n+radius]/(radius*2+1));
                     //Console.Write(" ");
                     S += GaussMatrix[k + radius, n + radius];
                     

                 }
                 //Console.WriteLine();

             }
             float S1 = 0;
             for (int k = -radius; k < radius + 1; k++)
             {
                 for (int n = -radius; n < radius + 1; n++)
                 {
                     GaussMatrix[k + radius, n + radius] = GaussMatrix[k + radius, n + radius] / S;
                     //Console.Write(GaussMatrix[k + radius, n + radius] / (radius * 2 + 1));
                     //Console.Write(" ");
                     S1 += GaussMatrix[k + radius, n + radius];
                     

                 }
                 //Console.WriteLine();

             }

             for (int j = radius; j < imageCopyExtended.Height - radius; j++)
             {
                 for (int i = radius; i < imageCopyExtended.Width - radius; i++)
                 {
                     sumR = 0;
                     sumG = 0;
                     sumB = 0;

                     for (int k = -radius; k < radius + 1; k++)
                     {
                         for (int n = -radius; n < radius + 1; n++)
                         {
                             sumR = sumR + GaussMatrix[k + radius, n + radius] * imageCopyExtended[i + k, j + n].r;
                             sumG = sumG + GaussMatrix[k + radius, n + radius] * imageCopyExtended[i + k, j + n].g;
                             sumB = sumB + GaussMatrix[k + radius, n + radius] * imageCopyExtended[i + k, j + n].b;

                         }

                     }


                     ColorFloatPixel pixel = new ColorFloatPixel();
                     pixel.r = sumR;
                     pixel.g = sumG;
                     pixel.b = sumB;
                     image[i - radius, j - radius] = pixel;
                 }
             }



         }

         static void GradientFromGaussImage(GrayscaleFloatImage image, int radius, float sigma)
         {


             GrayscaleFloatImage imageCopyExtended = new GrayscaleFloatImage(image.Width + 2 * radius, image.Height + 2 * radius);
             GrayscaleFloatImage derivativeX = new GrayscaleFloatImage(image.Width, image.Height);
             GrayscaleFloatImage derivativeY = new GrayscaleFloatImage(image.Width, image.Height);
             //for (int k = 0; k < image.Width; k++)
             //{
             //    for (int l = 0; l < image.Height; l++)
             //    {
             //        derivativeX[k, l] = image[k, l];
             //        derivativeY[k, l] = image[k, l];
             //    }
             //}


             //mirrow edges
             //copy center
             //копируем в середину
             for (int k = 0; k < image.Width; k++)
             {
                 for (int l = 0; l < image.Height; l++)
                 {
                     imageCopyExtended[k + radius, l + radius] = image[k, l];
                 }
             }
             //копируем края
             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[k, radius - l - 1] = imageCopyExtended[k, l + radius];

                 }
             }

             for (int k = radius; k < image.Width + radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[k, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k, l - radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, l] = imageCopyExtended[k + radius, l];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = radius; l < image.Height + radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, l] = imageCopyExtended[k - radius, l];

                 }
             }
             /////////////////////////////////////////////////////
             //копируем уголки
             for (int k = 0; k < radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, radius - l - 1] = imageCopyExtended[k + radius, l + radius];

                 }
             }

             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k - radius - 1, l - radius - 1];

                 }
             }
             for (int k = image.Width + radius; k < image.Width + 2 * radius; k++)
             {
                 for (int l = 0; l < radius; l++)
                 {

                     imageCopyExtended[image.Width + radius + image.Width + 2 * radius - k - 1, radius - l - 1] = imageCopyExtended[k - radius, l + radius];

                 }
             }
             for (int k = 0; k < radius; k++)
             {
                 for (int l = image.Height + radius; l < image.Height + 2 * radius; l++)
                 {

                     imageCopyExtended[radius - k - 1, image.Height + radius + image.Height + 2 * radius - l - 1] = imageCopyExtended[k + radius, l - radius];

                 }
             }
             float[] Kernel = new float[radius * 2 + 1];
             float SumKernel = 0;
             //Console.WriteLine("GausKernel_without norm");
             for (int k = -radius; k < radius + 1; k++)
             {

                 Kernel[k + radius] = (float)GaussFunctionDerivative(sigma, k);
                 SumKernel += Kernel[k + radius];
                 //Console.WriteLine(Kernel[k + radius] + "|");
             }
             float sum = 0;

             for (int j = radius; j < imageCopyExtended.Height - radius; j++)
             {
                 for (int i = radius; i < imageCopyExtended.Width - radius; i++)
                 {
                     sum = 0;
                     for (int k = -radius; k < radius + 1; k++)
                     {
                         sum += imageCopyExtended[i, j + k] * Kernel[k + radius];

                     }
                     derivativeX[i - radius, j - radius] = sum; //*100


                 }
             }

             //y

             for (int j = radius; j < imageCopyExtended.Height - radius; j++)
             {
                 for (int i = radius; i < imageCopyExtended.Width - radius; i++)
                 {
                     sum = 0;
                     for (int k = -radius; k < radius + 1; k++)
                     {
                         sum += imageCopyExtended[i + k, j] * Kernel[k + radius];
                     }
                     derivativeY[i - radius, j - radius] = sum; //*100

                 }
             }



             float maxX = derivativeX[0, 0];
             float minX = derivativeX[0, 0];
             float maxY = derivativeY[0, 0];
             float minY = derivativeY[0, 0];
             for (int j = 0; j < image.Height; j++)
             {
                 for (int i = 0; i < image.Width; i++)
                 {
                     if (derivativeX[i, j] > maxX) { maxX = derivativeX[i, j]; }
                     if (derivativeX[i, j] < minX) { minX = derivativeY[i, j]; }
                     if (derivativeY[i, j] > maxY) { maxY = derivativeY[i, j]; }
                     if (derivativeY[i, j] < minY) { minY = derivativeY[i, j]; }
                 }
             }

             for (int j = 0; j < image.Height; j++)
             {
                 for (int i = 0; i < image.Width; i++)
                 {
                     derivativeX[i, j] = (derivativeX[i, j] / (maxX)) * 256;
                     derivativeY[i, j] = (derivativeY[i, j] / (maxY)) * 256;
                 }
             }



             for (int x = 0; x < image.Width; x++)
             {
                 for (int y = 0; y < image.Height; y++)
                 {
                     image[x, y] = (float)Math.Sqrt(derivativeX[x, y] * derivativeX[x, y] + derivativeY[x, y] * derivativeY[x, y]);
                 }
             }

           //  ImageIO.ImageToFile(derivativeX, "InputImage/DerivativeX.bmp");
          //   ImageIO.ImageToFile(derivativeY, "InputImage/DerivativeY.bmp");


         }
        static void eqhist(GrayscaleFloatImage image)
         {
            float[] histogram = new float[256];
            //hist
            for (int i=0;i<image.Width;i++)
            {
                for (int j=0;j<image.Height;j++)
                {
                    histogram[(int)image[i, j]]++;  
                }
            }
             // cumulative distribution function
            float[] cdf = new float[256];
            cdf[0] = histogram[0];
            int L = 256;
            for (int i=1; i<cdf.Length; i++ )
            {
                cdf[i] = cdf[i - 1] + histogram[i];
            }
            float minNonZero = 0;
            for (int i = 0; i < cdf.Length; i++)
            {
                if (cdf[i]>0)
                {
                    minNonZero = cdf[i];
                    break;
                }
            }
            for (int i = 0; i < cdf.Length; i++)
            {
                histogram[i] = (float)Math.Round((double)((cdf[i] - minNonZero) * (L - 2) / (image.Height * image.Width - minNonZero))) + 1;
            }
            for (int i = 0; i < image.Width; i++)
            {
                for (int j = 0; j < image.Height; j++)
                {
                    image[i, j] = histogram[(int)image[i, j]];
                }
            }

         }


        static void dcci(GrayscaleFloatImage image, int radius, string fileOutName)
        {
            radius = 2;
            int radius_more = 4;
            GrayscaleFloatImage imageCopyExtended = new GrayscaleFloatImage(image.Width + 2 * radius_more, image.Height + 2 * radius_more);
            GrayscaleFloatImage imageCopyExtended_DCIM = new GrayscaleFloatImage((imageCopyExtended.Width * radius), (imageCopyExtended.Height * radius));//((image.Width * radius), (image.Height * radius));
            GrayscaleFloatImage imageResult = new GrayscaleFloatImage((image.Width * radius), (image.Height * radius));
            //mirrow edges
            //copy center
            //копируем в середину
            for (int k = 0; k < image.Width; k++)
            {
                for (int l = 0; l < image.Height; l++)
                {
                    imageCopyExtended[k + radius_more, l + radius_more] = image[k, l];
                }
            }
            //копируем края
            for (int k = radius_more; k < image.Width + radius_more; k++)
            {
                for (int l = 0; l < radius_more; l++)
                {

                    imageCopyExtended[k, radius_more - l - 1] = imageCopyExtended[k, l + radius_more];

                }
            }

            for (int k = radius_more; k < image.Width + radius_more; k++)
            {
                for (int l = image.Height + radius_more; l < image.Height + 2 * radius_more; l++)
                {

                    imageCopyExtended[k, image.Height + radius_more + image.Height + 2 * radius_more - l - 1] = imageCopyExtended[k, l - radius_more];

                }
            }
            for (int k = 0; k < radius_more; k++)
            {
                for (int l = radius_more; l < image.Height + radius_more; l++)
                {

                    imageCopyExtended[radius_more - k - 1, l] = imageCopyExtended[k + radius_more, l];

                }
            }
            for (int k = image.Width + radius_more; k < image.Width + 2 * radius_more; k++)
            {
                for (int l = radius_more; l < image.Height + radius_more; l++)
                {

                    imageCopyExtended[image.Width + radius_more + image.Width + 2 * radius_more - k - 1, l] = imageCopyExtended[k - radius_more, l];

                }
            }
            /////////////////////////////////////////////////////
            //копируем уголки
            for (int k = 0; k < radius_more; k++)
            {
                for (int l = 0; l < radius_more; l++)
                {

                    imageCopyExtended[radius_more - k - 1, radius_more - l - 1] = imageCopyExtended[k + radius_more, l + radius_more];

                }
            }

            for (int k = image.Width + radius_more; k < image.Width + 2 * radius_more; k++)
            {
                for (int l = image.Height + radius_more; l < image.Height + 2 * radius_more; l++)
                {

                    imageCopyExtended[image.Width + radius_more + image.Width + 2 * radius_more - k - 1, image.Height + radius_more + image.Height + 2 * radius_more - l - 1] = imageCopyExtended[k - radius_more - 1, l - radius_more - 1];

                }
            }
            for (int k = image.Width + radius_more; k < image.Width + 2 * radius_more; k++)
            {
                for (int l = 0; l < radius_more; l++)
                {

                    imageCopyExtended[image.Width + radius_more + image.Width + 2 * radius_more - k - 1, radius_more - l - 1] = imageCopyExtended[k - radius_more, l + radius_more];

                }
            }
            for (int k = 0; k < radius_more; k++)
            {
                for (int l = image.Height + radius_more; l < image.Height + 2 * radius_more; l++)
                {

                    imageCopyExtended[radius_more - k - 1, image.Height + radius_more + image.Height + 2 * radius_more - l - 1] = imageCopyExtended[k + radius_more, l - radius_more];

                }
            }
            
            
            for (int i=0;i<imageCopyExtended.Width;i++)
            {
                for (int j=0;j<imageCopyExtended.Height;j++)
                {
                    imageCopyExtended_DCIM[i * radius, j * radius] = imageCopyExtended[i, j];
                }
            }
            /////////////////////////////////////////////////
            for (int i2 = radius; i2 < imageCopyExtended.Width-radius;i2++ )
            {
                for (int j2=radius;j2<imageCopyExtended.Height-radius;j2++)
                {
                    int i = i2*radius + 1;
                    int j = j2*radius + 1;
                    double G1 = 0;
                    double G2 = 0;


                    G1 = Math.Abs(imageCopyExtended_DCIM[i + 3, j - 3] - imageCopyExtended_DCIM[i + 1, j - 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1] - imageCopyExtended_DCIM[i + 1, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1] - imageCopyExtended_DCIM[i + 1, j + 3])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 3] - imageCopyExtended_DCIM[i - 1, j - 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1] - imageCopyExtended_DCIM[i - 1, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1] - imageCopyExtended_DCIM[i - 1, j + 3])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 3] - imageCopyExtended_DCIM[i - 3, j - 1])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1] - imageCopyExtended_DCIM[i - 3, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1] - imageCopyExtended_DCIM[i - 3, j + 3]);

                    G2 = Math.Abs(imageCopyExtended_DCIM[i + 3, j + 3] - imageCopyExtended_DCIM[i + 1, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1] - imageCopyExtended_DCIM[i + 1, j - 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1] - imageCopyExtended_DCIM[i + 1, j - 3])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 3] - imageCopyExtended_DCIM[i - 1, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1] - imageCopyExtended_DCIM[i - 1, j - 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1] - imageCopyExtended_DCIM[i - 1, j - 3])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 3] - imageCopyExtended_DCIM[i - 3, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1] - imageCopyExtended_DCIM[i - 3, j - 1])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1] - imageCopyExtended_DCIM[i - 3, j - 3]);

                    double T = 1.15;
                    double p=0, p1=0, p2=0;

                    //if ((1 + G1) / (1 + G2) > T)
                    //{
                    //    p = p1;
                    //    p = (-1 * imageCopyExtended[i2, j2] + 9 * imageCopyExtended[i2 + 1, j2 + 1] + 9 * imageCopyExtended[i2 + 2, j2 + 2] - 1 * imageCopyExtended[i2 + 3, j2 + 3]) / 16;
                    //}
                    //else if ((1 + G2) / (1 + G1) > T)
                    //{
                    //    p = p2;
                    //    p = (-1 * imageCopyExtended[i2 + 3, j2] + 9 * imageCopyExtended[i2 + 2, j2 + 1] + 9 * imageCopyExtended[i2 + 1, j2 + 2] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //}
                    //else
                    //{
                    //    p1 = (-1 * imageCopyExtended[i2, j2] + 9 * imageCopyExtended[i2 + 1, j2 + 1] + 9 * imageCopyExtended[i2 + 2, j2 + 2] - 1 * imageCopyExtended[i2 + 3, j2 + 3]) / 16;
                    //    p2 = (-1 * imageCopyExtended[i2 + 3, j2] + 9 * imageCopyExtended[i2 + 2, j2 + 1] + 9 * imageCopyExtended[i2 + 1, j2 + 2] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //    double w1, w2;
                    //    int k = 5;
                    //    w1 = 1 / (1 + Math.Pow(G1, k));
                    //    w2 = 1 / (1 + Math.Pow(G2, k));
                    //    p = (w1 * p1 + w2 * p2) / (w1 + w2);
                    //}

                    if ((1 + G1) / (1 + G2) > T)
                    {
                        p = p1;
                        p = (-1 * imageCopyExtended_DCIM[i-3, j-3] + 9 * imageCopyExtended_DCIM[i - 1, j - 1] + 9 * imageCopyExtended_DCIM[i + 1, j + 1] - 1 * imageCopyExtended_DCIM[i + 3, j + 3]) / 16;
                    }
                    else if ((1 + G2) / (1 + G1) > T)
                    {
                        p = p2;
                        p = (-1 * imageCopyExtended_DCIM[i + 3, j - 3] + 9 * imageCopyExtended_DCIM[i + 1, j - 1] + 9 * imageCopyExtended_DCIM[i - 1, j + 1] - 1 * imageCopyExtended_DCIM[i - 3, j + 3]) / 16;
                    }
                    else
                    {
                        p1 = (-1 * imageCopyExtended_DCIM[i - 3, j - 3] + 9 * imageCopyExtended_DCIM[i - 1, j - 1] + 9 * imageCopyExtended_DCIM[i + 1, j + 1] - 1 * imageCopyExtended_DCIM[i + 3, j + 3]) / 16;
                        p2 = (-1 * imageCopyExtended_DCIM[i + 3, j - 3] + 9 * imageCopyExtended_DCIM[i + 1, j - 1] + 9 * imageCopyExtended_DCIM[i - 1, j + 1] - 1 * imageCopyExtended_DCIM[i - 3, j + 3]) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1, k));
                        w2 = 1 / (1 + Math.Pow(G2, k));
                        p = (w1 * p1 + w2 * p2) / (w1 + w2);
                    }


                    imageCopyExtended_DCIM[i, j] = (float)p;

                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////step 2
                    ///////////////////////////////////////////////////////////////////////////////////
                    //////G1 = Math.Abs(imageCopyExtended_DCIM[i + 1, j] - imageCopyExtended_DCIM[i + 1, j + 2])
                    //////    + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2] - imageCopyExtended_DCIM[i + 1, j])
                    //////    + Math.Abs(imageCopyExtended_DCIM[i - 1, j] - imageCopyExtended_DCIM[i - 1, j + 2])
                    //////    + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2] - imageCopyExtended_DCIM[i + 1, j])
                    //////    //horizontal
                    //////    + Math.Abs(imageCopyExtended_DCIM[i, j - 1] - imageCopyExtended_DCIM[i, j + 1])
                    //////    + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1] - imageCopyExtended_DCIM[i + 2, j + 1])
                    //////    + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1] - imageCopyExtended_DCIM[i - 2, j + 1]);

                    //////G2 = Math.Abs(imageCopyExtended_DCIM[i , j + 1] - imageCopyExtended_DCIM[i + 2, j + 1])
                    //////  + Math.Abs(imageCopyExtended_DCIM[i , j - 1] - imageCopyExtended_DCIM[i + 2, j - 1])
                    //////  + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1] - imageCopyExtended_DCIM[i , j + 1])
                    //////  + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1] - imageCopyExtended_DCIM[i , j - 1])
                    //////    //vertical
                    //////  + Math.Abs(imageCopyExtended_DCIM[i - 1, j] - imageCopyExtended_DCIM[i + 1, j])
                    //////  + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2] - imageCopyExtended_DCIM[i + 1, j + 2])
                    //////  + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2] - imageCopyExtended_DCIM[i + 1, j - 2]);

                    //////if ((1 + G1) / (1 + G2) > T)
                    //////{
                    //////    p = p1;
                    //////    // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                    //////    p = (-1 * imageCopyExtended[i2, j2 - 3] + 9 * imageCopyExtended[i2, j2 - 1] + 9 * imageCopyExtended[i2, j2 + 1] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //////}
                    //////else if ((1 + G2) / (1 + G1) > T)
                    //////{
                    //////    p = p2;
                    //////    //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                    //////    p = (-1 * imageCopyExtended[i2 - 3, j2] + 9 * imageCopyExtended[i2 - 1, j2] + 9 * imageCopyExtended[i2 + 1, j2] - 1 * imageCopyExtended[i2 + 3, j2 ]) / 16;
                    //////}
                    //////else
                    //////{
                    //////    p1 = (-1 * imageCopyExtended[i2, j2 - 3] + 9 * imageCopyExtended[i2, j2 - 1] + 9 * imageCopyExtended[i2, j2 + 1] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //////    p2 = (-1 * imageCopyExtended[i2 - 3, j2] + 9 * imageCopyExtended[i2 - 1, j2] + 9 * imageCopyExtended[i2 + 1, j2] - 1 * imageCopyExtended[i2 + 3, j2]) / 16;
                    //////    double w1, w2;
                    //////    int k = 5;
                    //////    w1 = 1 / (1 + Math.Pow(G1, k));
                    //////    w2 = 1 / (1 + Math.Pow(G2, k));
                    //////    p = (w1 * p1 + w2 * p2) / (w1 + w2);
                    //////}


                }
            }

            for (int i2 = radius+1; i2 < imageCopyExtended.Width - radius - 1; i2++)
            {
                for (int j2 = radius+1; j2 < imageCopyExtended.Height - radius - 1; j2++)
                {
                    int i = i2*radius+1;
                    int j = j2*radius;
                    double G1 = 0;
                    double G2 = 0;
                    double T = 1.15;
                    double p = 0, p1 = 0, p2 = 0;

                    G1 = Math.Abs(imageCopyExtended_DCIM[i + 1, j] - imageCopyExtended_DCIM[i + 1, j + 2])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2] - imageCopyExtended_DCIM[i + 1, j])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j] - imageCopyExtended_DCIM[i - 1, j + 2])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2] - imageCopyExtended_DCIM[i + 1, j])
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1] - imageCopyExtended_DCIM[i, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1] - imageCopyExtended_DCIM[i + 2, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1] - imageCopyExtended_DCIM[i - 2, j + 1]);

                    G2 = Math.Abs(imageCopyExtended_DCIM[i, j + 1] - imageCopyExtended_DCIM[i + 2, j + 1])
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1] - imageCopyExtended_DCIM[i + 2, j - 1])
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1] - imageCopyExtended_DCIM[i, j + 1])
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1] - imageCopyExtended_DCIM[i, j - 1])
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j] - imageCopyExtended_DCIM[i + 1, j])
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2] - imageCopyExtended_DCIM[i + 1, j + 2])
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2] - imageCopyExtended_DCIM[i + 1, j - 2]);

                    //if ((1 + G1) / (1 + G2) > T)
                    //{
                    //    p = p1;
                    //    // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                    //    p = (-1 * imageCopyExtended[i2, j2 - 3] + 9 * imageCopyExtended[i2, j2 - 1] + 9 * imageCopyExtended[i2, j2 + 1] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //}
                    //else if ((1 + G2) / (1 + G1) > T)
                    //{
                    //    p = p2;
                    //    //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                    //    p = (-1 * imageCopyExtended[i2 - 3, j2] + 9 * imageCopyExtended[i2 - 1, j2] + 9 * imageCopyExtended[i2 + 1, j2] - 1 * imageCopyExtended[i2 + 3, j2]) / 16;
                    //}
                    //else
                    //{
                    //    p1 = (-1 * imageCopyExtended[i2, j2 - 3] + 9 * imageCopyExtended[i2, j2 - 1] + 9 * imageCopyExtended[i2, j2 + 1] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //    p2 = (-1 * imageCopyExtended[i2 - 3, j2] + 9 * imageCopyExtended[i2 - 1, j2] + 9 * imageCopyExtended[i2 + 1, j2] - 1 * imageCopyExtended[i2 + 3, j2]) / 16;
                    //    double w1, w2;
                    //    int k = 5;
                    //    w1 = 1 / (1 + Math.Pow(G1, k));
                    //    w2 = 1 / (1 + Math.Pow(G2, k));
                    //    p = (w1 * p1 + w2 * p2) / (w1 + w2);
                    //}
                    if ((1 + G1) / (1 + G2) > T)
                    {
                        p = p1;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p = (-1 * imageCopyExtended_DCIM[i, j - 3] + 9 * imageCopyExtended_DCIM[i, j - 1] + 9 * imageCopyExtended_DCIM[i, j + 1] - 1 * imageCopyExtended_DCIM[i, j + 3]) / 16;
                    }
                    else if ((1 + G2) / (1 + G1) > T)
                    {
                        p = p2;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p = (-1 * imageCopyExtended_DCIM[i - 3, j] + 9 * imageCopyExtended_DCIM[i - 1, j] + 9 * imageCopyExtended_DCIM[i + 1, j] - 1 * imageCopyExtended_DCIM[i + 3, j]) / 16;
                    }
                    else
                    {
                        p1 = (-1 * imageCopyExtended_DCIM[i, j - 3] + 9 * imageCopyExtended_DCIM[i, j - 1] + 9 * imageCopyExtended_DCIM[i, j + 1] - 1 * imageCopyExtended_DCIM[i, j + 3]) / 16;
                        p2 = (-1 * imageCopyExtended_DCIM[i - 3, j] + 9 * imageCopyExtended_DCIM[i - 1, j] + 9 * imageCopyExtended_DCIM[i + 1, j] - 1 * imageCopyExtended_DCIM[i + 3, j]) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1, k));
                        w2 = 1 / (1 + Math.Pow(G2, k));
                        p = (w1 * p1 + w2 * p2) / (w1 + w2);
                    }

                    imageCopyExtended_DCIM[i, j] = (float)p;


                }
            }

            for (int i2 = radius+1; i2 < imageCopyExtended.Width - radius - 1; i2 ++)
            {
                for (int j2 = radius+1; j2 < imageCopyExtended.Height - radius - 1; j2 ++)
                {
                    int i = i2*radius;
                    int j = j2*radius+1;
                    double G1 = 0;
                    double G2 = 0;
                    double T = 1.15;
                    double p = 0, p1 = 0, p2 = 0;

                    G1 = Math.Abs(imageCopyExtended_DCIM[i + 1, j] - imageCopyExtended_DCIM[i + 1, j + 2])
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2] - imageCopyExtended_DCIM[i + 1, j])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j] - imageCopyExtended_DCIM[i - 1, j + 2])
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2] - imageCopyExtended_DCIM[i + 1, j])
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1] - imageCopyExtended_DCIM[i, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1] - imageCopyExtended_DCIM[i + 2, j + 1])
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1] - imageCopyExtended_DCIM[i - 2, j + 1]);


                    G2 = Math.Abs(imageCopyExtended_DCIM[i, j + 1] - imageCopyExtended_DCIM[i + 2, j + 1])
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1] - imageCopyExtended_DCIM[i + 2, j - 1])
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1] - imageCopyExtended_DCIM[i, j + 1])
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1] - imageCopyExtended_DCIM[i, j - 1])
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j] - imageCopyExtended_DCIM[i + 1, j])
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2] - imageCopyExtended_DCIM[i + 1, j + 2])
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2] - imageCopyExtended_DCIM[i + 1, j - 2]);


                    //if ((1 + G1) / (1 + G2) > T)
                    //{
                    //    p = p1;
                    //    // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                    //    p = (-1 * imageCopyExtended[i2, j2 - 3] + 9 * imageCopyExtended[i2, j2 - 1] + 9 * imageCopyExtended[i2, j2 + 1] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //}
                    //else if ((1 + G2) / (1 + G1) > T)
                    //{
                    //    p = p2;
                    //    //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                    //    p = (-1 * imageCopyExtended[i2 - 3, j2] + 9 * imageCopyExtended[i2 - 1, j2] + 9 * imageCopyExtended[i2 + 1, j2] - 1 * imageCopyExtended[i2 + 3, j2]) / 16;
                    //}
                    //else
                    //{
                    //    p1 = (-1 * imageCopyExtended[i2, j2 - 3] + 9 * imageCopyExtended[i2, j2 - 1] + 9 * imageCopyExtended[i2, j2 + 1] - 1 * imageCopyExtended[i2, j2 + 3]) / 16;
                    //    p2 = (-1 * imageCopyExtended[i2 - 3, j2] + 9 * imageCopyExtended[i2 - 1, j2] + 9 * imageCopyExtended[i2 + 1, j2] - 1 * imageCopyExtended[i2 + 3, j2]) / 16;
                    //    double w1, w2;
                    //    int k = 5;
                    //    w1 = 1 / (1 + Math.Pow(G1, k));
                    //    w2 = 1 / (1 + Math.Pow(G2, k));
                    //    p = (w1 * p1 + w2 * p2) / (w1 + w2);
                    //}

                    if ((1 + G1) / (1 + G2) > T)
                    {
                        p = p1;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p = (-1 * imageCopyExtended_DCIM[i, j - 3] + 9 * imageCopyExtended_DCIM[i, j - 1] + 9 * imageCopyExtended_DCIM[i, j + 1] - 1 * imageCopyExtended_DCIM[i, j + 3]) / 16;
                    }
                    else if ((1 + G2) / (1 + G1) > T)
                    {
                        p = p2;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p = (-1 * imageCopyExtended_DCIM[i - 3, j] + 9 * imageCopyExtended_DCIM[i - 1, j] + 9 * imageCopyExtended_DCIM[i + 1, j] - 1 * imageCopyExtended_DCIM[i + 3, j]) / 16;
                    }
                    else
                    {
                        p1 = (-1 * imageCopyExtended_DCIM[i, j - 3] + 9 * imageCopyExtended_DCIM[i, j - 1] + 9 * imageCopyExtended_DCIM[i, j + 1] - 1 * imageCopyExtended_DCIM[i, j + 3]) / 16;
                        p2 = (-1 * imageCopyExtended_DCIM[i - 3, j] + 9 * imageCopyExtended_DCIM[i - 1, j] + 9 * imageCopyExtended_DCIM[i + 1, j] - 1 * imageCopyExtended_DCIM[i + 3, j]) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1, k));
                        w2 = 1 / (1 + Math.Pow(G2, k));
                        p = (w1 * p1 + w2 * p2) / (w1 + w2);
                    }

                    imageCopyExtended_DCIM[i, j] = (float)p;

                    //imageCopyExtended_DCIM[i, j] = 255;


                }
            }


            for (int i = 0; i < imageResult.Width; i++)
            {
                for (int j = 0; j < imageResult.Height; j++)
                {
                    imageResult[i,j]=imageCopyExtended_DCIM[i+2*radius_more, j+2*radius_more];
                    
                }
            }
                ImageIO.ImageToFile(imageResult, fileOutName);
        }


        static void dcci(ColorFloatImage image, int radius, string fileOutName)
        {
            radius = 2;
            int radius_more = 4;
            ColorFloatImage imageCopyExtended = new ColorFloatImage(image.Width + 2 * radius_more, image.Height + 2 * radius_more);
            ColorFloatImage imageCopyExtended_DCIM = new ColorFloatImage((imageCopyExtended.Width * radius), (imageCopyExtended.Height * radius));//((image.Width * radius), (image.Height * radius));
            ColorFloatImage imageResult = new ColorFloatImage((image.Width * radius), (image.Height * radius));
            //mirrow edges
            //copy center
            //копируем в середину
            for (int k = 0; k < image.Width; k++)
            {
                for (int l = 0; l < image.Height; l++)
                {
                    imageCopyExtended[k + radius_more, l + radius_more] = image[k, l];
                }
            }
            //копируем края
            for (int k = radius_more; k < image.Width + radius_more; k++)
            {
                for (int l = 0; l < radius_more; l++)
                {

                    imageCopyExtended[k, radius_more - l - 1] = imageCopyExtended[k, l + radius_more];

                }
            }

            for (int k = radius_more; k < image.Width + radius_more; k++)
            {
                for (int l = image.Height + radius_more; l < image.Height + 2 * radius_more; l++)
                {

                    imageCopyExtended[k, image.Height + radius_more + image.Height + 2 * radius_more - l - 1] = imageCopyExtended[k, l - radius_more];

                }
            }
            for (int k = 0; k < radius_more; k++)
            {
                for (int l = radius_more; l < image.Height + radius_more; l++)
                {

                    imageCopyExtended[radius_more - k - 1, l] = imageCopyExtended[k + radius_more, l];

                }
            }
            for (int k = image.Width + radius_more; k < image.Width + 2 * radius_more; k++)
            {
                for (int l = radius_more; l < image.Height + radius_more; l++)
                {

                    imageCopyExtended[image.Width + radius_more + image.Width + 2 * radius_more - k - 1, l] = imageCopyExtended[k - radius_more, l];

                }
            }
            /////////////////////////////////////////////////////
            //копируем уголки
            for (int k = 0; k < radius_more; k++)
            {
                for (int l = 0; l < radius_more; l++)
                {

                    imageCopyExtended[radius_more - k - 1, radius_more - l - 1] = imageCopyExtended[k + radius_more, l + radius_more];

                }
            }

            for (int k = image.Width + radius_more; k < image.Width + 2 * radius_more; k++)
            {
                for (int l = image.Height + radius_more; l < image.Height + 2 * radius_more; l++)
                {

                    imageCopyExtended[image.Width + radius_more + image.Width + 2 * radius_more - k - 1, image.Height + radius_more + image.Height + 2 * radius_more - l - 1] = imageCopyExtended[k - radius_more - 1, l - radius_more - 1];

                }
            }
            for (int k = image.Width + radius_more; k < image.Width + 2 * radius_more; k++)
            {
                for (int l = 0; l < radius_more; l++)
                {

                    imageCopyExtended[image.Width + radius_more + image.Width + 2 * radius_more - k - 1, radius_more - l - 1] = imageCopyExtended[k - radius_more, l + radius_more];

                }
            }
            for (int k = 0; k < radius_more; k++)
            {
                for (int l = image.Height + radius_more; l < image.Height + 2 * radius_more; l++)
                {

                    imageCopyExtended[radius_more - k - 1, image.Height + radius_more + image.Height + 2 * radius_more - l - 1] = imageCopyExtended[k + radius_more, l - radius_more];

                }
            }


            for (int i = 0; i < imageCopyExtended.Width; i++)
            {
                for (int j = 0; j < imageCopyExtended.Height; j++)
                {
                    imageCopyExtended_DCIM[i * radius, j * radius] = imageCopyExtended[i, j];
                }
            }
            /////////////////////////////////////////////////
            for (int i2 = radius; i2 < imageCopyExtended.Width - radius; i2++)
            {
                for (int j2 = radius; j2 < imageCopyExtended.Height - radius; j2++)
                {
                    int i = i2 * radius + 1;
                    int j = j2 * radius + 1;
                    double G1_r = 0;
                    double G2_r = 0;
                    double G1_g = 0;
                    double G2_g = 0;
                    double G1_b = 0;
                    double G2_b = 0;


                    G1_r = Math.Abs(imageCopyExtended_DCIM[i + 3, j - 3].r - imageCopyExtended_DCIM[i + 1, j - 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1].r - imageCopyExtended_DCIM[i + 1, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1].r - imageCopyExtended_DCIM[i + 1, j + 3].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 3].r - imageCopyExtended_DCIM[i - 1, j - 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1].r - imageCopyExtended_DCIM[i - 1, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1].r - imageCopyExtended_DCIM[i - 1, j + 3].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 3].r - imageCopyExtended_DCIM[i - 3, j - 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1].r - imageCopyExtended_DCIM[i - 3, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1].r - imageCopyExtended_DCIM[i - 3, j + 3].r);

                    G2_r = Math.Abs(imageCopyExtended_DCIM[i + 3, j + 3].r - imageCopyExtended_DCIM[i + 1, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1].r - imageCopyExtended_DCIM[i + 1, j - 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1].r - imageCopyExtended_DCIM[i + 1, j - 3].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 3].r - imageCopyExtended_DCIM[i - 1, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1].r - imageCopyExtended_DCIM[i - 1, j - 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1].r - imageCopyExtended_DCIM[i - 1, j - 3].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 3].r - imageCopyExtended_DCIM[i - 3, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1].r - imageCopyExtended_DCIM[i - 3, j - 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1].r - imageCopyExtended_DCIM[i - 3, j - 3].r);

                    G1_g = Math.Abs(imageCopyExtended_DCIM[i + 3, j - 3].g - imageCopyExtended_DCIM[i + 1, j - 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1].g - imageCopyExtended_DCIM[i + 1, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1].g - imageCopyExtended_DCIM[i + 1, j + 3].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 3].g - imageCopyExtended_DCIM[i - 1, j - 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1].g - imageCopyExtended_DCIM[i - 1, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1].g - imageCopyExtended_DCIM[i - 1, j + 3].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 3].g - imageCopyExtended_DCIM[i - 3, j - 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1].g - imageCopyExtended_DCIM[i - 3, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1].g - imageCopyExtended_DCIM[i - 3, j + 3].g);

                    G2_g = Math.Abs(imageCopyExtended_DCIM[i + 3, j + 3].g - imageCopyExtended_DCIM[i + 1, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1].g - imageCopyExtended_DCIM[i + 1, j - 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1].g - imageCopyExtended_DCIM[i + 1, j - 3].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 3].g - imageCopyExtended_DCIM[i - 1, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1].g - imageCopyExtended_DCIM[i - 1, j - 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1].g - imageCopyExtended_DCIM[i - 1, j - 3].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 3].g - imageCopyExtended_DCIM[i - 3, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1].g - imageCopyExtended_DCIM[i - 3, j - 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1].g - imageCopyExtended_DCIM[i - 3, j - 3].g);

                    G1_b = Math.Abs(imageCopyExtended_DCIM[i + 3, j - 3].b - imageCopyExtended_DCIM[i + 1, j - 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1].b - imageCopyExtended_DCIM[i + 1, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1].b - imageCopyExtended_DCIM[i + 1, j + 3].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 3].b - imageCopyExtended_DCIM[i - 1, j - 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1].b - imageCopyExtended_DCIM[i - 1, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1].b - imageCopyExtended_DCIM[i - 1, j + 3].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 3].b - imageCopyExtended_DCIM[i - 3, j - 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1].b - imageCopyExtended_DCIM[i - 3, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1].b - imageCopyExtended_DCIM[i - 3, j + 3].b);

                    G2_b = Math.Abs(imageCopyExtended_DCIM[i + 3, j + 3].b - imageCopyExtended_DCIM[i + 1, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j + 1].b - imageCopyExtended_DCIM[i + 1, j - 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 3, j - 1].b - imageCopyExtended_DCIM[i + 1, j - 3].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 3].b - imageCopyExtended_DCIM[i - 1, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j + 1].b - imageCopyExtended_DCIM[i - 1, j - 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 1].b - imageCopyExtended_DCIM[i - 1, j - 3].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 3].b - imageCopyExtended_DCIM[i - 3, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 1].b - imageCopyExtended_DCIM[i - 3, j - 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 1].b - imageCopyExtended_DCIM[i - 3, j - 3].b);

                    double T = 1.15;
                    double p_r = 0, p1_r = 0, p2_r = 0;
                    double p_g = 0, p1_g = 0, p2_g = 0;
                    double p_b = 0, p1_b = 0, p2_b = 0;


                    if ((1 + G1_r) / (1 + G2_r) > T)
                    {
                        p_r = p1_r;
                        p_r = (-1 * imageCopyExtended_DCIM[i - 3, j - 3].r + 9 * imageCopyExtended_DCIM[i - 1, j - 1].r + 9 * imageCopyExtended_DCIM[i + 1, j + 1].r - 1 * imageCopyExtended_DCIM[i + 3, j + 3].r) / 16;
                    }
                    else if ((1 + G2_r) / (1 + G1_r) > T)
                    {
                        p_r = p2_r;
                        p_r = (-1 * imageCopyExtended_DCIM[i + 3, j - 3].r + 9 * imageCopyExtended_DCIM[i + 1, j - 1].r + 9 * imageCopyExtended_DCIM[i - 1, j + 1].r - 1 * imageCopyExtended_DCIM[i - 3, j + 3].r) / 16;
                    }
                    else
                    {
                        p1_r = (-1 * imageCopyExtended_DCIM[i - 3, j - 3].r + 9 * imageCopyExtended_DCIM[i - 1, j - 1].r + 9 * imageCopyExtended_DCIM[i + 1, j + 1].r - 1 * imageCopyExtended_DCIM[i + 3, j + 3].r) / 16;
                        p2_r = (-1 * imageCopyExtended_DCIM[i + 3, j - 3].r + 9 * imageCopyExtended_DCIM[i + 1, j - 1].r + 9 * imageCopyExtended_DCIM[i - 1, j + 1].r - 1 * imageCopyExtended_DCIM[i - 3, j + 3].r) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_r, k));
                        w2 = 1 / (1 + Math.Pow(G2_r, k));
                        p_r = (w1 * p1_r + w2 * p2_r) / (w1 + w2);
                    }



                    if ((1 + G1_g) / (1 + G2_g) > T)
                    {
                        p_g = p1_g;
                        p_g = (-1 * imageCopyExtended_DCIM[i - 3, j - 3].g + 9 * imageCopyExtended_DCIM[i - 1, j - 1].g + 9 * imageCopyExtended_DCIM[i + 1, j + 1].g - 1 * imageCopyExtended_DCIM[i + 3, j + 3].g) / 16;
                    }
                    else if ((1 + G2_g) / (1 + G1_g) > T)
                    {
                        p_g = p2_g;
                        p_g = (-1 * imageCopyExtended_DCIM[i + 3, j - 3].g + 9 * imageCopyExtended_DCIM[i + 1, j - 1].g + 9 * imageCopyExtended_DCIM[i - 1, j + 1].g - 1 * imageCopyExtended_DCIM[i - 3, j + 3].g) / 16;
                    }
                    else
                    {
                        p1_g = (-1 * imageCopyExtended_DCIM[i - 3, j - 3].g + 9 * imageCopyExtended_DCIM[i - 1, j - 1].g + 9 * imageCopyExtended_DCIM[i + 1, j + 1].g - 1 * imageCopyExtended_DCIM[i + 3, j + 3].g) / 16;
                        p2_g = (-1 * imageCopyExtended_DCIM[i + 3, j - 3].g + 9 * imageCopyExtended_DCIM[i + 1, j - 1].g + 9 * imageCopyExtended_DCIM[i - 1, j + 1].g - 1 * imageCopyExtended_DCIM[i - 3, j + 3].g) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_g, k));
                        w2 = 1 / (1 + Math.Pow(G2_g, k));
                        p_g = (w1 * p1_g + w2 * p2_g) / (w1 + w2);
                    }



                    if ((1 + G1_b) / (1 + G2_b) > T)
                    {
                        p_b = p1_b;
                        p_b = (-1 * imageCopyExtended_DCIM[i - 3, j - 3].b + 9 * imageCopyExtended_DCIM[i - 1, j - 1].b + 9 * imageCopyExtended_DCIM[i + 1, j + 1].b - 1 * imageCopyExtended_DCIM[i + 3, j + 3].b) / 16;
                    }
                    else if ((1 + G2_b) / (1 + G1_b) > T)
                    {
                        p_b = p2_b;
                        p_b = (-1 * imageCopyExtended_DCIM[i + 3, j - 3].b + 9 * imageCopyExtended_DCIM[i + 1, j - 1].b + 9 * imageCopyExtended_DCIM[i - 1, j + 1].b - 1 * imageCopyExtended_DCIM[i - 3, j + 3].b) / 16;
                    }
                    else
                    {
                        p1_b = (-1 * imageCopyExtended_DCIM[i - 3, j - 3].b + 9 * imageCopyExtended_DCIM[i - 1, j - 1].b + 9 * imageCopyExtended_DCIM[i + 1, j + 1].b - 1 * imageCopyExtended_DCIM[i + 3, j + 3].b) / 16;
                        p2_b = (-1 * imageCopyExtended_DCIM[i + 3, j - 3].b + 9 * imageCopyExtended_DCIM[i + 1, j - 1].b + 9 * imageCopyExtended_DCIM[i - 1, j + 1].b - 1 * imageCopyExtended_DCIM[i - 3, j + 3].b) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_b, k));
                        w2 = 1 / (1 + Math.Pow(G2_b, k));
                        p_b = (w1 * p1_b + w2 * p2_b) / (w1 + w2);
                    }
                    ColorFloatPixel p = new ColorFloatPixel();
                    p.r = (float)p_r;
                    p.g = (float)p_g;
                    p.b = (float)p_b;
                    imageCopyExtended_DCIM[i, j] = p;

                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////step 2
                    ///////////////////////////////////////////////////////////////////////////////////


                }
            }

            for (int i2 = radius + 1; i2 < imageCopyExtended.Width - radius - 1; i2++)
            {
                for (int j2 = radius + 1; j2 < imageCopyExtended.Height - radius - 1; j2++)
                {
                    int i = i2 * radius + 1;
                    int j = j2 * radius;
                    double G1_r = 0;
                    double G2_r = 0;
                    double G1_g = 0;
                    double G2_g = 0;
                    double G1_b = 0;
                    double G2_b = 0;
                    double T = 1.15;
                    double p_r = 0, p1_r = 0, p2_r = 0;
                    double p_g = 0, p1_g = 0, p2_g = 0;
                    double p_b = 0, p1_b = 0, p2_b = 0;

                    G1_r = Math.Abs(imageCopyExtended_DCIM[i + 1, j].r - imageCopyExtended_DCIM[i + 1, j + 2].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2].r - imageCopyExtended_DCIM[i + 1, j].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j].r - imageCopyExtended_DCIM[i - 1, j + 2].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].r - imageCopyExtended_DCIM[i + 1, j].r)
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1].r - imageCopyExtended_DCIM[i, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1].r - imageCopyExtended_DCIM[i + 2, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].r - imageCopyExtended_DCIM[i - 2, j + 1].r);

                    G2_r = Math.Abs(imageCopyExtended_DCIM[i, j + 1].r - imageCopyExtended_DCIM[i + 2, j + 1].r)
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1].r - imageCopyExtended_DCIM[i + 2, j - 1].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1].r - imageCopyExtended_DCIM[i, j + 1].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].r - imageCopyExtended_DCIM[i, j - 1].r)
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j].r - imageCopyExtended_DCIM[i + 1, j].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2].r - imageCopyExtended_DCIM[i + 1, j + 2].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].r - imageCopyExtended_DCIM[i + 1, j - 2].r);


                    G1_g = Math.Abs(imageCopyExtended_DCIM[i + 1, j].g - imageCopyExtended_DCIM[i + 1, j + 2].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2].g - imageCopyExtended_DCIM[i + 1, j].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j].g - imageCopyExtended_DCIM[i - 1, j + 2].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].g - imageCopyExtended_DCIM[i + 1, j].g)
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1].g - imageCopyExtended_DCIM[i, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1].g - imageCopyExtended_DCIM[i + 2, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].g - imageCopyExtended_DCIM[i - 2, j + 1].g);

                    G2_g = Math.Abs(imageCopyExtended_DCIM[i, j + 1].g - imageCopyExtended_DCIM[i + 2, j + 1].g)
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1].g - imageCopyExtended_DCIM[i + 2, j - 1].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1].g - imageCopyExtended_DCIM[i, j + 1].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].g - imageCopyExtended_DCIM[i, j - 1].g)
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j].g - imageCopyExtended_DCIM[i + 1, j].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2].g - imageCopyExtended_DCIM[i + 1, j + 2].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].g - imageCopyExtended_DCIM[i + 1, j - 2].g);

                    G1_b = Math.Abs(imageCopyExtended_DCIM[i + 1, j].b - imageCopyExtended_DCIM[i + 1, j + 2].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2].b - imageCopyExtended_DCIM[i + 1, j].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j].b - imageCopyExtended_DCIM[i - 1, j + 2].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].b - imageCopyExtended_DCIM[i + 1, j].b)
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1].b - imageCopyExtended_DCIM[i, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1].b - imageCopyExtended_DCIM[i + 2, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].b - imageCopyExtended_DCIM[i - 2, j + 1].b);

                    G2_b = Math.Abs(imageCopyExtended_DCIM[i, j + 1].b - imageCopyExtended_DCIM[i + 2, j + 1].b)
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1].b - imageCopyExtended_DCIM[i + 2, j - 1].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1].b - imageCopyExtended_DCIM[i, j + 1].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].b - imageCopyExtended_DCIM[i, j - 1].b)
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j].b - imageCopyExtended_DCIM[i + 1, j].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2].b - imageCopyExtended_DCIM[i + 1, j + 2].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].b - imageCopyExtended_DCIM[i + 1, j - 2].b);


                    if ((1 + G1_r) / (1 + G2_r) > T)
                    {
                        p_r = p1_r;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p_r = (-1 * imageCopyExtended_DCIM[i, j - 3].r + 9 * imageCopyExtended_DCIM[i, j - 1].r + 9 * imageCopyExtended_DCIM[i, j + 1].r - 1 * imageCopyExtended_DCIM[i, j + 3].r) / 16;
                    }
                    else if ((1 + G2_r) / (1 + G1_r) > T)
                    {
                        p_r = p2_r;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p_r = (-1 * imageCopyExtended_DCIM[i - 3, j].r + 9 * imageCopyExtended_DCIM[i - 1, j].r + 9 * imageCopyExtended_DCIM[i + 1, j].r - 1 * imageCopyExtended_DCIM[i + 3, j].r) / 16;
                    }
                    else
                    {
                        p1_r = (-1 * imageCopyExtended_DCIM[i, j - 3].r + 9 * imageCopyExtended_DCIM[i, j - 1].r + 9 * imageCopyExtended_DCIM[i, j + 1].r - 1 * imageCopyExtended_DCIM[i, j + 3].r) / 16;
                        p2_r = (-1 * imageCopyExtended_DCIM[i - 3, j].r + 9 * imageCopyExtended_DCIM[i - 1, j].r + 9 * imageCopyExtended_DCIM[i + 1, j].r - 1 * imageCopyExtended_DCIM[i + 3, j].r) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_r, k));
                        w2 = 1 / (1 + Math.Pow(G2_r, k));
                        p_r = (w1 * p1_r + w2 * p2_r) / (w1 + w2);
                    }

                    if ((1 + G1_g) / (1 + G2_g) > T)
                    {
                        p_g = p1_g;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p_g = (-1 * imageCopyExtended_DCIM[i, j - 3].g + 9 * imageCopyExtended_DCIM[i, j - 1].g + 9 * imageCopyExtended_DCIM[i, j + 1].g - 1 * imageCopyExtended_DCIM[i, j + 3].g) / 16;
                    }
                    else if ((1 + G2_g) / (1 + G1_g) > T)
                    {
                        p_g = p2_g;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p_g = (-1 * imageCopyExtended_DCIM[i - 3, j].g + 9 * imageCopyExtended_DCIM[i - 1, j].g + 9 * imageCopyExtended_DCIM[i + 1, j].g - 1 * imageCopyExtended_DCIM[i + 3, j].g) / 16;
                    }
                    else
                    {
                        p1_g = (-1 * imageCopyExtended_DCIM[i, j - 3].g + 9 * imageCopyExtended_DCIM[i, j - 1].g + 9 * imageCopyExtended_DCIM[i, j + 1].g - 1 * imageCopyExtended_DCIM[i, j + 3].g) / 16;
                        p2_g = (-1 * imageCopyExtended_DCIM[i - 3, j].g + 9 * imageCopyExtended_DCIM[i - 1, j].g + 9 * imageCopyExtended_DCIM[i + 1, j].g - 1 * imageCopyExtended_DCIM[i + 3, j].g) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_g, k));
                        w2 = 1 / (1 + Math.Pow(G2_g, k));
                        p_g = (w1 * p1_g + w2 * p2_g) / (w1 + w2);
                    }

                    if ((1 + G1_b) / (1 + G2_b) > T)
                    {
                        p_b = p1_b;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p_b = (-1 * imageCopyExtended_DCIM[i, j - 3].b + 9 * imageCopyExtended_DCIM[i, j - 1].b + 9 * imageCopyExtended_DCIM[i, j + 1].b - 1 * imageCopyExtended_DCIM[i, j + 3].b) / 16;
                    }
                    else if ((1 + G2_b) / (1 + G1_b) > T)
                    {
                        p_b = p2_b;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p_b = (-1 * imageCopyExtended_DCIM[i - 3, j].b + 9 * imageCopyExtended_DCIM[i - 1, j].b + 9 * imageCopyExtended_DCIM[i + 1, j].b - 1 * imageCopyExtended_DCIM[i + 3, j].b) / 16;
                    }
                    else
                    {
                        p1_b = (-1 * imageCopyExtended_DCIM[i, j - 3].b + 9 * imageCopyExtended_DCIM[i, j - 1].b + 9 * imageCopyExtended_DCIM[i, j + 1].b - 1 * imageCopyExtended_DCIM[i, j + 3].b) / 16;
                        p2_b = (-1 * imageCopyExtended_DCIM[i - 3, j].b + 9 * imageCopyExtended_DCIM[i - 1, j].b + 9 * imageCopyExtended_DCIM[i + 1, j].b - 1 * imageCopyExtended_DCIM[i + 3, j].b) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_b, k));
                        w2 = 1 / (1 + Math.Pow(G2_b, k));
                        p_b = (w1 * p1_b + w2 * p2_b) / (w1 + w2);
                    }

                    ColorFloatPixel p = new ColorFloatPixel();
                    p.r = (float)p_r;
                    p.g = (float)p_g;
                    p.b = (float)p_b;
                    imageCopyExtended_DCIM[i, j] = p;
                    


                }
            }

            for (int i2 = radius + 1; i2 < imageCopyExtended.Width - radius - 1; i2++)
            {
                for (int j2 = radius + 1; j2 < imageCopyExtended.Height - radius - 1; j2++)
                {
                    int i = i2 * radius;
                    int j = j2 * radius + 1;

                    double T = 1.15;


                    double G1_r = 0;
                    double G2_r = 0;
                    double G1_g = 0;
                    double G2_g = 0;
                    double G1_b = 0;
                    double G2_b = 0;
                    double p_r = 0, p1_r = 0, p2_r = 0;
                    double p_g = 0, p1_g = 0, p2_g = 0;
                    double p_b = 0, p1_b = 0, p2_b = 0;

                    G1_r = Math.Abs(imageCopyExtended_DCIM[i + 1, j].r - imageCopyExtended_DCIM[i + 1, j + 2].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2].r - imageCopyExtended_DCIM[i + 1, j].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j].r - imageCopyExtended_DCIM[i - 1, j + 2].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].r - imageCopyExtended_DCIM[i + 1, j].r)
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1].r - imageCopyExtended_DCIM[i, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1].r - imageCopyExtended_DCIM[i + 2, j + 1].r)
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].r - imageCopyExtended_DCIM[i - 2, j + 1].r);

                    G2_r = Math.Abs(imageCopyExtended_DCIM[i, j + 1].r - imageCopyExtended_DCIM[i + 2, j + 1].r)
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1].r - imageCopyExtended_DCIM[i + 2, j - 1].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1].r - imageCopyExtended_DCIM[i, j + 1].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].r - imageCopyExtended_DCIM[i, j - 1].r)
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j].r - imageCopyExtended_DCIM[i + 1, j].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2].r - imageCopyExtended_DCIM[i + 1, j + 2].r)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].r - imageCopyExtended_DCIM[i + 1, j - 2].r);
                    G1_g = Math.Abs(imageCopyExtended_DCIM[i + 1, j].g - imageCopyExtended_DCIM[i + 1, j + 2].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2].g - imageCopyExtended_DCIM[i + 1, j].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j].g - imageCopyExtended_DCIM[i - 1, j + 2].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].g - imageCopyExtended_DCIM[i + 1, j].g)
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1].g - imageCopyExtended_DCIM[i, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1].g - imageCopyExtended_DCIM[i + 2, j + 1].g)
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].g - imageCopyExtended_DCIM[i - 2, j + 1].g);

                    G2_g = Math.Abs(imageCopyExtended_DCIM[i, j + 1].g - imageCopyExtended_DCIM[i + 2, j + 1].g)
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1].g - imageCopyExtended_DCIM[i + 2, j - 1].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1].g - imageCopyExtended_DCIM[i, j + 1].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].g - imageCopyExtended_DCIM[i, j - 1].g)
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j].g - imageCopyExtended_DCIM[i + 1, j].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2].g - imageCopyExtended_DCIM[i + 1, j + 2].g)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].g - imageCopyExtended_DCIM[i + 1, j - 2].g);

                    G1_b = Math.Abs(imageCopyExtended_DCIM[i + 1, j].b - imageCopyExtended_DCIM[i + 1, j + 2].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 1, j - 2].b - imageCopyExtended_DCIM[i + 1, j].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j].b - imageCopyExtended_DCIM[i - 1, j + 2].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].b - imageCopyExtended_DCIM[i + 1, j].b)
                        //horizontal
                        + Math.Abs(imageCopyExtended_DCIM[i, j - 1].b - imageCopyExtended_DCIM[i, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i + 2, j - 1].b - imageCopyExtended_DCIM[i + 2, j + 1].b)
                        + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].b - imageCopyExtended_DCIM[i - 2, j + 1].b);

                    G2_b = Math.Abs(imageCopyExtended_DCIM[i, j + 1].b - imageCopyExtended_DCIM[i + 2, j + 1].b)
                      + Math.Abs(imageCopyExtended_DCIM[i, j - 1].b - imageCopyExtended_DCIM[i + 2, j - 1].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j + 1].b - imageCopyExtended_DCIM[i, j + 1].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 2, j - 1].b - imageCopyExtended_DCIM[i, j - 1].b)
                        //vertical
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j].b - imageCopyExtended_DCIM[i + 1, j].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j + 2].b - imageCopyExtended_DCIM[i + 1, j + 2].b)
                      + Math.Abs(imageCopyExtended_DCIM[i - 1, j - 2].b - imageCopyExtended_DCIM[i + 1, j - 2].b);


                    if ((1 + G1_r) / (1 + G2_r) > T)
                    {
                        p_r = p1_r;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p_r = (-1 * imageCopyExtended_DCIM[i, j - 3].r + 9 * imageCopyExtended_DCIM[i, j - 1].r + 9 * imageCopyExtended_DCIM[i, j + 1].r - 1 * imageCopyExtended_DCIM[i, j + 3].r) / 16;
                    }
                    else if ((1 + G2_r) / (1 + G1_r) > T)
                    {
                        p_r = p2_r;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p_r = (-1 * imageCopyExtended_DCIM[i - 3, j].r + 9 * imageCopyExtended_DCIM[i - 1, j].r + 9 * imageCopyExtended_DCIM[i + 1, j].r - 1 * imageCopyExtended_DCIM[i + 3, j].r) / 16;
                    }
                    else
                    {
                        p1_r = (-1 * imageCopyExtended_DCIM[i, j - 3].r + 9 * imageCopyExtended_DCIM[i, j - 1].r + 9 * imageCopyExtended_DCIM[i, j + 1].r - 1 * imageCopyExtended_DCIM[i, j + 3].r) / 16;
                        p2_r = (-1 * imageCopyExtended_DCIM[i - 3, j].r + 9 * imageCopyExtended_DCIM[i - 1, j].r + 9 * imageCopyExtended_DCIM[i + 1, j].r - 1 * imageCopyExtended_DCIM[i + 3, j].r) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_r, k));
                        w2 = 1 / (1 + Math.Pow(G2_r, k));
                        p_r = (w1 * p1_r + w2 * p2_r) / (w1 + w2);
                    }

                    if ((1 + G1_g) / (1 + G2_g) > T)
                    {
                        p_g = p1_g;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p_g = (-1 * imageCopyExtended_DCIM[i, j - 3].g + 9 * imageCopyExtended_DCIM[i, j - 1].g + 9 * imageCopyExtended_DCIM[i, j + 1].g - 1 * imageCopyExtended_DCIM[i, j + 3].g) / 16;
                    }
                    else if ((1 + G2_g) / (1 + G1_g) > T)
                    {
                        p_g = p2_g;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p_g = (-1 * imageCopyExtended_DCIM[i - 3, j].g + 9 * imageCopyExtended_DCIM[i - 1, j].g + 9 * imageCopyExtended_DCIM[i + 1, j].g - 1 * imageCopyExtended_DCIM[i + 3, j].g) / 16;
                    }
                    else
                    {
                        p1_g = (-1 * imageCopyExtended_DCIM[i, j - 3].g + 9 * imageCopyExtended_DCIM[i, j - 1].g + 9 * imageCopyExtended_DCIM[i, j + 1].g - 1 * imageCopyExtended_DCIM[i, j + 3].g) / 16;
                        p2_g = (-1 * imageCopyExtended_DCIM[i - 3, j].g + 9 * imageCopyExtended_DCIM[i - 1, j].g + 9 * imageCopyExtended_DCIM[i + 1, j].g - 1 * imageCopyExtended_DCIM[i + 3, j].g) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_g, k));
                        w2 = 1 / (1 + Math.Pow(G2_g, k));
                        p_g = (w1 * p1_g + w2 * p2_g) / (w1 + w2);
                    }

                    if ((1 + G1_b) / (1 + G2_b) > T)
                    {
                        p_b = p1_b;
                        // (-1 * P(X, Y - 3) +                     9 * P(X, Y - 1) +                   9 * P(X, Y + 1) -                      1 * P(X, Y + 3)) / 16
                        p_b = (-1 * imageCopyExtended_DCIM[i, j - 3].b + 9 * imageCopyExtended_DCIM[i, j - 1].b + 9 * imageCopyExtended_DCIM[i, j + 1].b - 1 * imageCopyExtended_DCIM[i, j + 3].b) / 16;
                    }
                    else if ((1 + G2_b) / (1 + G1_b) > T)
                    {
                        p_b = p2_b;
                        //(-1 *            P(X - 3, Y) +              9 * P(X - 1, Y) +                9 * P(X + 1, Y) -                  1 * P(X + 3, Y)) / 16
                        p_b = (-1 * imageCopyExtended_DCIM[i - 3, j].b + 9 * imageCopyExtended_DCIM[i - 1, j].b + 9 * imageCopyExtended_DCIM[i + 1, j].b - 1 * imageCopyExtended_DCIM[i + 3, j].b) / 16;
                    }
                    else
                    {
                        p1_b = (-1 * imageCopyExtended_DCIM[i, j - 3].b + 9 * imageCopyExtended_DCIM[i, j - 1].b + 9 * imageCopyExtended_DCIM[i, j + 1].b - 1 * imageCopyExtended_DCIM[i, j + 3].b) / 16;
                        p2_b = (-1 * imageCopyExtended_DCIM[i - 3, j].b + 9 * imageCopyExtended_DCIM[i - 1, j].b + 9 * imageCopyExtended_DCIM[i + 1, j].b - 1 * imageCopyExtended_DCIM[i + 3, j].b) / 16;
                        double w1, w2;
                        int k = 5;
                        w1 = 1 / (1 + Math.Pow(G1_b, k));
                        w2 = 1 / (1 + Math.Pow(G2_b, k));
                        p_b = (w1 * p1_b + w2 * p2_b) / (w1 + w2);
                    }

                    ColorFloatPixel p = new ColorFloatPixel();
                    p.r = (float)p_r;
                    p.g = (float)p_g;
                    p.b = (float)p_b;
                    imageCopyExtended_DCIM[i, j] = p;

                    //imageCopyExtended_DCIM[i, j] = 255;


                }
            }


            for (int i = 0; i < imageResult.Width; i++)
            {
                for (int j = 0; j < imageResult.Height; j++)
                {
                    imageResult[i, j] = imageCopyExtended_DCIM[i + 2 * radius_more, j + 2 * radius_more];

                }
            }
            ImageIO.ImageToFile(imageResult, fileOutName);
        }

        static double rgb2grey(float r, float g, float b)
        {
            return 0.114f * b + 0.587f * g + 0.299f * r;
        }



        static void BillinearInterpolation(ColorFloatImage image, float radius, string fileOutName)
        {


            ColorFloatImage imageCopyExtended = new ColorFloatImage((int)(image.Width * radius), (int)(image.Height * radius));


            for (int x = 0; x < imageCopyExtended.Height ; x++)
            {
                for (int y = 0; y < imageCopyExtended.Width; y++)
                {

                    float x0 = (float)y * image.Width / imageCopyExtended.Width; //Red.GetLength(1) / Red1.GetLength(1);
                    float y0 = (float)x * image.Height / imageCopyExtended.Height; //Red.GetLength(0) / Red1.GetLength(0);


                    int y1 = 0, y2 = 0;
                    int x1 = 0, x2 = 0;
                    x1 = (int)Math.Floor((double)x0);
                    x2 = (int)Math.Min(x1 + 1, image.Width - 1);
                    y1 = (int)Math.Floor((double)y0);
                    y2 = (int)Math.Min(y1 + 1, image.Height - 1);
                    int[] rows = new int[2] { y1, y2 };
                    int[] cols = new int[2] { x1, x2 };

                    if (x2 >= image.Width)
                    {
                        x2 = image.Width - 1;
                    }
                    if (y2 >= image.Height)
                    {
                        y2 = image.Height - 1;
                    }
                    if (x1 < 0)
                    {
                        x1 = 0;
                    }
                    if (y1 < 0)
                    {
                        y1 = 0;
                    }

                    float[,] matrix = new float[2, 2] 
                        {
                            {1, -1},
                            {-1, 1}
                        };
                    ColorFloatPixel pixel = new ColorFloatPixel();

                    for (int i = 0; i < 2; ++i)
                    {
                        for (int j = 0; j < 2; ++j)
                        {
                            float kernel = matrix[i, j];
                            for (int k = 0; k < 2; ++k)
                            {
                                if (k != i)
                                    kernel *= (y0 - (int)(y0) - k );
                                if (k != j)
                                    kernel *= (x0 - (int)(x0) - k );
                            }
                            
                            pixel.r += kernel * image[cols[j], rows[i]].r;
                            pixel.g += kernel * image[cols[j], rows[i]].g;
                            pixel.b += kernel * image[cols[j], rows[i]].b;

                            

                        }


                    }
                    imageCopyExtended[y, x] = pixel;

                }
            }

                ImageIO.ImageToFile( imageCopyExtended , fileOutName);


        }


        static void BicubicInterpolation(ColorFloatImage image, float radius, string fileOutName)
        {


            ColorFloatImage imageCopyExtended = new ColorFloatImage((int)(image.Width * radius), (int)(image.Height * radius));

            float[,] Coeff = new float[4, 4]         
            {{(float)1 / 36, (float)-1 / 12, (float)1 / 12, (float)-1 / 36},
            {(float)-1 / 12, (float)1 / 4, (float)-1 / 4, (float)1 / 12},
            {(float)1 / 12, (float)-1 / 4, (float)1 / 4, (float)-1 / 12},
            {(float)-1 / 36, (float)1 / 12, (float)-1 / 12, (float)1 / 36}};


            for (int x0 = 0; x0 < imageCopyExtended.Height ; ++x0)
            {
                for (int y0 = 0; y0 < imageCopyExtended.Width ; ++y0)
                {
                    float y = (float)y0 * image.Width / imageCopyExtended.Width; //Red.GetLength(1) / Red1.GetLength(1);
                    float x = (float)x0 * image.Height / imageCopyExtended.Height;//Red.GetLength(0) / Red1.GetLength(0);
                    int DownX = (int)Math.Floor(x);
                    int x1 = (int)Math.Max(DownX - 1, 0),
                        x2 = (int)Math.Max(DownX, 0),
                        x3 = (int)Math.Min(DownX + 1, image.Height - 1),
                        x4 = (int)Math.Min(DownX + 2, image.Height - 1);
                    int[] rows = new int[4] { x1, x2, x3, x4 };
                    int DownY = (int)Math.Floor(y);
                    int y1 = (int)Math.Max(DownY - 1, 0),
                        y2 = (int)Math.Max(DownY, 0),
                        y3 = (int)Math.Min(DownY + 1, image.Width - 1),
                        y4 = (int)Math.Min(DownY + 2, image.Width - 1);
                    int[] cols = new int[4] { y1, y2, y3, y4 };

                    ColorFloatPixel pixel = new ColorFloatPixel();
                    for (int i = 0; i < 4; ++i)
                    {
                        for (int j = 0; j < 4; ++j)
                        {
                            float CoeffFunc = Coeff[i, j];
                            for (int k = 0; k < 4; ++k)
                            {
                                if (k != i)
                                    CoeffFunc = CoeffFunc * (x - (int)(x) - k + 1);
                                if (k != j)
                                    CoeffFunc = CoeffFunc * (y - (int)(y) - k + 1);
                            }


                            pixel.r += CoeffFunc * image[cols[j], rows[i]].r;
                            pixel.g += CoeffFunc * image[cols[j], rows[i]].g;
                            pixel.b += CoeffFunc * image[cols[j], rows[i]].b;

                        }
                    }

                    imageCopyExtended[y0, x0] = pixel;
                }
            }




            ImageIO.ImageToFile(imageCopyExtended, fileOutName);


        }

        static void DownSample(ColorFloatImage image, float radius, string fileOutName)
        {


            ColorFloatImage imageCopyExtended = new ColorFloatImage((int)(image.Width / radius), (int)(image.Height / radius));


            for (int x = 0; x < imageCopyExtended.Height; x++)
            {
                for (int y = 0; y < imageCopyExtended.Width; y++)
                {

                    float x0 = (float)y * radius;//image.Width * imageCopyExtended.Width; //Red.GetLength(1) / Red1.GetLength(1);
                    float y0 = (float)x * radius;//image.Height * imageCopyExtended.Height; //Red.GetLength(0) / Red1.GetLength(0);


                    int y1 = 0, y2 = 0;
                    int x1 = 0, x2 = 0;
                    x1 = (int)Math.Floor((double)x0);
                    x2 = (int)Math.Min(x1 + 1, image.Width - 1);
                    y1 = (int)Math.Floor((double)y0);
                    y2 = (int)Math.Min(y1 + 1, image.Height - 1);
                    int[] rows = new int[2] { y1, y2 };
                    int[] cols = new int[2] { x1, x2 };

                    if (x2 >= image.Width)
                    {
                        x2 = image.Width - 1;
                    }
                    if (y2 >= image.Height)
                    {
                        y2 = image.Height - 1;
                    }
                    if (x1 < 0)
                    {
                        x1 = 0;
                    }
                    if (y1 < 0)
                    {
                        y1 = 0;
                    }

                    float[,] matrix = new float[2, 2] 
                        {
                            {1, -1},
                            {-1, 1}
                        };
                    ColorFloatPixel pixel = new ColorFloatPixel();

                    for (int i = 0; i < 2; ++i)
                    {
                        for (int j = 0; j < 2; ++j)
                        {
                            float kernel = matrix[i, j];
                            for (int k = 0; k < 2; ++k)
                            {
                                if (k != i)
                                    kernel *= (y0 - (int)(y0) - k);
                                if (k != j)
                                    kernel *= (x0 - (int)(x0) - k);
                            }
                            pixel.r += kernel * image[cols[j], rows[i]].r;
                            pixel.g += kernel * image[cols[j], rows[i]].g;
                            pixel.b += kernel * image[cols[j], rows[i]].b;
                        }


                    }
                    imageCopyExtended[y, x] = pixel;

                }
            }

            ImageIO.ImageToFile(imageCopyExtended, fileOutName);


        }
        static double mse_metric(GrayscaleFloatImage image1,GrayscaleFloatImage image2)
        {
            double mectric_mse = 0.0;//1/(image1.Height*image1.Width*3);

            for (int y = 0; y < image1.Height;y++ )
            {
                for (int x=0;x<image1.Width;x++)
                {
                    mectric_mse += (image1[x, y] - image2[x, y]) * (image1[x, y] - image2[x, y]);

                }
            }
            mectric_mse = mectric_mse /(image1.Height*image1.Width);


                return mectric_mse;
        }
        static double psnr_metric(GrayscaleFloatImage image1, GrayscaleFloatImage image2)
        {
            double mectric_psnr = 0.0;//1/(image1.Height*image1.Width*3);

            mectric_psnr = 20 * Math.Log10(255) - 10 * Math.Log10(mse_metric(image1, image2));

            return mectric_psnr;
        }
        static double ssim_metric(GrayscaleFloatImage image1, GrayscaleFloatImage image2 )
        {
            double metric_ssim = 0.0;
            double x_avg = 0;
            double y_avg = 0;

            for (int i = 0; i < image1.Width;i++ )
            {
                for (int j=0; j<image1.Height; j++)
                {
                    x_avg += image1[i, j];
                    y_avg += image2[i, j];
                }
            }
            x_avg = x_avg / (image1.Height * image1.Width);
            y_avg = y_avg / (image2.Height * image2.Width);

            double sigma_x_2 = 0;
            double sigma_y_2 = 0;
            double sigma_xy_2 = 0;

            for (int i = 0; i < image1.Width; i++)
            {
                for (int j = 0; j < image1.Height; j++)
                {
                    sigma_x_2 += (image1[i, j] - x_avg) * (image1[i, j] - x_avg);
                    sigma_y_2 += (image2[i, j] - y_avg) * (image2[i, j] - y_avg);
                    sigma_xy_2 += (image1[i, j] - x_avg) * (image2[i, j] - y_avg);
                }
            }

            sigma_x_2 = sigma_x_2 / (image1.Height * image1.Width-1);
            sigma_y_2 = sigma_y_2 / (image2.Height * image2.Width - 1);
            sigma_xy_2 = sigma_xy_2 / (image1.Height * image1.Width - 1);

            double c1=(255*0.01)*(255*0.01);
            double c2 = (255*0.03)*(255*0.03);

            metric_ssim = ((2 * x_avg * y_avg + c1) * (2 * sigma_xy_2 + c2)) / 
                ((x_avg * x_avg + y_avg * y_avg + c1) * (sigma_x_2 + sigma_y_2 + c2));

            return metric_ssim;
        }
        static double mssim_metric(GrayscaleFloatImage image1, GrayscaleFloatImage image2)
        {
            double metric_ssim = 0.0;
            double metric_mssim = 0.0;
            int msim_count = 0;

            int N= 8;


            for (int x = 0; x < image1.Width; x+=8)
            {
                for (int y = 0; y < image1.Height; y+=8)
                {
                    msim_count++;
                    double x_avg = 0;
                    double y_avg = 0;
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            x_avg += image1[i+x, j+y];
                            y_avg += image2[i+x, j+y];
                        }
                    }
                    x_avg = x_avg / (N * N);
                    y_avg = y_avg / (N * N);

                    double sigma_x_2 = 0;
                    double sigma_y_2 = 0;
                    double sigma_xy_2 = 0;

                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            sigma_x_2 += (image1[i+x, j+y] - x_avg) * (image1[i+x, j+y] - x_avg);
                            sigma_y_2 += (image2[i+x, j+y] - y_avg) * (image2[i+x, j+y] - y_avg);
                            sigma_xy_2 += (image1[i+x, j+y] - x_avg) * (image2[i+x, j+y] - y_avg);
                        }
                    }

                    sigma_x_2 = sigma_x_2 / (N * N - 1);
                    sigma_y_2 = sigma_y_2 / (N * N - 1);
                    sigma_xy_2 = sigma_xy_2 / (N * N - 1);

                    double c1 = (255 * 0.01) * (255 * 0.01);
                    double c2 = (255 * 0.03) * (255 * 0.03);

                    metric_ssim = ((2 * x_avg * y_avg + c1) * (2 * sigma_xy_2 + c2)) /
                        ((x_avg * x_avg + y_avg * y_avg + c1) * (sigma_x_2 + sigma_y_2 + c2));

                    metric_mssim += metric_ssim;

                }
            }
            metric_mssim = metric_mssim / msim_count;
            return metric_mssim;
        }

        public static GrayscaleFloatImage ApplyOneDimentionalFilter(float sigma, GrayscaleFloatImage img, int radius, string direction, bool derivative)
        {
            GrayscaleFloatImage outImg = new GrayscaleFloatImage(img.Width, img.Height);
            float[,] imageTmpGray = new float[img.Width + 2 * radius, img.Height + 2 * radius];
            float[] Kernel = new float[radius * 2 + 1];
            float SumKernel = 0;


            //копируем в середину
            for (int x = 0; x < img.Width; x++)
            {
                for (int y = 0; y < img.Height; y++)
                {
                    imageTmpGray[x + radius, y + radius] = img[x, y];
                }
            }
            //копируем края
            for (int x = radius; x < img.Width + radius; x++)
            {
                for (int y = 0; y < radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x, y + radius];

                }
            }

            for (int x = radius; x < img.Width + radius; x++)
            {
                for (int y = img.Height + radius; y < img.Height + 2 * radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x, y - radius];

                }
            }
            for (int x = 0; x < radius; x++)
            {
                for (int y = radius; y < img.Height + radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x + radius, y];

                }
            }
            for (int x = img.Width + radius; x < img.Width + 2 * radius; x++)
            {
                for (int y = radius; y < img.Height + radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x - radius, y];

                }
            }
            /////////////////////////////////////////////////////
            //копируем уголки
            for (int x = 0; x < radius; x++)
            {
                for (int y = 0; y < radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x + radius, y + radius];

                }
            }

            for (int x = img.Width + radius; x < img.Width + 2 * radius; x++)
            {
                for (int y = img.Height + radius; y < img.Height + 2 * radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x - radius - 1, y - radius - 1];

                }
            }
            for (int x = img.Width + radius; x < img.Width + 2 * radius; x++)
            {
                for (int y = 0; y < radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x - radius, y + radius];

                }
            }
            for (int x = 0; x < radius; x++)
            {
                for (int y = img.Height + radius; y < img.Height + 2 * radius; y++)
                {

                    imageTmpGray[x, y] = imageTmpGray[x + radius, y - radius];

                }
            }

            ////////////////////////////////////////////////////
            //проверка

            GrayscaleFloatImage imageTmpGray_img = new GrayscaleFloatImage(img.Width + 2 * radius, img.Height + 2 * radius);
            for (int x = 0; x < img.Width + 2 * radius; x++)
            {
                for (int y = 0; y < img.Height + 2 * radius; y++)
                {
                    imageTmpGray_img[x, y] = imageTmpGray[x, y];
                }
            }
            ImageIO.ImageToFile(imageTmpGray_img, "Canny_tmp.png");
            ////////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////

            //Console.WriteLine("GausKernel_without norm");
            for (int k = -radius; k < radius + 1; k++)
            {
                if (derivative)
                {
                    Kernel[k + radius] = (float)GaussFunctionDerivative(sigma, k);
                }
                else
                {
                    Kernel[k + radius] = (float)GaussFunction(sigma, k);
                }

                SumKernel += Kernel[k + radius];
                //Console.WriteLine(Kernel[k + radius] + "|");
            }

            //Console.WriteLine("GausKernel+norm");
            if (!derivative)
            {
                for (int k = -radius; k < radius + 1; k++)
                {
                    Kernel[k + radius] /= SumKernel;

                    //Console.WriteLine(Kernel[k + radius] + "|");
                }
            }

            float sum = 0;
            if (direction == "x")
            {
                for (int j = radius; j < imageTmpGray.GetLength(1) - radius; j++)
                {
                    for (int i = radius; i < imageTmpGray.GetLength(0) - radius; i++)
                    {
                        sum = 0;
                        for (int k = -radius; k < radius + 1; k++)
                        {
                            sum += imageTmpGray[i, j + k] * Kernel[k + radius];
                        }
                        if (derivative) { outImg[i - radius, j - radius] = sum; }//*100
                        else { outImg[i - radius, j - radius] = sum; }

                    }
                }
            }
            else //y
            {
                for (int j = radius; j < imageTmpGray.GetLength(1) - radius; j++)
                {
                    for (int i = radius; i < imageTmpGray.GetLength(0) - radius; i++)
                    {
                        sum = 0;
                        for (int k = -radius; k < radius + 1; k++)
                        {
                            sum += imageTmpGray[i + k, j] * Kernel[k + radius];
                        }
                        if (derivative) { outImg[i - radius, j - radius] = sum; }//*100
                        else { outImg[i - radius, j - radius] = sum; }
                    }
                }
            }

            if (derivative)
            {
                float max = outImg[0, 0];
                float min = outImg[0, 0];
                for (int j = 0; j < outImg.Height; j++)
                {
                    for (int i = 0; i < outImg.Width; i++)
                    {
                        if (outImg[i, j] > max) { max = outImg[i, j]; }
                        if (outImg[i, j] < min) { min = outImg[i, j]; }
                    }
                }

                for (int j = 0; j < outImg.Height; j++)
                {
                    for (int i = 0; i < outImg.Width; i++)
                    {
                        outImg[i, j] = (outImg[i, j] / (max)) * 256;
                    }
                }
            }


            return outImg;
        }


        public static GrayscaleFloatImage DerivativeX(float sigma, GrayscaleFloatImage img, int radius)
        {
            return ApplyOneDimentionalFilter(sigma, ApplyOneDimentionalFilter(sigma, img, radius, "x", false), radius, "y", true);
        }

        public static GrayscaleFloatImage DerivativeY(float sigma, GrayscaleFloatImage img, int radius)
        {
            return ApplyOneDimentionalFilter(sigma, ApplyOneDimentionalFilter(sigma, img, radius, "y", false), radius, "x", true);
        }

        static void Canny(GrayscaleFloatImage img, float sigma, float t1, float t2 , string OutFileName)
        {
            GrayscaleFloatImage dx, dy, grad, nms, res;
            dx = DerivativeX(sigma, img, (int)(2 * sigma));// или собель по х
            dy = DerivativeY(sigma, img, (int)(2 * sigma));//собель по y

            //ImageIO.ImageToFile(dy, "cannyTest/dy.png");
            //ImageIO.ImageToFile(dx, "cannyTest/dx.png");

            grad = new GrayscaleFloatImage(img.Width, img.Height);
            res = new GrayscaleFloatImage(img.Width, img.Height);
            nms = new GrayscaleFloatImage(img.Width, img.Height);
            for (int x = 0; x < img.Width; x++)
            {
                for (int y = 0; y < img.Height; y++)
                {
                    grad[x, y] = (float)Math.Sqrt(dx[x, y] * dx[x, y] + dy[x, y] * dy[x, y]);
                }
            }
            //ImageIO.ImageToFile(grad, "cannyTest/Gradient.png");
            for (int x = 1; x < img.Width - 1; x++)
            {
                for (int y = 1; y < img.Height - 1; y++)//странно, почему то не работает без -1
                {
                    float valueNMS;
                    if (x < 1 || x >= img.Width || y < 1 || y >= img.Height || y == img.Height || x == img.Width)
                        valueNMS = 0.0f;

                    float gval = grad[x, y];
                    float tan = (dx[x, y] == 0.0f ? 10000.0f : dy[x, y] / dx[x, y]);

                    if (tan > 2.5f || tan < -2.5f)  //vert
                    {
                        if (gval > grad[x, y - 1] && gval > grad[x, y + 1])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    else if (tan > 0.4f && tan <= 2.5f)//21.8-68.2
                    {
                        if (gval > grad[x - 1, y - 1] && gval > grad[x + 1, y + 1])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    else if (tan > -0.4f && tan <= 0.4f)//21.8-(-21.8)  (vert)
                    {
                        if (gval > grad[x - 1, y] && gval > grad[x + 1, y])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    else
                    {
                        if (gval > grad[x - 1, y + 1] && gval > grad[x + 1, y - 1])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    //res[x, y] =
                        nms[x, y] = valueNMS;
                    //res[x, y] =
                    //    nms[x, y] = NonMaxSupPixel(x, y, img);
                }
            }
            //ImageIO.ImageToFile(nms, "cannyTest/NMS.png");
            //grad = (dx * dx + dy * dy).Sqrt().Cast<float>();// вычисление градинтов, что есть dx?dy?и квадрат из них попиксельно ?
            //res = nms = NonMaximumSuppression(dx, dy, grad);// подавление не максимумов
            /////////////////////////////////////////////////////////////////////////////////
            //aplyAutoTreshhold
            //res = new GrayscaleFloatImage(grad.Width, grad.Height);
            //float thr_tweak = t1;//1.0f;
            //float max_g2 = grad.Max();// ммаксимальное значение градиента
            float max_g2 = nms[0, 0];
            for (int x = 0; x < nms.Width; x++)
            {
                for (int y = 0; y < nms.Height; y++)
                {
                    if (nms[x, y] > max_g2) { max_g2 = nms[x, y]; }
                }
            }


            float thr1 = max_g2 * t1;// *0.2f;//0.2f * thr_tweak;     // The most strongest line segments
            //float thr2 = thr1 * 0.1f;                   // Other line segments, maximum canny edge criterion
            float thr3 = max_g2 * t2;// 0.1f;                   // Minimum canny edge criterion
            ///////////////////////////////////////////////////////////////////////////////////
            //hysteresis
            //float thr_strongest=thr1;
            float thr_strong=thr3; 
            float thr_weak=thr1;
            for (int x = 0; x < res.Width; x++)
            {
                for (int y = 0; y < res.Height; y++)
                {
                    if (nms[x, y] > thr_strong) { 
                        res[x, y] = 255.0f; }
                    else if (nms[x, y] > thr_weak) { 
                        res[x, y] = 128.0f; }
                    else { 
                        res[x, y] = 0.0f; }//{ res[x, y] = 1.0f; }
                    //else { res[x, y] = 0.0f; }
                }
            }

                Queue<Point> qp = new Queue<Point>();  //создаем список(очередь) точек, градиент в которых больше большого порога
                for (int j = 0; j < nms.Height; j++)
                    for (int i = 0; i < nms.Width; i++)
                        if (nms[i, j] > thr_strong)
                            qp.Enqueue(new Point(i, j));

                //ImageIO.ImageToFile(res, "cannyTest/TestTreshold_beforeH.png"); 
                // Гистерезис
                while (qp.Count > 0)
                {
                    Point point = qp.Dequeue(); //берем первую точку 
                    if (point.X > 0 && res[point.X - 1, point.Y] == 128.0f && nms[point.X - 1, point.Y] > thr_weak)
                    {
                        res[point.X - 1, point.Y] = 255.0f;
                        qp.Enqueue(new Point(point.X - 1, point.Y));
                    }
                    if (point.X < nms.Width - 1 && res[point.X + 1, point.Y] == 128.0f && nms[point.X + 1, point.Y] > thr_weak)
                    {
                        res[point.X + 1, point.Y] = 255.0f;
                        qp.Enqueue(new Point(point.X + 1, point.Y));
                    }
                    if (point.Y > 0 && res[point.X, point.Y - 1] == 128.0f && nms[point.X, point.Y - 1] > thr_weak)
                    {
                        res[point.X, point.Y - 1] = 255.0f;
                        qp.Enqueue(new Point(point.X, point.Y - 1));
                    }
                    if (point.Y < nms.Height - 1 && res[point.X, point.Y + 1] == 128.0f && nms[point.X, point.Y + 1] > thr_weak)
                    {
                        res[point.X, point.Y + 1] = 255.0f;
                        qp.Enqueue(new Point(point.X, point.Y + 1));
                    }


                    if (point.Y < nms.Height - 1 && point.X < nms.Width - 1 && res[point.X + 1, point.Y + 1] == 128.0f && nms[point.X + 1, point.Y + 1] > thr_weak)
                    {
                        res[point.X+1, point.Y + 1] = 255.0f;
                        qp.Enqueue(new Point(point.X+1, point.Y + 1));
                    }
                    if (point.Y < nms.Height - 1 && point.X > 0 && res[point.X-1, point.Y + 1] == 128.0f && nms[point.X-1, point.Y + 1] > thr_weak)
                    {
                        res[point.X-1, point.Y + 1] = 255.0f;
                        qp.Enqueue(new Point(point.X-1, point.Y + 1));
                    }
                    if (point.X < nms.Width - 1 && point.Y > 0 && res[point.X + 1, point.Y - 1] == 128.0f && nms[point.X+1, point.Y - 1] > thr_weak)
                    {
                        res[point.X+1, point.Y - 1] = 255.0f;
                        qp.Enqueue(new Point(point.X+1, point.Y - 1));
                    }
                    if (point.Y > 0 && point.X > 0 && res[point.X-1, point.Y - 1] == 128.0f && nms[point.X-1, point.Y - 1] > thr_weak)
                    {
                        res[point.X -1, point.Y - 1] = 255.0f;
                        qp.Enqueue(new Point(point.X -1, point.Y - 1));
                    }


                }
                //ImageIO.ImageToFile(res, "cannyTest/TestTresholdmmm.png"); 

                for (int x = 0; x < res.Width; x++)
                {
                    for (int y = 0; y < res.Height; y++)
                    {
                        if (res[x, y] !=255)
                        {
                            res[x, y] = 0.0f;
                        }
                    }
                }

            ImageIO.ImageToFile(res, OutFileName);
        }
        static GrayscaleFloatImage cannyGetImg(GrayscaleFloatImage img, float sigma, float t1, float t2)
        {
            GrayscaleFloatImage dx, dy, grad, nms, res;
            dx = DerivativeX(sigma, img, (int)(2 * sigma));// или собель по х
            dy = DerivativeY(sigma, img, (int)(2 * sigma));//собель по y

            //ImageIO.ImageToFile(dy, "cannyTest/dy.png");
            //ImageIO.ImageToFile(dx, "cannyTest/dx.png");

            grad = new GrayscaleFloatImage(img.Width, img.Height);
            res = new GrayscaleFloatImage(img.Width, img.Height);
            nms = new GrayscaleFloatImage(img.Width, img.Height);
            for (int x = 0; x < img.Width; x++)
            {
                for (int y = 0; y < img.Height; y++)
                {
                    grad[x, y] = (float)Math.Sqrt(dx[x, y] * dx[x, y] + dy[x, y] * dy[x, y]);
                }
            }
            //ImageIO.ImageToFile(grad, "cannyTest/Gradient.png");
            for (int x = 1; x < img.Width - 1; x++)
            {
                for (int y = 1; y < img.Height - 1; y++)//странно, почему то не работает без -1
                {
                    float valueNMS;
                    if (x < 1 || x >= img.Width || y < 1 || y >= img.Height || y == img.Height || x == img.Width)
                        valueNMS = 0.0f;

                    float gval = grad[x, y];
                    float tan = (dx[x, y] == 0.0f ? 10000.0f : dy[x, y] / dx[x, y]);

                    if (tan > 2.5f || tan < -2.5f)  //vert
                    {
                        if (gval > grad[x, y - 1] && gval > grad[x, y + 1])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    else if (tan > 0.4f && tan <= 2.5f)//21.8-68.2
                    {
                        if (gval > grad[x - 1, y - 1] && gval > grad[x + 1, y + 1])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    else if (tan > -0.4f && tan <= 0.4f)//21.8-(-21.8)  (vert)
                    {
                        if (gval > grad[x - 1, y] && gval > grad[x + 1, y])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    else
                    {
                        if (gval > grad[x - 1, y + 1] && gval > grad[x + 1, y - 1])
                            valueNMS = gval;
                        else
                            valueNMS = 0.0f;
                    }
                    //res[x, y] =
                        nms[x, y] = valueNMS;
                    //res[x, y] =
                    //    nms[x, y] = NonMaxSupPixel(x, y, img);
                }
            }
           // ImageIO.ImageToFile(nms, "cannyTest/NMS.png");
            //grad = (dx * dx + dy * dy).Sqrt().Cast<float>();// вычисление градинтов, что есть dx?dy?и квадрат из них попиксельно ?
            //res = nms = NonMaximumSuppression(dx, dy, grad);// подавление не максимумов
            /////////////////////////////////////////////////////////////////////////////////
            //aplyAutoTreshhold
            //res = new GrayscaleFloatImage(grad.Width, grad.Height);
            //float thr_tweak = t1;//1.0f;
            //float max_g2 = grad.Max();// ммаксимальное значение градиента
            float max_g2 = nms[0, 0];
            for (int x = 0; x < nms.Width; x++)
            {
                for (int y = 0; y < nms.Height; y++)
                {
                    if (nms[x, y] > max_g2) { max_g2 = nms[x, y]; }
                }
            }


            float thr1 = max_g2 * t1;// *0.2f;//0.2f * thr_tweak;     // The most strongest line segments
            //float thr2 = thr1 * 0.1f;                   // Other line segments, maximum canny edge criterion
            float thr3 = max_g2 * t2;// 0.1f;                   // Minimum canny edge criterion
            ///////////////////////////////////////////////////////////////////////////////////
            //hysteresis
            //float thr_strongest=thr1;
            float thr_strong=thr3; 
            float thr_weak=thr1;
            for (int x = 0; x < res.Width; x++)
            {
                for (int y = 0; y < res.Height; y++)
                {
                    if (nms[x, y] > thr_strong) { 
                        res[x, y] = 255.0f; }
                    else if (nms[x, y] > thr_weak) { 
                        res[x, y] = 128.0f; }
                    else { 
                        res[x, y] = 0.0f; }//{ res[x, y] = 1.0f; }
                    //else { res[x, y] = 0.0f; }
                }
            }

                Queue<Point> qp = new Queue<Point>();  //создаем список(очередь) точек, градиент в которых больше большого порога
                for (int j = 0; j < nms.Height; j++)
                    for (int i = 0; i < nms.Width; i++)
                        if (nms[i, j] > thr_strong)
                            qp.Enqueue(new Point(i, j));

                //ImageIO.ImageToFile(res, "cannyTest/TestTreshold_beforeH.png"); 
                // Гистерезис
                while (qp.Count > 0)
                {
                    Point point = qp.Dequeue(); //берем первую точку 
                    if (point.X > 0 && res[point.X - 1, point.Y] == 128.0f && nms[point.X - 1, point.Y] > thr_weak)
                    {
                        res[point.X - 1, point.Y] = 255.0f;
                        qp.Enqueue(new Point(point.X - 1, point.Y));
                    }
                    if (point.X < nms.Width - 1 && res[point.X + 1, point.Y] == 128.0f && nms[point.X + 1, point.Y] > thr_weak)
                    {
                        res[point.X + 1, point.Y] = 255.0f;
                        qp.Enqueue(new Point(point.X + 1, point.Y));
                    }
                    if (point.Y > 0 && res[point.X, point.Y - 1] == 128.0f && nms[point.X, point.Y - 1] > thr_weak)
                    {
                        res[point.X, point.Y - 1] = 255.0f;
                        qp.Enqueue(new Point(point.X, point.Y - 1));
                    }
                    if (point.Y < nms.Height - 1 && res[point.X, point.Y + 1] == 128.0f && nms[point.X, point.Y + 1] > thr_weak)
                    {
                        res[point.X, point.Y + 1] = 255.0f;
                        qp.Enqueue(new Point(point.X, point.Y + 1));
                    }


                    if (point.Y < nms.Height - 1 && point.X < nms.Width - 1 && res[point.X + 1, point.Y + 1] == 128.0f && nms[point.X + 1, point.Y + 1] > thr_weak)
                    {
                        res[point.X+1, point.Y + 1] = 255.0f;
                        qp.Enqueue(new Point(point.X+1, point.Y + 1));
                    }
                    if (point.Y < nms.Height - 1 && point.X > 0 && res[point.X-1, point.Y + 1] == 128.0f && nms[point.X-1, point.Y + 1] > thr_weak)
                    {
                        res[point.X-1, point.Y + 1] = 255.0f;
                        qp.Enqueue(new Point(point.X-1, point.Y + 1));
                    }
                    if (point.X < nms.Width - 1 && point.Y > 0 && res[point.X + 1, point.Y - 1] == 128.0f && nms[point.X+1, point.Y - 1] > thr_weak)
                    {
                        res[point.X+1, point.Y - 1] = 255.0f;
                        qp.Enqueue(new Point(point.X+1, point.Y - 1));
                    }
                    if (point.Y > 0 && point.X > 0 && res[point.X-1, point.Y - 1] == 128.0f && nms[point.X-1, point.Y - 1] > thr_weak)
                    {
                        res[point.X -1, point.Y - 1] = 255.0f;
                        qp.Enqueue(new Point(point.X -1, point.Y - 1));
                    }


                }
               // ImageIO.ImageToFile(res, "cannyTest/TestTresholdmmm.png"); 

                for (int x = 0; x < res.Width; x++)
                {
                    for (int y = 0; y < res.Height; y++)
                    {
                        if (res[x, y] !=255)
                        {
                            res[x, y] = 0.0f;
                        }
                    }
                }
               // ImageIO.ImageToFile(res, "cannyTest/TestTresholdmmmFINAL.png"); 
            return res;
        }

        static Point[] getAnglePoints(GrayscaleFloatImage img, GrayscaleFloatImage canny)
        {
            //Create the accu 
            double scale_r = 1;
            int h = img.Width;
            int w = img.Height;
            int hough_h = ((int)(Math.Sqrt(2.0) * (double)(h > w ? h : w)));///2
            int _accu_h = (int)(hough_h * 2 * scale_r + 1); // -r -> +r //*2 == увеличиваем расширение по р 
            int _accu_w = 360;//180 + 1;
            List<Point> angelAndCount = new List<Point>();


            //float[,] _accu = new float[_accu_h , _accu_w];
            //float[,] _accu_copy = new float[_accu_h , _accu_w];

            GrayscaleFloatImage _accu = new GrayscaleFloatImage(_accu_h, _accu_w);
            List<Point>[,] accuCoordinatesList = new List<Point>[_accu_h, _accu_w];
            for (int r = 0; r < _accu_h; r++)
            {
                for (int t = 0; t < _accu_w; t++)
                {
                    accuCoordinatesList[r, t] = new List<Point>();
                }
            }

            int center_x = w / 2;
            int center_y = h / 2;




            double DEG2RAD = ((Math.PI / 180.0)) * (180.0 / _accu_w);
            for (int y = 0; y < h; y++)
            {
                for (int x = 0; x < w; x++)
                {
                    if (canny[(y), x] >= 128)
                    {
                        //Console.WriteLine( "x= "+ x+ " y=  " + y );
                        for (int t = 0; t < _accu_w; t++)
                        {
                            double r = (((double)x) * Math.Cos((double)t * DEG2RAD)) + (((double)y) * Math.Sin((double)t * DEG2RAD));
                            int indexX = (int)((Math.Round(r + hough_h)) * scale_r);
                            int indexY = t;
                            _accu[indexX, indexY]++;
                            // _accu_copy[(int)((Math.Round(r + hough_h))) , t]++;

                            //Rd = (int)r;
                            //R = (int)Math.Round(r + hough_h);
                            accuCoordinatesList[indexX, indexY].Add(new Point(x, y));
                            //cout << "t= " << t << " R=  " << R << endl;
                        }
                    }
                }
            }

            for (int r = 0; r < _accu_h; r++)
            {
                for (int t = 0; t < _accu_w; t++)
                {
                    accuCoordinatesList[r, t].OrderBy(p => p.X).ThenBy(p => p.Y);
                }
            }

            float maxa = 0;
            for (int r = 0; r < _accu_h; r++)
            {
                for (int t = 0; t < _accu_w; t++)
                {
                    if ((int)_accu[r, t] > maxa)
                        maxa = _accu[r, t];
                }
            }
            float p1=0, p2=0, p3=0, p4=0;
            int r1=0, r2=0, r3=0, r4=0;
            int t1=0, t2=0, t3=0, t4=0;
            for (int r = 0; r < _accu_h; r++)
            {
                for (int t = 0; t < _accu_w; t++)
                {

                    float max = _accu[r, t];
                    for (int ly = -4; ly <= 4; ly++)
                    {
                        for (int lx = -4; lx <= 4; lx++)
                        {
                            if ((ly + r >= 0 && ly + r < _accu_h) && (lx + t >= 0 && lx + t < _accu_w))
                            {
                                if ((int)_accu[((r + ly)), (t + lx)] > max)
                                {
                                    max = _accu[((r + ly)), (t + lx)];
                                    ly = lx = 5;
                                }
                            }
                        }
                    }
                    if (max > (int)_accu[(r), t])
                        continue;
                    if (p1 < _accu[r, t])
                    {
                        p1 = _accu[r, t];
                        r1 = r;
                        t1 = t;
                    }
                    else if (p2 < _accu[r, t])
                    {
                        p2 = _accu[r, t];
                        r2 = r;
                        t2 = t;
                    }
                    //else if (p3 < _accu[r, t])
                    //{
                    //    p3 = _accu[r, t];
                    //    r3 = r;
                    //    t3 = t;
                    //}
                    //else if (p4 < _accu[r, t])
                    //{
                    //    p4 = _accu[r, t];
                    //    r4 = r;
                    //    t4 = t;
                    //}
                    
                }
            }


            int rTransp = 0;
            if (t1<180)
            {
                rTransp = t1 + 180;
            }
            else
            {
                rTransp = t1 - 180;
            }

            for (int r = 0; r < _accu_h; r++)
            {
                for (int t = rTransp-8; t < rTransp+8; t++)
                {

                    float max = _accu[r, t];
                    for (int ly = -4; ly <= 4; ly++)
                    {
                        for (int lx = -4; lx <= 4; lx++)
                        {
                            if ((ly + r >= 0 && ly + r < _accu_h) && (lx + t >= 0 && lx + t < _accu_w))
                            {
                                if ((int)_accu[((r + ly)), (t + lx)] > max)
                                {
                                    max = _accu[((r + ly)), (t + lx)];
                                    ly = lx = 5;
                                }
                            }
                        }
                    }
                    if (max > (int)_accu[(r), t])
                        continue;
                    //if (p1 < _accu[r, t])
                    //{
                    //    p1 = _accu[r, t];
                    //    r1 = r;
                    //    t1 = t;
                    //}
                    //else if (p2 < _accu[r, t])
                    //{
                    //    p2 = _accu[r, t];
                    //    r2 = r;
                    //    t2 = t;
                    //}
                    //else 
                    if (p3 < _accu[r, t])
                    {
                        p3 = _accu[r, t];
                        r3 = r;
                        t3 = t;
                    }
                    else if (p4 < _accu[r, t])
                    {
                        p4 = _accu[r, t];
                        r4 = r;
                        t4 = t;
                    }

                }
            }





            Console.WriteLine("point_max" + t1 + "  " + t2 + "  " + t3 + "  " + t4);
            ColorFloatImage imageColor = new ColorFloatImage(img.Width, img.Height);
            foreach (var OnePoint in accuCoordinatesList[r1, t1])
            {
                ColorFloatPixel pointColor = new ColorFloatPixel();
                pointColor.a = 256;
                pointColor.b = 0;
                pointColor.r = 256;
                pointColor.g = 0;

                imageColor[(int)(OnePoint.Y * scale_r), (int)(OnePoint.X * scale_r)] = pointColor;
                //ImageIO.ImageToFile(imageColor, "test/TMP.png");


            }
            foreach (var OnePoint in accuCoordinatesList[r2, t2])
            {
                ColorFloatPixel pointColor = new ColorFloatPixel();
                pointColor.a = 256;
                pointColor.b = 0;
                pointColor.r = 256;
                pointColor.g = 0;

                imageColor[(int)(OnePoint.Y * scale_r), (int)(OnePoint.X * scale_r)] = pointColor;
                //ImageIO.ImageToFile(imageColor, "test/TMP.png");


            }
            foreach (var OnePoint in accuCoordinatesList[r3, t3])
            {
                ColorFloatPixel pointColor = new ColorFloatPixel();
                pointColor.a = 256;
                pointColor.b = 0;
                pointColor.r = 256;
                pointColor.g = 0;

                imageColor[(int)(OnePoint.Y * scale_r), (int)(OnePoint.X * scale_r)] = pointColor;
                //ImageIO.ImageToFile(imageColor, "test/TMP.png");


            }
            foreach (var OnePoint in accuCoordinatesList[r4, t4])
            {
                ColorFloatPixel pointColor = new ColorFloatPixel();
                pointColor.a = 256;
                pointColor.b = 0;
                pointColor.r = 256;
                pointColor.g = 0;

                imageColor[(int)(OnePoint.Y * scale_r), (int)(OnePoint.X * scale_r)] = pointColor;
                


            }
            ImageIO.ImageToFile(imageColor, "test/TMP.png");
            Console.WriteLine("max Accu = " + maxa);
            double contrast = 1.0;
            double coef = 255.0 / (double)maxa * contrast;

            for (int r = 0; r < _accu_h; r++)
            {
                for (int t = 0; t < _accu_w; t++)
                {
                    double c = ((_accu[r, t] * coef) < 255.0) ? (double)_accu[r, t] * coef : 255.0;
                    _accu[r, t] = (float)(255 - c);
                    //_accu[r, t] = (float)( c);
                }
            }
            GrayscaleFloatImage retPoint = new GrayscaleFloatImage(img.Width, img.Height);
            Point[] arrayOfPoints = new Point[4];
            
            
            
            int i=accuCoordinatesList[r1, t1][0].X;
            int j=accuCoordinatesList[r1, t1][0].Y;
            retPoint[j, i] = 255;
            arrayOfPoints[0] = new Point(j, i);
            
            int cap = accuCoordinatesList[r1, t1].Count;
            i = accuCoordinatesList[r1, t1][cap-1].X;
            j = accuCoordinatesList[r1, t1][cap-1].Y;
            retPoint[j, i] = 255;
            arrayOfPoints[1] = new Point(j, i);

            i = accuCoordinatesList[r2, t2][0].X;
            j = accuCoordinatesList[r2, t2][0].Y;
            retPoint[j, i] = 255;
            arrayOfPoints[2] = new Point(j, i);

            cap = accuCoordinatesList[r2, t2].Count;
            i = accuCoordinatesList[r2, t2][cap - 1].X;
            j = accuCoordinatesList[r2, t2][cap - 1].Y;
            retPoint[j, i] = 255;
            arrayOfPoints[3] = new Point(j, i);
            
            
            ImageIO.ImageToFile(_accu, "test/Accu.bmp");
            ImageIO.ImageToFile(retPoint, "test/TMP2.bmp");


        
            return arrayOfPoints;

        }

        static int det(int a1, int a2, int b1, int b2)
        {
            return a1 * b2 - a2 * b1;
        }

        static Point[] getPointFromCanny(GrayscaleFloatImage img, GrayscaleFloatImage canny)
        {
            GrayscaleFloatImage retPoint = new GrayscaleFloatImage(img.Width, img.Height);
            Point[] arrayOfPoints = new Point[4];

            int colum = img.Width; // количество столбцов массива
            int row = img.Height;
            int x = 0, y = 0;  // Координаты текущего элемента массива
             // значение, которым заполняется массив
            bool endOfDiag = false;
            // зполнение первой половины массива по диагонали, зигзагом, начиная
            // слева и сверху, заканчивая  побочной диагональю
            for (int diag = 0; diag < colum; diag++) // выполняем проход по диагоналям
            {
                if (endOfDiag) break;
                if (diag % 2 == 0) // по четным диагоналям
                {
                    x = 0; // х-координата первого лемента массива на диагонали - diag
                    y = diag; // у-координата элемента массива на диагонали - diag

                    while (y >= 0) // пока y-координата находится в верхней части диагонали
                    {
                        if (canny[x, y] == 255)
                        { endOfDiag = true; break; }// записать значение в массив
                        img[x, y] = 255;
                        x++;     // по горизонтали, смещаемся влево
                        y--;    // по вертикали, смещаемся вниз
                    }
                }
                else // по нечетным диагоналям
                {
                    x = diag; // х-координата элемента массива на диагонали - diag
                    y = 0; // у-координата первого элемента массива на диагонали - diag

                    while (x >= 0) // пока x-координата находится в левой части диагонали
                    {
                        if (canny[x, y] == 255)
                        { endOfDiag = true; break; } // записать значение в массив
                        img[x, y] = 255;

                        x -= 1;  // по горизонтали, смещаемся вправо
                        y += 1; // по вертикали, смещаемся вверх
                    }
                }
            }

            arrayOfPoints[0] = new Point(x, y);
            // конец for
            endOfDiag = false;

            for (int diag = img.Height; diag > 0 ; diag--) // выполняем проход по диагоналям
            {
                if (endOfDiag) break;
                 // по четным диагоналям
                
                    x = 0; // х-координата первого лемента массива на диагонали - diag
                    y = diag-1; // у-координата элемента массива на диагонали - diag

                    while ((y<img.Height)) // пока y-координата находится в верхней части диагонали
                    {
                        if (canny[x, y] == 255)
                        { endOfDiag = true; break; }// записать значение в массив
                        img[x, y] = 255;
                        x++;     // по горизонтали, смещаемся влево
                        y++;    // по вертикали, смещаемся вниз
                    }
                                
            }

            arrayOfPoints[1] = new Point(x, y);


            endOfDiag = false;

            for (int diag = img.Height; diag > 0; diag--) // выполняем проход по диагоналям
            {
                if (endOfDiag) break;
                // по четным диагоналям

                x = img.Width-1; // х-координата первого лемента массива на диагонали - diag
                y = diag - 1; // у-координата элемента массива на диагонали - diag

                while ((y < img.Height)) // пока y-координата находится в верхней части диагонали
                {
                    if (canny[x, y] == 255)
                    { endOfDiag = true; break; }// записать значение в массив
                    img[x, y] = 255;
                    x--;     // по горизонтали, смещаемся влево
                    y++;    // по вертикали, смещаемся вниз
                }

            }
            arrayOfPoints[2] = new Point(x, y);






            endOfDiag = false;

            for (int diag = 0; diag < img.Height; diag++) // выполняем проход по диагоналям
            {
                if (endOfDiag) break;
                // по четным диагоналям

                x = img.Width - 1; // х-координата первого лемента массива на диагонали - diag
                y = diag ; // у-координата элемента массива на диагонали - diag

                while ((y > 0)) // пока y-координата находится в верхней части диагонали
                {
                    if (canny[x, y] == 255)
                    { endOfDiag = true; break; }// записать значение в массив
                    img[x, y] = 255;
                    x--;     // по горизонтали, смещаемся влево
                    y--;    // по вертикали, смещаемся вниз
                }

            }
            arrayOfPoints[3] = new Point(x, y);






            //for (int i = 0; i < img.Width;i++ )
            //{
            //    for (int j=0;j<img.Height;j++)
            //    {
            //        if (canny[i,j]==255)
            //        {

            //        }

            //    }
            //}
                //int h_mid = img.Height / 2;
                //int w_mid = img.Width / 2;
                //bool end=false;
                //int i=0;
                //int x1 = 0, y1 = 0;
                //int x2 = 0, y2 = 0;
                //int x5 = 0, y5 = 0;
                //int x6 = 0, y6 = 0;
                //int x7 = 0, y7 = 0;
                //int x8 = 0, y8 = 0;
                //i = 0;
                //int step = 5;
                //while (i<img.Height)
                //{
                //    if (canny[w_mid-step, i]==255)
                //    {
                //        x1 = w_mid - step;
                //        y1 = i;
                //        break;
                //    }
                //    i++;
                //}

                //i = 0;
                //while (i < img.Height)
                //{
                //    if (canny[w_mid+step, i] == 255)
                //    {
                //        x2 = w_mid + step;
                //        y2 = i;
                //        break;
                //    }
                //    i++;
                //}

                //int x3 = 0, y3 = 0;
                //int x4 = 0, y4 = 0;
                //i = 0;
                //while (i < img.Width)
                //{
                //    if (canny[i, h_mid-step] == 255)
                //    {
                //        x3 = i;
                //        y3 = h_mid - step;
                //        break;
                //    }
                //    i++;
                //}
                //end = false;
                //i = 0;
                //while (i < img.Width)
                //{
                //    if (canny[i, h_mid+step] == 255)
                //    {
                //        x4 = i;
                //        y4 = h_mid + step;
                //        break;
                //    }
                //    i++;
                //}



                //int a_l1, b_l1, c_l1;
                //a_l1 = y2 - y1;//y1 - y2;
                //b_l1 = x1 - x2;//x2 - x1;
                //c_l1 = -(x1 * y2 - y1 * x2);// -((x2 - x1) * y1 - (y2 - y1) * x1); //( x1 * y2 -  y1 * x2);

                //int a_l2, b_l2, c_l2;
                //a_l2 = y4 - y3;
                //b_l2 = x3 - x4;
                //c_l2 = -(x3 * y4 - y3 * x4);//-((x4 - x3) * y3 - (y4 - y3) * x3);//-( x3 * y4 - y3 * x4);

                //int d, d1, d2;
                //d = det(a_l1, b_l1, a_l2, b_l2);
                //d1 = -det(c_l1, b_l1, c_l2, b_l2);
                //d2 = -det(a_l1, c_l1, a_l2, c_l2);

                //arrayOfPoints[0] = new Point((int)((double)d1 / d), (int)((double)d2 / d));

                //i = img.Height-1;
                //while (i >0)
                //{
                //    if (canny[w_mid - step, i] == 255)
                //    {
                //        x5 = w_mid - step;
                //        y5 = i;
                //        break;
                //    }
                //    i--;
                //}

                //i = img.Height-1;
                //while (i < img.Height)
                //{
                //    if (canny[w_mid + step, i] == 255)
                //    {
                //        x6 = w_mid + step;
                //        y6 = i;
                //        break;
                //    }
                //    i--;
                //}

                //i = img.Width-1;
                //while (i < img.Width)
                //{
                //    if (canny[i, h_mid - step] == 255)
                //    {
                //        x7 = i;
                //        y7 = h_mid - step;
                //        break;
                //    }
                //    i--;
                //}
                //end = false;
                //i = img.Width-1;
                //while (i < img.Width)
                //{
                //    if (canny[i, h_mid + step] == 255)
                //    {
                //        x8 = i;
                //        y8 = h_mid + step;
                //        break;
                //    }
                //    i--;
                //}


                //int a_l3, b_l3, c_l3;
                //a_l3 = y6 - y5;
                //b_l3 = x5 - x6;
                //c_l3 = -(x5 * y6 - y5 * x6);//-((x6 - x5) * y5 - (y6 - y5) * x5);//-( x3 * y4 - y3 * x4);

                //d = det(a_l3, b_l3, a_l2, b_l2);
                //d1 = -det(c_l3, b_l3, c_l2, b_l2);
                //d2 = -det(a_l3, c_l3, a_l2, c_l2);

                //arrayOfPoints[1] = new Point((int)((double)d1 / d), (int)((double)d2 / d));

                ////int a_l4, b_l4, c_l4;
                //int a_l4 = y8 - y7;
                //int b_l4 = x7 - x8;
                //int c_l4 = -(x7 * y8 - y7 * x8);//-((x8 - x7) * y7 - (y8 - y7) * x7);//-((x4 - x3) * y3 - (y4 - y3) * x3);

                //d = det(a_l3, b_l3, a_l4, b_l4);
                //d1 = -det(c_l3, b_l3, c_l4, b_l4);
                //d2 = -det(a_l3, c_l3, a_l4, c_l4);

                //arrayOfPoints[2] = new Point((int)(Math.Floor((double)d1 / d)), (int)Math.Floor( (double)d2 / d ));

                //d = det(a_l1, b_l1, a_l4, b_l4);
                //d1 = -det(c_l1, b_l1, c_l4, b_l4);
                //d2 = -det(a_l1, c_l1, a_l4, c_l4);

                //arrayOfPoints[3] = new Point((int)((double)d1 / d), (int)((double)d2 / d));

            retPoint[arrayOfPoints[0].X, arrayOfPoints[0].Y] = 255;
            retPoint[arrayOfPoints[1].X, arrayOfPoints[1].Y] = 255;
            retPoint[arrayOfPoints[2].X, arrayOfPoints[2].Y] = 255;
            retPoint[arrayOfPoints[3].X, arrayOfPoints[3].Y] = 255;
            //ImageIO.ImageToFile(canny, "test/TMP1.bmp");
            //ImageIO.ImageToFile(retPoint, "test/TMP2.bmp");
            //ImageIO.ImageToFile(img, "test/TMP3.bmp");

            return arrayOfPoints;
        }

        static ColorFloatImage getRectangle(ColorFloatImage img1Color, Point[] points, int width,int height, bool isFirst)
        {
            //List<Point> sorted = points.OrderBy(p => p.X).ThenBy(p => p.X).ToList();
            List<Point> sorted = points.OrderBy(p => p.X).ToList();

            //for (int i = 0; i < sorted.Count; i++)
            //{
            //    Console.WriteLine(points[i].ToString() + "        " + sorted[i].ToString());
            //    //Console.WriteLine(sorted[i].ToString());
            //}

            Point leftHight = new Point(sorted[0].X, sorted[0].Y);
            Point leftDown = new Point(sorted[1].X, sorted[1].Y);
            Point rightHight = new Point(sorted[2].X, sorted[2].Y);
            Point rightDown = new Point(sorted[3].X, sorted[3].Y);
            if (sorted[0].Y>sorted[1].Y)
            {
                leftHight = new Point(sorted[1].X, sorted[1].Y);
                leftDown = new Point(sorted[0].X, sorted[0].Y);
            }
            if (sorted[2].Y > sorted[3].Y)
            {
                rightHight = new Point(sorted[3].X, sorted[3].Y);
                rightDown = new Point(sorted[2].X, sorted[2].Y);
            }

            //Console.WriteLine(leftHight + "------------" + rightHight);
            //Console.WriteLine(leftDown + "------------" + rightDown);

            int finalNormSizeW = width;
            int finalNormSizeH = height;

            if (isFirst)
            {
                //Getting size
                double widthHight = Math.Sqrt(Math.Pow((leftHight.X - rightHight.X), 2) + Math.Pow((leftHight.Y - rightHight.Y), 2));
                double widthDown = Math.Sqrt(Math.Pow((leftDown.X - rightDown.X), 2) + Math.Pow((leftDown.Y - rightDown.Y), 2));

                double heightLeft = Math.Sqrt(Math.Pow((leftHight.X - leftDown.X), 2) + Math.Pow((leftHight.Y - leftDown.Y), 2));
                double heightRight = Math.Sqrt(Math.Pow((rightHight.X - rightDown.X), 2) + Math.Pow((rightHight.Y - rightDown.Y), 2));

                finalNormSizeW = (int)(Math.Max(widthHight, widthDown));
                finalNormSizeH = (int)(Math.Max(heightLeft, heightRight));

                //Console.Write("f");

                //ColorFloatImage normalImg = new ColorFloatImage((int)(Math.Max(widthHight, widthDown)), (int)(Math.Max(heightLeft, heightRight)));
            }
            
            ColorFloatImage normalImg = new ColorFloatImage(finalNormSizeW, finalNormSizeH);
         


            for (int x = 0; x < normalImg.Height; x++)
            {
                for (int y = 0; y < normalImg.Width; y++)
                {


                    float lambdaWidth = (float)y / (normalImg.Width - y); //здесь олучаем пропорцию
                    float lambdaHeight = (float)x / (normalImg.Height - x);//и здесь....только как ?
                    if (x == normalImg.Height)
                    {
                        lambdaWidth = normalImg.Width;
                        continue;
                    }
                    if (y == normalImg.Width)
                    {
                        lambdaHeight = normalImg.Height;
                        continue;
                    }

                    float x_hight = (float)(leftHight.X + lambdaWidth * rightHight.X) / (1 + lambdaWidth);
                    float y_hight = (float)(leftHight.Y + lambdaWidth * rightHight.Y) / (1 + lambdaWidth);

                    float x_down = (float)(leftDown.X + lambdaWidth * rightDown.X) / (1 + lambdaWidth);
                    float y_down = (float)(leftDown.Y + lambdaWidth * rightDown.Y) / (1 + lambdaWidth);


                    float x_left = (float)(leftHight.X + lambdaHeight * leftDown.X) / (1 + lambdaHeight);
                    float y_left = (float)(leftHight.Y + lambdaHeight * leftDown.Y) / (1 + lambdaHeight);

                    float x_right = (float)(rightHight.X + lambdaHeight * rightDown.X) / (1 + lambdaHeight);
                    float y_right = (float)(rightHight.Y + lambdaHeight * rightDown.Y) / (1 + lambdaHeight);

                    float u = ((x_right - x_left) * (y_down - y_left) - (y_right - y_left) * (x_down - x_left)) /
                        ((y_right - y_left) * (x_hight - x_down) - (x_right - x_left) * (y_hight - y_down));

                    float x0 = x_down + u * (x_hight - x_down);
                    float y0 = y_down + u * (y_hight - y_down);


                    //float x0 = (float)y / (normalImg.Width - y); //Red.GetLength(1) / Red1.GetLength(1);
                    //float y0 = (float)x / (normalImg.Height - x); //Red.GetLength(0) / Red1.GetLength(0);


                    int y1 = 0, y2 = 0;
                    int x1 = 0, x2 = 0;
                    x1 = (int)Math.Floor((double)x0);
                    x2 = (int)Math.Min(x1 + 1, img1Color.Width - 1);
                    y1 = (int)Math.Floor((double)y0);
                    y2 = (int)Math.Min(y1 + 1, img1Color.Height - 1);
                    int[] rows = new int[2] { y1, y2 };
                    int[] cols = new int[2] { x1, x2 };

                    if (x2 >= img1Color.Width)
                    {
                        x2 = img1Color.Width - 1;
                    }
                    if (y2 >= img1Color.Height)
                    {
                        y2 = img1Color.Height - 1;
                    }
                    if (x1 < 0)
                    {
                        x1 = 0;
                    }
                    if (y1 < 0)
                    {
                        y1 = 0;
                    }

                    float[,] matrix = new float[2, 2] 
                        {
                            {1, -1},
                            {-1, 1}
                        };
                    ColorFloatPixel pixel = new ColorFloatPixel();

                    for (int i = 0; i < 2; ++i)
                    {
                        for (int j = 0; j < 2; ++j)
                        {
                            float kernel = matrix[i, j];
                            for (int k = 0; k < 2; ++k)
                            {
                                if (k != i)
                                    kernel *= (y0 - (int)(y0) - k);
                                if (k != j)
                                    kernel *= (x0 - (int)(x0) - k);
                            }

                            pixel.r += kernel * img1Color[cols[j], rows[i]].r;
                            pixel.g += kernel * img1Color[cols[j], rows[i]].g;
                            pixel.b += kernel * img1Color[cols[j], rows[i]].b;



                        }


                    }
                    normalImg[y, x] = pixel;

                }
            }
            return normalImg;
        }

        static void HightlightDe(string [] names, string outputFileName)
        {
            GrayscaleFloatImage img1= ImageIO.FileToGrayscaleFloatImage(names[0]);
            //GrayscaleFloatImage img2_c= ImageIO.FileToGrayscaleFloatImage(image2Name);
            GrayscaleFloatImage canny1 = cannyGetImg(img1, 5, 0.1f, 0.7f);
            //GrayscaleFloatImage canny2 = cannyGetImg(img2_c, 5, 0.01f, 0.3f);
            ColorFloatImage img1Color = ImageIO.FileToColorFloatImage(names[0]);
            //ColorFloatImage img1Color_onlyForme = ImageIO.FileToColorFloatImage(image1Name);


            //Point[] points = getAnglePoints(img1, canny1);

            Point[] points = getPointFromCanny(img1, canny1);

            List<ColorFloatImage> listOfGoodImages = new List<ColorFloatImage>();

            listOfGoodImages.Add(getRectangle(img1Color, points, 0, 0, true));
            
            ////
            for (int i=1;i<names.Length;i++)
            {
                if (!File.Exists(names[i]))
                {
                    Console.WriteLine("Some File doesn't exist, sorry");
                    return;
                }
                GrayscaleFloatImage img_i = ImageIO.FileToGrayscaleFloatImage(names[i]);
                GrayscaleFloatImage canny_i = cannyGetImg(img_i, 5, 0.1f, 0.7f);
                ColorFloatImage img1Color_i = ImageIO.FileToColorFloatImage(names[i]);
                Point[] points_i = getPointFromCanny(img_i, canny_i);
                listOfGoodImages.Add(getRectangle(img1Color_i, points_i, listOfGoodImages[0].Width, listOfGoodImages[0].Height, false));

            }
            ColorFloatImage finaleDe = new ColorFloatImage(listOfGoodImages[0].Width, listOfGoodImages[0].Height);
            for (int x=0;x<listOfGoodImages[0].Width;x++)
            {
                for (int y=0;y<listOfGoodImages[0].Height;y++)
                {
                    ColorFloatPixel p = new ColorFloatPixel();
                    float r = 0;
                    float g = 0;
                    float b = 0;
                    float theshold= 240;
                    if (listOfGoodImages[0][x, y].r > theshold && listOfGoodImages[0][x, y].g > theshold && listOfGoodImages[0][x, y].b > theshold)
                    {
                        

                        int numberOfImages = listOfGoodImages.Count;
                        for (int n = 1; n < numberOfImages; n++)
                        {
                            r += listOfGoodImages[n][x, y].r;
                            g += listOfGoodImages[n][x, y].g;
                            b += listOfGoodImages[n][x, y].b;
                        }
                        p.r = r / (numberOfImages-1);
                        p.g = g / (numberOfImages-1);
                        p.b = b / (numberOfImages-1);

                        
                    }
                    else
                    {
                        p.r = listOfGoodImages[0][x, y].r;
                        p.g = listOfGoodImages[0][x, y].g;
                        p.b = listOfGoodImages[0][x, y].b;
                    }
                    finaleDe[x, y] = p;
                }
            }
            for (int n = 0; n < listOfGoodImages.Count; n++)
            {
                ImageIO.ImageToFile(listOfGoodImages[n], "test/TMP_" + n + ".jpg");
            }

            ImageIO.ImageToFile(finaleDe, outputFileName);
            
        }

        static void Main(string[] args)
        {
            try
            {
                if (args.Length < 3)
                {
                    Console.WriteLine("Not enough arguments.");
                    Console.ReadLine();
                    return;
                }
                string InputFileName = args[0];
                string OutputFileName = args[1];
                string command = args[2];
                if (!File.Exists(InputFileName))
                {
                    Console.WriteLine(" File doesn't exist");
                    Console.ReadLine();
                    return;
                }
                ColorFloatImage image = ImageIO.FileToColorFloatImage(InputFileName);
                GrayscaleFloatImage imageGray = ImageIO.FileToGrayscaleFloatImage(InputFileName);
                switch (command)
                {
                    case "invert":
                        if (args.Length > 3)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        InvertImage(image);
                        ImageIO.ImageToFile(image, OutputFileName);
                        break;
                    case "eqhist":
                        if (args.Length > 3)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        eqhist(imageGray);
                        ImageIO.ImageToFile(imageGray, OutputFileName);
                        break;
                    case "mirror":// {x|y}"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        switch (args[3])
                        {
                            case "x":
                                MirrowXImage(image);
                                break;
                            case "y":
                                MirrowYImage(image);
                                break;
                            default:
                                Console.WriteLine("Wrong parametr.");
                                return;
                                break;

                        }
                        ImageIO.ImageToFile(image, OutputFileName);

                        break;
                    case "rotate":// {cw|ccw} (angle)"
                        if (args.Length < 5)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 5)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        float number;
                        if (!float.TryParse(args[4], out number))
                        {
                            Console.WriteLine("Strange number");
                        }
                        number = float.Parse(args[4]);
                        switch (args[3])
                        {

                            case "cw":
                                rotateImage(image, true, number, OutputFileName);
                                break;
                            case "ccw":
                                rotateImage(image, true, (-1 * number), OutputFileName);
                                break;
                            default:
                                Console.WriteLine("Wrong parametr.");
                                return;
                                break;

                        }

                        break;
                    case "prewitt":// {x|y}"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        switch (args[3])
                        {
                            case "x":
                                PrewittImage(imageGray, true);
                                break;
                            case "y":
                                PrewittImage(imageGray, false);
                                break;
                            default:
                                Console.WriteLine("Wrong parametr.");
                                return;
                                break;

                        }
                        ImageIO.ImageToFile(imageGray, OutputFileName);
                        break;
                    case "sobel":// {x|y}"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        switch (args[3])
                        {
                            case "x":
                                SobelImage(imageGray, true);
                                break;
                            case "y":
                                SobelImage(imageGray, false);
                                break;
                            default:
                                Console.WriteLine("Wrong parametr.");
                                return;
                                break;

                        }
                        ImageIO.ImageToFile(imageGray, OutputFileName);
                        break;
                    case "roberts":// {1|2}":
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        switch (args[3])
                        {
                            case "1":
                                RobertsImage(imageGray, true);
                                break;
                            case "2":
                                RobertsImage(imageGray, false);
                                break;
                            default:
                                Console.WriteLine("Wrong parametr.");
                                return;
                                break;

                        }
                        ImageIO.ImageToFile(imageGray, OutputFileName);
                        break;
                    case "median":// (rad)"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        int rad;
                        if (!int.TryParse(args[3], out rad))
                        {
                            Console.WriteLine("Strange number");
                        }
                        rad = int.Parse(args[3]);
                        MedianImage(image, rad);
                        ImageIO.ImageToFile(image, OutputFileName);
                        break;
                    case "gauss":// (sigma)"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        float sigma;
                        if (!float.TryParse(args[3], out sigma))
                        {
                            Console.WriteLine("Strange number");
                        }
                        sigma = float.Parse(args[3]);
                        GaussImage(image, (int)(3 * sigma), sigma);
                        ImageIO.ImageToFile(image, OutputFileName);
                        break;
                    case "gradient":// (sigma)"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        float sigma_grad;
                        if (!float.TryParse(args[3], out sigma_grad))
                        {
                            Console.WriteLine("Strange number");
                        }
                        sigma_grad = float.Parse(args[3]);
                        GradientFromGaussImage(imageGray, (int)(3 * sigma_grad), sigma_grad);
                        ImageIO.ImageToFile(imageGray, OutputFileName);
                        break;
                    case "up_billinear":// (rad)"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        float number_up;
                        if (!float.TryParse(args[3], out number_up))
                        {
                            Console.WriteLine("Strange number");
                        }
                        number_up = float.Parse(args[3]);
                        BillinearInterpolation(image, number_up, OutputFileName);
                        break;
                    case "up_bicubic":// (rad)"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        float number_up_bic;
                        if (!float.TryParse(args[3], out number_up_bic))
                        {
                            Console.WriteLine("Strange number");
                        }
                        number_up_bic = float.Parse(args[3]);
                        BicubicInterpolation(image, number_up_bic, OutputFileName);
                        break;
                    case "downsample":// (rad)"
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 4)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        float number_down;
                        if (!float.TryParse(args[3], out number_down))
                        {
                            Console.WriteLine("Strange number");
                        }
                        number_down = float.Parse(args[3]);
                        float sigmaDown = (float)Math.Sqrt(0.1f * (number_down * number_down - 1));
                        GaussImage(image, (int)(3 * sigmaDown), sigmaDown);

                        DownSample(image, number_down, OutputFileName);
                        break;
                    case "metric":
                        if (args.Length < 4)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (!File.Exists(OutputFileName))
                        {
                            Console.WriteLine(" File doesn't exist");
                            Console.ReadLine();
                            return;
                        }
                        GrayscaleFloatImage image2 = ImageIO.FileToGrayscaleFloatImage(OutputFileName);

                        switch (args[3])
                        {
                            case "mse":
                                Console.WriteLine("MSE= " + mse_metric(imageGray, image2));
                                break;
                            case "psnr":
                                Console.WriteLine("PSNR= " + psnr_metric(imageGray, image2));
                                break;
                            case "ssim":
                                Console.WriteLine("SSIM= " + ssim_metric(imageGray, image2));
                                break;
                            case "mssim":
                                Console.WriteLine("MSSIM= " + mssim_metric(imageGray, image2));
                                break;
                            default:
                                Console.WriteLine("Wrong parametr.");
                                return;
                                break;

                        }
                        break;
                    case "dcci":
                        if (args.Length > 3)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        dcci(image, 2, OutputFileName);
                        //ImageIO.ImageToFile(imageGray, OutputFileName);
                        break;
                    case "canny":// (rad)"
                        if (args.Length < 6)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        if (args.Length > 6)
                        {
                            Console.WriteLine("Too many arguments.");
                            Console.ReadLine();
                            return;
                        }
                        float sigmaCanny;
                        float h1;
                        float h2;
                        if ((!float.TryParse(args[3], out sigmaCanny)) || ((!float.TryParse(args[4], out h1)) || ((!float.TryParse(args[5], out h2)))))
                        {
                            Console.WriteLine("Strange number");
                        }
                        sigmaCanny = float.Parse(args[3]);
                        h1 = float.Parse(args[4]);
                        h2 = float.Parse(args[5]);
                        Canny(imageGray, sigmaCanny, h1, h2, OutputFileName);
                        break;
                    case "dehighlight":
                        if (args.Length < 3)
                        {
                            Console.WriteLine("Not enough arguments.");
                            Console.ReadLine();
                            return;
                        }
                        List<string> listOfNames = new List<string>();
                        listOfNames.Add(InputFileName);
                        for (int i=3;i<args.Length;i++)
                        {
                            listOfNames.Add(args[i]);
                        }
                        try
                        {
                            HightlightDe(listOfNames.ToArray(), OutputFileName);
                        }
                        catch (Exception ex)
                        {
                            Console.WriteLine("Something went wrong...Maybe input images was not good enough for Canny with parameters 5, 0.01f, 0.3f");
                        }

                        break;
                    default:
                        Console.WriteLine("Wrong command");
                        return;
                        break;

                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Something went wrong...");
            }


            // Console.WriteLine("Not enough arguments.");
             
            

            //string parametrs=args[3];
            //if (!File.Exists(InputFileName))
            //{
            //    Console.WriteLine("Doesn't exist");
            //    Console.ReadLine();
            //    return;
            //}

            //ColorFloatImage image = ImageIO.FileToColorFloatImage(InputFileName);
            //GrayscaleFloatImage imageGray = ImageIO.FileToGrayscaleFloatImage(InputFileName);

            



            //rotateImage(image,true,30);
            //PrewittImage(imageGray, true);
            //SobelImage(imageGray, true);
            //RobertsImage(imageGray, true);
            //MedianImage(image, 5);
            //GaussImage(image, 5, 3.5f);
            //GradientFromGaussImage(imageGray, 1, 1);
            //eqhist(imageGray);


            //ImageIO.ImageToFile(imageGray, OutputFileName);
            //ImageIO.ImageToFile(image, OutputFileName);

            Console.WriteLine("Done");
            Console.ReadLine();

        }
    }
}
