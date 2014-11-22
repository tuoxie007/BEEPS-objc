//
//  ViewController.m
//  BEEPS
//
//  Created by Jason Hsu on 14/11/22.
//  Copyright (c) 2014年 Jason Hsu. All rights reserved.
//

#import "ViewController.h"
#import "ImageHelper.h"

void copy_char_to_float(float *dst, unsigned char *src, size_t size)
{
    for (size_t i=0; i<size; i++) {
        dst[i] = src[i];
    }
}

void copy_float_to_char(unsigned char *dst, float *src, size_t size)
{
    for (size_t i=0; i<size; i++) {
        dst[i] = src[i];
    }
}

void P(float *rawData, int startIndex, int length, double photometricStandardDeviation, double spatialContraDecay) {
    double mu = 0.0;
    double rho = 1.0 + spatialContraDecay;
    double c = -0.5 / (photometricStandardDeviation * photometricStandardDeviation);
    for (int k=0; k<3; k++) {
        rawData[startIndex*4+k] /= rho;
        for (int i=startIndex+1; i<startIndex+length; i++) {
            mu = rawData[i*4+k] - rho * rawData[(i - 1)*4+k];
            mu = spatialContraDecay * exp(c * mu * mu);
            rawData[i*4+k] = rawData[(i - 1)*4+k] * mu + rawData[i*4+k] * (1.0 - mu) / rho;
        }
    }
}

void G(float *rawData, int startIndex, int length, double spatialContraDecay) {
    double mu = (1.0 - spatialContraDecay) / (1.0 + spatialContraDecay);
    for (int k=0; k<3; k++) {
        for (int i=startIndex; i<startIndex+length; i++) {
            rawData[i*4+k] *= mu;
        }
    }
}

void R(float *rawData, int startIndex, int length, double photometricStandardDeviation, double spatialContraDecay) {
    double mu = 0.0;
    double rho = 1.0 + spatialContraDecay;
    double c = -0.5 / (photometricStandardDeviation * photometricStandardDeviation);
    for (int k=0; k<3; k++) {
        rawData[(startIndex + length - 1) * 4 + k] /= rho;
        for (int i=startIndex + length - 2; i>=startIndex; i--) {
            mu = rawData[i*4+k] - rho * rawData[(i + 1)*4+k];
            mu = spatialContraDecay * exp(c * mu * mu);
            rawData[i*4+k] = rawData[(i + 1)*4+k] * mu + rawData[i*4+k] * (1.0 - mu) / rho;
        }
    }
}

void HV(unsigned char *rawData, int width, int height, double photometricStandardDeviation, double spatialDecay) {
    
    float *gData = malloc(width*height*4*sizeof(float));
    float *pData = malloc(width*height*4*sizeof(float));
    float *rData = malloc(width*height*4*sizeof(float));
    
    copy_char_to_float(gData, rawData, width*height*4);
    
    memcpy(pData, gData, width*height*4*sizeof(float));
    memcpy(rData, gData, width*height*4*sizeof(float));
    for (int k1=0; k1<height; k1++) {
        P(pData, k1*width, width, photometricStandardDeviation, 1.0 - spatialDecay);
        G(gData, k1*width, width, 1.0 - spatialDecay);
        R(rData, k1*width, width, photometricStandardDeviation, 1.0 - spatialDecay);
    }
    
    for (int k=0; k<width*height; k++) {
        for (int offset=0; offset<3; offset++) {
            rData[k*4+offset] += pData[k*4+offset] - gData[k*4+offset];
        }
    }
    
    for (int k=0; k<3; k++) {
        int m = 0;
        for (int k1 = 0; k1 < height; k1++) {
            int n = k1;
            for (int k2 = 0; k2 < width; k2++) {
                gData[n*4+k] = rData[m*4+k];
                m++;
                n += height;
            }
        }
    }
    
    memcpy(pData, gData, width*height*4*sizeof(float));
    memcpy(rData, gData, width*height*4*sizeof(float));
    for (int k2=0; k2<width; k2++) {
        P(pData, k2*height, height, photometricStandardDeviation, 1.0 - spatialDecay);
        G(gData, k2*height, height, 1.0 - spatialDecay);
        R(rData, k2*height, height, photometricStandardDeviation, 1.0 - spatialDecay);
    }
    
    for (int k=0; k<width*height; k++) {
        for (int offset=0; offset<3; offset++) {
            rData[k*4+offset] += pData[k*4+offset] - gData[k*4+offset];
        }
    }
    
    for (int k=0; k<3; k++) {
        int m = 0;
        for (int k1 = 0; k1 < width; k1++) {
            int n = k1;
            for (int k2 = 0; k2 < height; k2++) {
                rawData[n*4+k] = rData[m*4+k];
                m++;
                n += width;
            }
        }
    }
}

void VH(unsigned char *rawData, int width, int height, double photometricStandardDeviation, double spatialDecay) {
    
    float *gData = malloc(width*height*4*sizeof(float));
    float *pData = malloc(width*height*4*sizeof(float));
    float *rData = malloc(width*height*4*sizeof(float));
    
    for (int k=0; k<3; k++) {
        int m = 0;
        for (int k2=0; k2<height; k2++) {
            int n = k2;
            for (int k1=0; k1<width; k1++) {
                gData[n*4+k] = rawData[m*4+k];
                m++;
                n += height;
            }
        }
    }
    
    memcpy(pData, gData, width*height*4*sizeof(float));
    memcpy(rData, gData, width*height*4*sizeof(float));
    for (int k1=0; k1<width; k1++) {
        P(pData, k1*height, height, photometricStandardDeviation, 1.0 - spatialDecay);
        G(gData, k1*height, height, 1.0 - spatialDecay);
        R(rData, k1*height, height, photometricStandardDeviation, 1.0 - spatialDecay);
    }
    
    for (int k=0; k<width*height; k++) {
        for (int offset=0; offset<3; offset++) {
            rData[k*4+offset] += pData[k*4+offset] - gData[k*4+offset];
        }
    }
    
    for (int k=0; k<3; k++) {
        int m = 0;
        for (int k1 = 0; k1 < width; k1++) {
            int n = k1;
            for (int k2 = 0; k2 < height; k2++) {
                gData[n*4+k] = rData[m*4+k];
                m++;
                n += width;
            }
        }
    }
    
    memcpy(pData, gData, width*height*4*sizeof(float));
    memcpy(rData, gData, width*height*4*sizeof(float));
    for (int k2=0; k2<height; k2++) {
        P(pData, k2*width, width, photometricStandardDeviation, 1.0 - spatialDecay);
        G(gData, k2*width, width, 1.0 - spatialDecay);
        R(rData, k2*width, width, photometricStandardDeviation, 1.0 - spatialDecay);
    }
    
    for (int k=0; k<width*height; k++) {
        for (int offset=0; offset<1; offset++) {
            int i = k*4+offset;
            rawData[i] = pData[i]- gData[i] + rData[i];
        }
    }
}

@interface ViewController ()

@end

@implementation ViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    
    self.view.backgroundColor = [UIColor blackColor];
    
    UIImage *image = [UIImage imageNamed:@"mopi4.jpg"];
    UIImageView *imageView = [[UIImageView alloc] initWithImage:image];
    [self.view addSubview:imageView];
    imageView.contentMode = UIViewContentModeScaleAspectFit;
    imageView.frame = CGRectMake(0, 10, self.view.bounds.size.width, self.view.bounds.size.height*0.5-15);
    
    unsigned char *rawData = [self copyRGBAsFromImage:image];
    int width = (int)CGImageGetWidth(image.CGImage);
    int height = (int)CGImageGetHeight(image.CGImage);
    int pixCount = width * height;
    
    unsigned char *hvData = malloc(pixCount * 4);
    memcpy(hvData, rawData, pixCount * 4);
    
    unsigned char *vhData = malloc(pixCount * 4);
    memcpy(vhData, rawData, pixCount * 4);
    
    HV(hvData, width, height, 30, 0.1);
    VH(vhData, width, height, 30, 0.1);
    
    for (int i=0; i<pixCount*4; i++) {
        rawData[i] = (hvData[i] + vhData[i]) / 2;
    }
    
    image = [self createImageFromRGBAs:rawData width:width height:height];
    imageView = [[UIImageView alloc] initWithImage:image];
    [self.view addSubview:imageView];
    imageView.contentMode = UIViewContentModeScaleAspectFit;
    imageView.frame = CGRectMake(0, self.view.bounds.size.height/2+15, self.view.bounds.size.width, self.view.bounds.size.height*0.5-15);
}

- (unsigned char *)copyRGBAsFromImage:(UIImage*)image
{
    CGImageRef imageRef = [image CGImage];
    int width = (int)CGImageGetWidth(imageRef);
    int height = (int)CGImageGetHeight(imageRef);
    CGColorSpaceRef colorSpace = CGColorSpaceCreateDeviceRGB();
    unsigned char *rawData = (unsigned char*) calloc(height * width * 4, sizeof(unsigned char));
    NSUInteger bytesPerPixel = 4;
    NSUInteger bytesPerRow = bytesPerPixel * width;
    NSUInteger bitsPerComponent = 8;
    CGContextRef context = CGBitmapContextCreate(rawData, width, height,
                                                 bitsPerComponent, bytesPerRow, colorSpace,
                                                 kCGImageAlphaPremultipliedLast | kCGBitmapByteOrder32Big);
    NSAssert(context != nil, @"context not created");
    CGColorSpaceRelease(colorSpace);
    
    CGContextDrawImage(context, CGRectMake(0, 0, width, height), imageRef);
    CGContextRelease(context);
    
    return rawData;
}

- (UIImage *)createImageFromRGBAs:(unsigned char *)rawData width:(int)width height:(int)height
{
    size_t bufferLength = width * height * 4;
    CGDataProviderRef provider = CGDataProviderCreateWithData(NULL, rawData, bufferLength, NULL);
    size_t bitsPerComponent = 8;
    size_t bitsPerPixel = 32;
    size_t bytesPerRow = 4 * width;
    CGColorSpaceRef colorSpaceRef = CGColorSpaceCreateDeviceRGB();
    CGBitmapInfo bitmapInfo = kCGBitmapByteOrderDefault | kCGImageAlphaPremultipliedLast;
    CGColorRenderingIntent renderingIntent = kCGRenderingIntentDefault;
    
    CGImageRef iref = CGImageCreate(width,
                                    height,
                                    bitsPerComponent,
                                    bitsPerPixel,
                                    bytesPerRow,
                                    colorSpaceRef,
                                    bitmapInfo,
                                    provider,   // data provider
                                    NULL,       // decode
                                    YES,        // should interpolate
                                    renderingIntent);
    
    return [[UIImage alloc] initWithCGImage:iref];
}

- (BOOL)prefersStatusBarHidden
{
    return YES;
}

@end
