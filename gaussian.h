#include <vector>

typedef std::vector<double> Array;
typedef std::vector<Array> Matrix;
typedef std::vector<Matrix> Image;

inline Matrix get_gaussian(int height, int width, double sigma)
{
    Matrix kernel(height, Array(width));
    double sum=0.0;
    int i,j;

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] = exp(-(i*i+j*j)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
            sum += kernel[i][j];
        }
    }

    for (i=0 ; i<height ; i++) {
        for (j=0 ; j<width ; j++) {
            kernel[i][j] /= sum;
        }
    }

    return kernel;
}

inline std::vector<double> apply_filter(const std::vector<double>& image, const int width, const int height, const Matrix &filter){
    int filterHeight   = filter.size();
    int filterWidth    = filter[0].size();
    int newImageHeight = height-filterHeight+1;
    int newImageWidth  = width-filterWidth+1;
    int d,i,j,h,w;

    std::vector<double> newImage;
    newImage.resize(width * height);

    for (i=0 ; i<newImageHeight ; i++) {
        for (j=0 ; j<newImageWidth ; j++) {
            for (h=i ; h<i+filterHeight ; h++) {
                for (w=j ; w<j+filterWidth ; w++) {
                    newImage[i * width + j] += filter[h-i][w-j]*image[h * width + w];
                }
            }
        }
    }

    return newImage;
}