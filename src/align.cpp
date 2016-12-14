#include "align.h"
#include <string>
#include <cmath>
#include <limits.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <tuple> 
#define MAX_SHIFT 15

using std::string;
using std::cout;
using std::endl;
using std::vector;
using std::tie;
Image mirror(Image, int);
Image dis_mirror(Image, int);
double check_pi(double val) {
    if(val > 255) 
        return 255;
    if(val < 0)
        return 0;
    return val;
}
int check(int x, int n) 
{
    return (x >= 0) && (x < n);
}
Image make(Image b, Image g, Image r, int shift_i_bg, int shift_j_bg, int shift_i_gr, int shift_j_gr) 
{
    int n = b.n_rows;
    int m = b.n_cols;
    Image res(n, m);
    int left_i = n;
    int left_j = m;
    int right_i = 0;
    int right_j = 0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            int b_pos_row = i;
            int b_pos_col = j;
            int g_pos_row = b_pos_row + shift_i_bg;
            int g_pos_col = b_pos_col + shift_j_bg;
            int r_pos_row = g_pos_row + shift_i_gr;
            int r_pos_col = g_pos_col + shift_j_gr; 
            if(check(r_pos_row, n) && check(g_pos_row, n)) {
                if(check(r_pos_col, m) && check(g_pos_col, m)) {
                    left_i = std::min(left_i, i);
                    left_j = std::min(left_j, j);
                    right_i = std::max(right_i, i);
                    right_j = std::max(right_j, j);
                    res(i, j) = std::make_tuple(std::get<0>(r(r_pos_row, r_pos_col)), std::get<1>(g(g_pos_row, g_pos_col)), std::get<2>(b(b_pos_row, b_pos_col)));
                }
            } 
        }
    }
    return res.submatrix(left_i, left_j, right_i - left_i + 1, right_j - left_j + 1);
}
// it makes shift
Image make_shift(Image r, Image g, Image b, int x, int y, int n, int m, int pos_r, int pos_g, int pos_b, int shift_r, int shift_g, int shift_b) 
{   
    Image res(n, m);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            int r_pos_row = i + shift_r * x;
            int r_pos_col = j + shift_r * y;

            int g_pos_row = i + shift_g * x;
            int g_pos_col = j + shift_g * y;
            
            int b_pos_row = i + shift_b * x;
            int b_pos_col = j + shift_b * y;
            
            if(check(r_pos_row, n) && check(g_pos_row, n) && check(b_pos_row, n)) {
                if(check(r_pos_col, m) && check(g_pos_col, m) && check(b_pos_col, m)) {
                    int val_r = pos_r * (std::get<0>(r(r_pos_row, r_pos_col)));
                    int val_g = pos_g * (std::get<1>(g(g_pos_row, g_pos_col))); 
                    int val_b = pos_b * (std::get<2>(b(b_pos_row, b_pos_col)));
                    res(i, j) = std::make_tuple(val_r, val_g, val_b);
                }
            } else {
                res(i, j) = std::make_tuple(0, 0, 0);
            }           
        }
    }
    return res;
}
long double mse(Image src, double pos_r, double pos_g, double pos_b, int shift_i, int shift_j) 
{
    long double tmp = 0;
    int n = src.n_rows;
    int m = src.n_cols;
    int N = 0.05 * n;
    int M = 0.05 * m;
    for(int i = 0 + N; i < n - N; i++) {
        for(int j = 0 + M; j < m - M; j++) {
            double preprocess = (pos_r * (std::get<0>(src(i, j))) + pos_g * (std::get<1>(src(i, j)))+ pos_b * (std::get<2>(src(i, j))));
            tmp += preprocess * preprocess;
        }
    }
    tmp /= ((src.n_rows - std::abs(shift_i) - 2 * N) * (src.n_cols - std::abs(shift_j) - 2 * N ));
    return tmp;
}
long double corr(Image src, int pos_r, int pos_g, int pos_b, int shift_i, int shift_j) 
{
    long double tmp = 0;
    int n = src.n_rows/3;
    int m = src.n_cols;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            long double preprocess = ((std::get<0>(src(i, j))) * (std::get<1>(src(i, j))) * (std::get<2>(src(i, j))));
            tmp += preprocess;
        }
    }
    return tmp;
}
Image sub_pixel_bigger(Image src, int k) {
    int n = src.n_rows;
    int m = src.n_cols;
    Image res(k * n, m * k);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            for(int t = i * k; t < i * k + k; t++) {
                for(int p = j * k; p < j * k + k; p++) {
                    res(t, p) = src(i, j);
                }
            }
        }
    }
    return res;
}
Image sub_pixel_smaller(Image src, int k) {
    int n = src.n_rows/k;
    int m = src.n_cols/k;
    Image res(n, m);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res(i, j) = src(i * k, j * k);
        }
    }
    return res;
}
// r - 0 g - 1 b - 2  
Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    if(isSubpixel) {
        srcImage = sub_pixel_bigger(srcImage, subScale);
    }
    
    Image b = srcImage.submatrix(0, 0, srcImage.n_rows/3, srcImage.n_cols);
    Image g = srcImage.submatrix(srcImage.n_rows / 3, 0, srcImage.n_rows / 3, srcImage.n_cols);
    Image r = srcImage.submatrix(2 * srcImage.n_rows / 3, 0, srcImage.n_rows / 3, srcImage.n_cols);
    // b & g
    long double ans_bg = 0;
    int shift_i_bg = 0;
    int shift_j_bg = 0;
    for(int i = -MAX_SHIFT; i <= MAX_SHIFT; i++){
        for(int j = -MAX_SHIFT; j <= MAX_SHIFT; j++){
            Image img = make_shift(r, g, b, i, j, srcImage.n_rows/3, srcImage.n_cols, 0, 1, 1, 0, 1, 0);
            long double tmp_mse = mse(img, 0, 1, -1, i, j);
            if(ans_bg > tmp_mse  || (j == -MAX_SHIFT && i == -MAX_SHIFT)) {
                ans_bg = tmp_mse;
                shift_i_bg = i;
                shift_j_bg = j;
                
            }
        }
    }
    // g && r
    long double ans_gr = 0;
    int shift_i_gr = 0;
    int shift_j_gr = 0;
    for(int i = -MAX_SHIFT; i <= MAX_SHIFT; i++){
        for(int j =-MAX_SHIFT; j <= MAX_SHIFT; j++){
            Image img = make_shift(r, g, b, i, j, srcImage.n_rows/3, srcImage.n_cols, 1, 1, 0, 1, 0, 0);
            long double tmp_mse = mse(img, 1, -1, 0, i, j);
            if(ans_gr > tmp_mse || (i == -MAX_SHIFT && j == -MAX_SHIFT)) {
                ans_gr = tmp_mse;
                shift_i_gr = i;
                shift_j_gr = j;
            }
        }
    }
    Image res = make(b, g, r, shift_i_bg, shift_j_bg, shift_i_gr, shift_j_gr);
    if(isSubpixel) {
        res = sub_pixel_smaller(res, subScale);
    }
    if(isPostprocessing) {
        if(postprocessingType =="--gray-world") {
            res = gray_world(res);
        } else if(postprocessingType == "--unsharp"){
            if(isMirror) {
                res = dis_mirror(unsharp(mirror(res, 1)), 1); 
            } else {
                res = unsharp(res);
            }
            
        } else if(postprocessingType == "--autocontrast") {
            res = autocontrast(res, fraction);
        }
    }

    return res;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    int n = src_image.n_rows;
    int m = src_image.n_cols;
    double prod[3];
    Image tmp(n, m);
    tmp = src_image.deep_copy();
    for(int i = 1; i < n - 1; i++) {
        for(int j = 1; j < m - 1; j++) {
                prod[0] = 0;
                prod[0] += (-1.0/6) * std::get<0>(src_image(i - 1, j - 1));
                prod[0] += (-2.0/3) * std::get<0>(src_image(i - 1, j));
                prod[0] += (-1.0/6) * std::get<0>(src_image(i - 1, j + 1));
                prod[0] += (-2.0/3) * std::get<0>(src_image(i, j - 1));
                prod[0] += (13.0/3) * std::get<0>(src_image(i, j));
                prod[0] += (-2.0/3) * std::get<0>(src_image(i, j + 1));
                prod[0] += (-1.0/6) * std::get<0>(src_image(i + 1, j - 1));
                prod[0] += (-2.0/3) * std::get<0>(src_image(i + 1, j));
                prod[0] += (-1.0/6) * std::get<0>(src_image(i + 1, j + 1));
                
                if(prod[0] < 0) prod[0] = 0;
                prod[1] = 0;
                prod[1] += (-1.0/6) * std::get<1>(src_image(i - 1, j - 1));
                prod[1] += (-2.0/3) * std::get<1>(src_image(i - 1, j));
                prod[1] += (-1.0/6) * std::get<1>(src_image(i - 1, j + 1));
                prod[1] += (-2.0/3) * std::get<1>(src_image(i, j - 1));
                prod[1] += (13.0/3) * std::get<1>(src_image(i, j));
                prod[1] += (-2.0/3) * std::get<1>(src_image(i, j + 1));
                prod[1] += (-1.0/6) * std::get<1>(src_image(i + 1, j - 1));
                prod[1] += (-2.0/3) * std::get<1>(src_image(i + 1, j));
                prod[1] += (-1.0/6) * std::get<1>(src_image(i + 1, j + 1));
                
                if(prod[1] < 0) prod[1] = 0;
                prod[2] = 0;
                prod[2] += (-1.0/6) * std::get<2>(src_image(i - 1, j - 1));
                prod[2] += (-2.0/3) * std::get<2>(src_image(i - 1, j));
                prod[2] += (-1.0/6) * std::get<2>(src_image(i - 1, j + 1));
                prod[2] += (-2.0/3) * std::get<2>(src_image(i, j - 1));
                prod[2] += (13.0/3) * std::get<2>(src_image(i, j));
                prod[2] += (-2.0/3) * std::get<2>(src_image(i, j + 1));
                prod[2] += (-1.0/6) * std::get<2>(src_image(i + 1, j - 1));
                prod[2] += (-2.0/3) * std::get<2>(src_image(i + 1, j));
                prod[2] += (-1.0/6) * std::get<2>(src_image(i + 1, j + 1));
                if(prod[2] < 0) prod[2] = 0;
               
            tmp(i, j) = std::make_tuple(check_pi(prod[0]), check_pi(prod[1]), check_pi(prod[2]));
        }
    }
    return tmp;
}

Image gray_world(Image src_image) {
    int n = src_image.n_rows;
    int m = src_image.n_cols;
    Image tmp(n, m);
    tmp = src_image.deep_copy();
    double s_r = 0;
    double s_g = 0;
    double s_b = 0;
    for(int i = 1; i < n - 1; i++) {
        for(int j = 1; j < m - 1; j++) {
        
                    s_r += std::get<0>(src_image(i, j));
                    s_g += std::get<1>(src_image(i, j));
                    s_b += std::get<2>(src_image(i, j));
        }

    }
    s_r /= n*m;
    s_g /= n*m;
    s_b /= n*m;
    double s = (s_r + s_g + s_b)/3;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            tmp(i, j) = std::make_tuple(check_pi(s * std::get<0>(src_image(i, j))/s_r),check_pi(s * std::get<1>(src_image(i, j))/s_g),check_pi( s * std::get<2>(src_image(i, j))/s_b));
        }
    }
    return tmp;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.
    return src_image;
}

Image autocontrast(Image src_image, double fraction) {
    int n = src_image.n_rows;
    int m = src_image.n_cols;
    vector <double> data;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            data.push_back(0.2125 * std::get<0>(src_image(i, j)) + 0.7154 * std::get<1>(src_image(i, j)) + 0.0721 * std::get<2>(src_image(i, j)));
        }
    }
    sort(data.begin(), data.end());
    int part = fraction * n * m;
    double y_min = data[part + 1];
    double y_max = data[n*m - part - 1];
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++) {
            double val_r = (int(std::get<0>(src_image(i, j))) - y_min) * 255.0 /(y_max - y_min); 
            double val_g = (int(std::get<1>(src_image(i, j))) - y_min) * 255.0 /(y_max - y_min);
            double val_b = (int(std::get<2>(src_image(i, j))) - y_min) * 255.0 /(y_max - y_min);
            src_image(i, j) = std::make_tuple(check_pi(val_r), check_pi(val_g), check_pi(val_b));
        }
    }
    return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    return src_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

int find(unsigned int *a, int k) {
    int sum = 0;
    int i = 0;
    for(i = 0; i < 256; i++) {
        sum += a[i];
        if(sum >= k) {
            break;
        }
    }
    return i;
}
Image mirror(Image src_image, int k) {
    int n = src_image.n_rows;
    int m = src_image.n_cols;
    Image res(2 * k + n, 2 * k + m);
    //base
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            res(i + k, j + k) = src_image(i, j);
        }
    }

    int new_i = k;
    int new_j = k - 1;
    for(int i = 0; i < n; i++) {
        new_j = k - 1;
        for(int j = 0; j < k; j++) {
            res(new_i, new_j) = src_image(i, j);
            new_j--;
        }
        new_i++;
    }
    new_i = k;
    for(int i = 0; i < n; i++) {
        new_j = k + m;
        for(int j = m - 1; j >= m - k; j--) {
            res(new_i, new_j) = src_image(i, j);
            new_j++;
        }
        new_i++;
    }
    new_i = k - 1;
    for(int i = 0; i < k; i++) {
        new_j = k;
        for(int j = 0; j < m; j++) {
            res(new_i, new_j) = src_image(i, j);
            new_j++;
        }
        new_i--;
    }
    new_i = n + k;
    for(int i = n - 1; i >= n - k; i--) {
        new_j = k;
        for(int j = 0; j < m; j++) {
            res(new_i, new_j) = src_image(i, j);
            new_j++;
        }
        new_i++;
    }
    new_i = k - 1;
    for(int i = 0; i < k; i++) {
        new_j = k - 1;
        for(int j = 0; j < k; j++){
            res(new_i, new_j) = src_image(i, j);
            new_j--;
        }
        new_i--;
    }
    new_i =  k - 1;
    for(int i = 0; i < k; i++) {
        new_j = k + m;
        for(int j = m - 1; j >= m - k; j--){
            res(new_i, new_j) = src_image(i, j);
            new_j++;
        }
        new_i--;
    }
    new_i = n + k;
    for(int i = n - 1; i >= n - k; i--) {
        new_j = k - 1;
        for(int j = 0; j < k; j++){
            res(new_i, new_j) = src_image(i, j);
            new_j--;
        }
        new_i++;
    }
    new_i = n + k;
    for(int i = n - 1; i >= n - k; i--) {
        new_j = m +  k;
        for(int j = 0; j < k; j++){
            res(new_i, new_j) = src_image(i, j);
            new_j++;
        }
        new_i++;
    }
    return res;
}
Image dis_mirror(Image src, int k) 
{
    return src.submatrix(k, k, src.n_rows - 2 * k, src.n_cols - 2 * k);
}
Image median(Image src_image, int radius) {
    int n = src_image.n_rows;
    int m = src_image.n_cols;
    Image tmp(n, m);
    tmp = src_image.deep_copy();
    unsigned int datar[256];
    unsigned int datag[256];
    unsigned int datab[256];
    for(int i = 0; i < 256; i++) {
        datar[i] = 0;
        datag[i] = 0;
        datab[i] = 0;
    }
    int val_r = 0, val_g = 0, val_b = 0;
    int ind = (2 * radius + 1) * (2 * radius + 1)/2 + 1;
    vector <unsigned int> data;
    for(int i = radius; i < n - radius; i++) {
        for(int j = radius; j < m - radius; j++) {
            for(int k = -radius + i; k <= i + radius; k++) {
                for(int t = -radius + j; t <= j + radius; t++) {
                    tie(val_r, val_g, val_b) = src_image(k, t);
                            datar[val_r]++;
                            datag[val_g]++;
                            datab[val_b]++;
                }
            }
            tmp(i, j) = std::make_tuple(check_pi(find(datar, ind)),check_pi(find(datag, ind)), check_pi(find(datab, ind)));
            for(int k = -radius + i; k <= i + radius; k++) {
                for(int t = -radius + j; t <= j + radius; t++) {
                    tie(val_r, val_g, val_b) = src_image(k, t);
                            datar[val_r]--;
                            datag[val_g]--;
                            datab[val_b]--;
                }
            }
        }
    }
    return tmp;
}

Image median_linear(Image src_image, int radius) {
    int n = src_image.n_rows;
    int m = src_image.n_cols;
    unsigned int datar[256];
    unsigned int datag[256];
    unsigned int datab[256];
    for(int i = 0; i < 256; i++) {
        datar[i] = 0;
        datag[i] = 0;
        datab[i] = 0;
    }
    Image tmp(n, m);
    tmp = src_image.deep_copy();
    int flag = 0;
    int val_r = 0, val_g = 0, val_b = 0;
    
    int ind = (2 * radius + 1) * (2 * radius + 1)/2 + 1;
    int i = radius;
    int flag2 = 1;
    while(i < n - radius) {
        for(int j = radius; j < m - radius; j++) {
            if(!flag) {
                for(int t = i - radius; t <= i + radius;t++) {
                    for(int k = j - radius; k <= j + radius; k++) {
                        tie(val_r, val_g, val_b) = src_image(t, k);
                        datar[val_r]++;
                        datag[val_g]++;
                        datab[val_b]++;
                    }
                }
                tmp(i, j) = std::make_tuple(check_pi(find(datar, ind)),check_pi(find(datag, ind)), check_pi(find(datab, ind)));
                flag = 1;
            } else {
                if (flag2) {
                
                    for(int t = i - radius; t <= i + radius; t++) {
                        tie(val_r, val_g, val_b) = src_image(t, j + radius);
                            datar[val_r]++;
                            datag[val_g]++;
                            datab[val_b]++;
                    }
                    for(int t = i - radius; t <= i + radius; t++) {
                        tie(val_r, val_g, val_b) = src_image(t, j - radius - 1);
                            datar[val_r]--;
                            datag[val_g]--;
                            datab[val_b]--;
                    }
                    tmp(i, j) = std::make_tuple(check_pi(find(datar, ind)),check_pi(find(datag, ind)), check_pi(find(datab, ind)));
                    if(j == m - radius - 1) {
                        if(i != n - radius - 1) {
                            for(int k = j - radius; k <= j + radius; k++) {
                                tie(val_r, val_g, val_b) = src_image(i - radius, k);
                                datar[val_r]--;
                                datag[val_g]--;
                                datab[val_b]--;
                            }
                            for(int k = j - radius; k <= j + radius; k++) {
                                tie(val_r, val_g, val_b) = src_image(i + radius + 1, k);
                                datar[val_r]++;
                                datag[val_g]++;
                                datab[val_b]++;
                            }
                            tmp(i+1, j) = std::make_tuple(check_pi(find(datar, ind)),check_pi(find(datag, ind)), check_pi(find(datab, ind)));
                            i++;
                        } else i++;
                    }
                }
                flag2 = 1;
            }
        }
            if(i < n - radius) {
            for(int j = m - radius - 2; j >= radius; j--) {
                for(int t = i - radius; t <= i + radius; t++) {
                    tie(val_r, val_g, val_b) = src_image(t, j - radius);
                        datar[val_r]++;
                        datag[val_g]++;
                        datab[val_b]++;
                }
                for(int t = i - radius; t <= i + radius; t++) {
                    tie(val_r, val_g, val_b) = src_image(t, j + radius + 1);
                        datar[val_r]--;
                        datag[val_g]--;
                        datab[val_b]--;
                }
                tmp(i, j) = std::make_tuple(check_pi(find(datar, ind)),check_pi(find(datag, ind)), check_pi(find(datab, ind)));
                if(j == radius) {
                    if(i != n - radius - 1) {
                        for(int k = j - radius; k <= j + radius; k++) {
                            tie(val_r, val_g, val_b) = src_image(i - radius, k);
                            datar[val_r]--;
                            datag[val_g]--;
                            datab[val_b]--;
                        }
                        for(int k = j - radius; k <= j + radius; k++) {
                            tie(val_r, val_g, val_b) = src_image(i + radius + 1, k);
                            datar[val_r]++;
                            datag[val_g]++;
                            datab[val_b]++;
                        }
                        tmp(i+1, j) = std::make_tuple(check_pi(find(datar, ind)),check_pi(find(datag, ind)), check_pi(find(datab, ind)));
                        i++;

                    } else i++;
                }
                flag2 = 0;
            }
        }
    }        
    return tmp;
}
Image median_const(Image src_image, int radius) {
    return src_image;
}
Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}
