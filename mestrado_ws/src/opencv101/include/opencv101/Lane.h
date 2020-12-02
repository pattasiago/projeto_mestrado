#include <ros/ros.h>
/*****************************************************************************
** Include
*****************************************************************************/
#include <ros/ros.h>
#include <std_msgs/Float32.h>
#include <ros/ros.h>
#include <geometry_msgs/PoseArray.h>
#include <geometry_msgs/Pose2D.h>
#include <geometry_msgs/PoseWithCovariance.h>
#include <geometry_msgs/Twist.h>
#include <nav_msgs/Odometry.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <sstream>
#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/QR>
#include <opencv101/optimizerParams.h>
#include "opencv101/optimizerResult.h"
//namespace enc = sensor_msgs::image_encodings;
#define pi 3.141592

using namespace cv;
using namespace std;
using namespace Eigen;

void curveFitRansac(const std::vector<cv::Point> &points, std::vector<float> &coeffy, int order, int rows, int cols);
void curveFitMMQ(const std::vector<cv::Point> &points, std::vector<float> &coeffy, int order, int rows, int cols);
void getLaneline(const cv::Mat &imgBGR, cv::Mat &laneLineBW, const std::vector<cv::Scalar> &laneThreshHSV);
void getLaneLineCurve(const cv::Mat &laneLineBW,std::vector<float> &curve,int order);
float getLaneLineCurvatureRadius_pxl(const std::vector<float> &curve,float y_coordinate);
float getOrientation(const cv::Point &cameraCenter,const cv::Point &mainCoordinate);
int getLaneHeight_pxl(const std::vector<float> &curve, int rows, int cols);
void getMainCoordinates(const std::vector<float> &curve,const int height,int qtd_points, std::vector<cv::Point> &mainCoordinates);
void getMainCurvaturesRadius(const std::vector<float> &curve, const std::vector<cv::Point> &mainCoordinates,std::vector<float> &mainCurvatures);
void getColorHSVThreshRange(std::vector<cv::Scalar> &range,int i);
void getLaneParameters(const cv::Mat &I,int color,std::vector<float> &curve,int laneHeight,int qtd_coordinates,float ori,std::vector<cv::Point> &mainCoordinates,std::vector<float> &mainCurvatures);


