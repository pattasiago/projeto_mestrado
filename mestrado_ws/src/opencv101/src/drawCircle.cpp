#include <ros/ros.h>
/*****************************************************************************
** Include
*****************************************************************************/
#include <ros/ros.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Header.h>
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
#include <opencv101/optimizerResult.h>
#include <opencv101/imageParams.h>
#include <opencv101/Lane.h>
#include <sensor_msgs/image_encodings.h>
#include <image_transport/image_transport.h>
#include <sensor_msgs/Image.h>
//namespace enc = sensor_msgs::image_encodings;
#define pi 3.141592

using namespace cv;
using namespace std;
using namespace Eigen;

ros::Publisher imageParams;

image_transport::Publisher image_pub;
cv_bridge::CvImage img_bridge;
sensor_msgs::Image img_msg;

void setImageParams(const cv::Mat &I){
		
	//----------------------IMAGE VARIABLES-------------------------------------------------------------------------------
	int color = 1;
	std::vector<float> curve;
	int laneHeight;
	int qtd_coordinates = 8;
	std::vector<cv::Point> mainCoordinates;
	std::vector<float> mainCurvatures;
	float ori;
	opencv101::imageParams params;
	
	getLaneParameters(I,color,curve,laneHeight,qtd_coordinates,ori,mainCoordinates,mainCurvatures);
	
	cv::flip(I,I,0);
	polylines(I, mainCoordinates, false, Scalar(0, 0, 255), 2);
	for(int k = 0; k < qtd_coordinates; k++){
		circle(I, Point2f(mainCoordinates[k].x,mainCoordinates[k].y), 2, (0,0,255),2);
	}
	cv::flip(I,I,0);
	cv::imshow("Image", I);
	cv::waitKey(3);
	
	int i_coord = int(0.5f+mainCoordinates.size()/2)-1; //central coordinate

	std_msgs::Header header = std_msgs::Header();
	header.stamp = ros::Time::now();
	
	img_bridge = cv_bridge::CvImage(header, sensor_msgs::image_encodings::BGR8, I);
			img_bridge.toImageMsg(img_msg);
	image_pub.publish(img_msg);
	
	params.header = header;
	
	params.circParams[0] = mainCoordinates[i_coord].x;//float32[3] circParams
	params.circParams[1] = mainCoordinates[i_coord].y;//I.rows - mainCoordinates[i_coord].y;//float32[3] circParams
	params.circParams[2] = mainCurvatures[i_coord];//
	
	params.poseori[0] = mainCoordinates[i_coord].x;
	params.poseori[1] = mainCoordinates[i_coord].y;//I.rows - mainCoordinates[i_coord].y; //float32
	params.poseori[2] = -(atan2(mainCoordinates[i_coord].x,mainCoordinates[i_coord].y) * 180)/pi;

	params.ponto2[0] = int(mainCoordinates[i_coord-1].x); //int32
	params.ponto2[1] = int(mainCoordinates[i_coord-1].y);//I.rows - int(mainCoordinates[i_coord-1].y);

	params.ponto3[0] = int(mainCoordinates[i_coord+1].x);
	params.ponto3[1] = int(mainCoordinates[i_coord+1].y);//I.rows -int(mainCoordinates[i_coord-1].y);

	imageParams.publish(params);
}


void imageCallback(const sensor_msgs::ImageConstPtr& msg){
	cv::Mat image;
    try
      {
        image = cv_bridge::toCvCopy(msg, "bgr8")->image;
      }
    catch (cv_bridge::Exception& e)
      {
        ROS_ERROR("cv_bridge exception: %s", e.what());
        return;
      }
    cv::resize(image,image,cv::Size(640,480), 0, 0, CV_INTER_LINEAR);
    setImageParams(image);
}


int main(int argc, char **argv)
{
    
	//------------------------------------------ROS---------------------------------
    ros::init(argc, argv, "image_converter");
    ros::NodeHandle n;
    //ros::Subscriber cam = n.subscribe("/usb_cam/image_raw", 1, imageCallback);

    image_transport::ImageTransport *it;
    it = new image_transport::ImageTransport(n);
    image_pub = it->advertise("/image_raw/result",1);

    imageParams = n.advertise<opencv101::imageParams>("/camera_params", 10);	
    ros::Subscriber cam = n.subscribe("/usb_cam/image_raw", 1, imageCallback);	///camera/rgb/image_raw

    ros::Rate loop_rate(10);  
    //while(ros::ok()) {
	//ros::Time time = ros::Time::now();
    //    ros::spinOnce();
    //    loop_rate.sleep();
    //}
    ros::spin();
    return 0;
}
