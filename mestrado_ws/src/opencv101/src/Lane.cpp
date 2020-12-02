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
#include <opencv101/Lane.h>
#include <opencv2/ximgproc/ridgefilter.hpp>
//namespace enc = sensor_msgs::image_encodings;
#define pi 3.141592

using namespace cv;
using namespace std;
using namespace Eigen;
using namespace cv::ximgproc;

void curveFitRansac(const std::vector<cv::Point> &points, std::vector<float> &coeffy, int order, int rows, int cols)
{
    
    order = 2;
    coeffy.resize(order+1);


    if(points.size() < 3){
        for (int k = 0; k < order+1; k++){
            coeffy[k] = 0;   
         }
        return;
    }

    MatrixXf A(points.size(), order+1);
    VectorXf B(points.size());
    VectorXf best_model(3);
    VectorXf X; 
    VectorXf X2;    

    int n_data = points.size();
    int divisor = n_data;
    int N = 55;//65;	//iterations 
    double T = 1*36;   // residual threshold

    int n_sample = 3;
    int max_cnt = 0;

    int nThreads = std::max(1, omp_get_max_threads());
    omp_set_dynamic(0); // Explicitly disable dynamic teams
    omp_set_num_threads(nThreads);


    int64_t start = cv::getTickCount();

    //----------------------- Construindo Matriz de Pontos---------------------//
    for (int i = 0; i < points.size(); i++){
        B(i) = points[i].x;
        A(i,0) = 1;        
        A(i,1) = rows - points[i].y ;         
        A(i,2) = std::pow(rows - points[i].y, 2);                 
    }

    #pragma omp parallel for private(X) shared(max_cnt,best_model)
    for( int i=0 ; i<N ; i++ )
    {
        //random sampling - 3 point  
        int k[3] = {-1, } ;
        k[0] = floor((rand()%divisor+1))+1;

        do
        {
            k[1] = floor((rand()%divisor+1))+1;
        }while(k[1]==k[0] || k[1]<0) ;

        do
        {
            k[2] = floor((rand()%divisor+1))+1;
        }while(k[2]==k[0] || k[2]==k[1] || k[2]<0) ;


        //model estimation
        MatrixXf AA(3,3);
        VectorXf BB(3);    
        
        for( int j=0 ; j<3 ; j++ ){
            BB(j) = points[k[j]].x;
            AA(j,0) = 1;        
            AA(j,1) = rows - points[k[j]].y ;                     
            AA(j,2) = std::pow(rows - points[k[j]].y, 2);         
        }

        X.resize(order+1);        
        X = AA.colPivHouseholderQr().solve(BB);

        //evaluation 
        VectorXf residual;    
        residual = A*X - B;
        residual = residual.cwiseAbs();
        int cnt = 0;

        for( int j=0 ; j<n_data ; j++ ){

            double data = residual(j) ;
            if( data < T ){
                cnt++ ;
            }
        }

        #pragma omp critical
        if( cnt > max_cnt ) {
            best_model = X ;
            max_cnt = cnt ;
        }
    }

  //-------------------------------------------------------------------LS fitting 
    VectorXf residual;    	
    residual = A*best_model - B;
    residual = residual.cwiseAbs();
    std::vector<int> vec_index ;

	for( int k=0 ; k<n_data ; k++ ){

		double data = residual(k) ;
		if( data < T ) {
	 		vec_index.push_back(k) ;
	 	}
	}

    MatrixXf A2(vec_index.size(),3);    	
    VectorXf B2(vec_index.size());

    //#pragma omp parallel for
    for( int j=0 ; j<vec_index.size() ; j++ ){
        B2(j) = points[vec_index[j]].x;

        A2(j,0) = 1;       
        A2(j,1) = rows - points[vec_index[j]].y ;                 
        A2(j,2) = std::pow(rows - points[vec_index[j]].y, 2);         
    }
     
     X2.resize(order+1);
     X2 = A2.colPivHouseholderQr().solve(B2);
    
    best_model = X2;
    
    //std::cout << "Best Model2: " << best_model(0) << " " << best_model(1) <<" "<< best_model(2)  << endl;
    int64_t end = cv::getTickCount();    
    //std::cout << "RANSAC took: " << ((end - start) / cv::getTickFrequency()) * 1000.0 << " ms." << std::endl;	

    for (int k = 0; k < order+1; k++){
       coeffy[k] =  best_model(k);   
    }
}

void curveFitMMQ(const std::vector<cv::Point> &points, std::vector<float> &coeffy, int order, int rows, int cols)
{
    order = 2;
    coeffy.resize(order+1);

	if(points.size() < 3){
        for (int k = 0; k < order+1; k++){
            coeffy[k] = 0;   
         }
        return;
    }

    MatrixXf A(points.size(),order+1);
    VectorXf B(points.size());
    VectorXf coefs;
    int64_t start = cv::getTickCount();
    for (int i = 0; i < points.size(); i++){
        B(i) = points[i].x;
        for (int j = 0; j < order+1; j++){
            A(i,j) = std::pow(rows - points[i].y, j);         
        }    
    }

    coefs = A.colPivHouseholderQr().solve(B);
    int64_t end = cv::getTickCount();
    //std::cout << "MMQ took: " << ((end - start) / cv::getTickFrequency()) * 1000.0 << " ms." << std::endl;	
    coefs.resize(order+1);

    for (int k = 0; k < order+1; k++){
        coeffy[k] =  coefs(k);    
    }
}

void getLaneline(const cv::Mat &imgBGR, cv::Mat &laneLineBW, const std::vector<cv::Scalar> &laneThreshHSV){

    cv::Mat imgHSV, maskHSV;
    cv::Scalar minHSV = laneThreshHSV[0]; 
    cv::Scalar maxHSV = laneThreshHSV[1];  
    cv::cvtColor(imgBGR, imgHSV, CV_BGR2HSV);
    cv::inRange(imgHSV, minHSV, maxHSV, maskHSV);
    // cv::imshow("SegmentedFinal",maskHSV);
    // cv::waitKey(3);
    laneLineBW = maskHSV;
};

void getLaneLineCurve(const cv::Mat &laneLineBW,std::vector<float> &curve,int order){
    std::vector<cv::Point> locations; 
    cv:: Mat img;  
    int row_length = laneLineBW.rows;
    int col_length = laneLineBW.cols;
    //cv::flip(laneLineBW,img,0);
    cv::findNonZero(laneLineBW, locations); // returns non-zero pixels position
    curveFitMMQ(locations,curve,order,row_length,col_length); // return quadratic function
}

float getLaneLineCurvatureRadius_pxl(const std::vector<float> &curve,float y_coordinate){
    float curvature = 0;
    float dx_dy = (2*curve[2]*y_coordinate) + curve[1];
    float d2x_dy2 = 2*curve[2];
    curvature = sqrt(pow((1+pow(dx_dy,2)),3))/(d2x_dy2+0.000001); // to avoid division by 0
    curvature = 1/curvature;
    //return curvature;
    return abs(curvature);
}

void getOrientation(const cv::Point &cameraCenter,const cv::Point &mainCoordinate, float orientacao){
    float ori = 0;
    ori = (atan2(mainCoordinate.x - cameraCenter.x,mainCoordinate.y) * 180)/pi ;
    orientacao = ori;
}

int getLaneHeight_pxl(const std::vector<float> &curve, int rows, int cols){
    float ans = 0;
    for(int i = rows; i > 0; i--){
        ans = curve[2]*i*i + curve[1]*i + curve[0];
        if(int(0.5f + ans) > 0 && int(0.5f + ans) < cols){ 
            return i;
        }   
    }
    return 0;
}

void getMainCoordinates(const std::vector<float> &curve,const int height,int qtd_points, std::vector<cv::Point> &mainCoordinates){
    float x_coordinate;
    mainCoordinates.resize(qtd_points);
    float height_space = height/(qtd_points-1);
    float y_coordinate = height;
    for(int i=0; i < qtd_points; i++){
        x_coordinate = curve[2]*y_coordinate*y_coordinate + curve[1]*y_coordinate + curve[0];
        //std::cout << "Main Coord: (" << x_coordinate << "," << y_coordinate<< ")"<< std::endl;
        mainCoordinates[i] = Point2f(x_coordinate,y_coordinate);
        y_coordinate = y_coordinate - height_space;
    } 
}

void getMainCurvaturesRadius(const std::vector<float> &curve, const std::vector<cv::Point> &mainCoordinates,std::vector<float> &mainCurvatures){
    mainCurvatures.resize(mainCoordinates.size());
    float curvature;
    float y_coordinate;
    for(int i=0; i < mainCoordinates.size(); i++){
        y_coordinate = mainCoordinates[i].y;
        curvature = getLaneLineCurvatureRadius_pxl(curve,y_coordinate);
        //std::cout << "Main Coord: " << curvature << std::endl;
        mainCurvatures[i] = curvature;
    } 
}

void getColorHSVThreshRange(std::vector<cv::Scalar> &range,int i){
    range.resize(2);
    switch (i){
            case 1://yellow line
                range[0] = cv::Scalar(20,90,0); //minHsv 20 60 20
                range[1] = cv::Scalar(40,255,255); //maxHsv 40 250 255
		        //range[0] = cv::Scalar(0, 61, 106);
		        //range[1] = cv::Scalar(45, 255, 194); // Calibration First Approach
                break;

            default:
                range[0] = cv::Scalar(20,90,0); //minHsv
                range[1] = cv::Scalar(40,255,255); //maxHsv
                break;
    }
}


void getLaneParameters(const cv::Mat &I,int color,std::vector<float> &curve,int laneHeight,int qtd_coordinates,float ori,std::vector<cv::Point> &mainCoordinates,std::vector<float> &mainCurvatures){
    std::vector<cv::Scalar> rangeColor;
    cv::Mat laneLineBW;
    std::vector<float> curva;
    std::vector<cv::Point> mainCord;
    std::vector<float> mainCurv;
    float orientation;
    int lanesize;

    getColorHSVThreshRange(rangeColor,color);//yellow range
    getLaneline(I,laneLineBW,rangeColor);
    getLaneLineCurve(laneLineBW,curva,2);
    lanesize = getLaneHeight_pxl(curva,laneLineBW.rows,laneLineBW.cols);
    getMainCoordinates(curva,lanesize,qtd_coordinates,mainCord);
    getOrientation(Point2f(I.cols/2,1),mainCord[int(0.5f+mainCord.size()/2)-1],orientation);
    getMainCurvaturesRadius(curva,mainCord,mainCurv);
    
    curve = curva;
    laneHeight = lanesize;
    ori = orientation;
    mainCoordinates = mainCord;
    mainCurvatures = mainCurv;

}


// int main(int argc, char **argv)
// {
//     cv::Mat I, yuv, channels[3], merged;
//     int color = 1;
//     std::vector<float> curve;
//     int laneHeight;
//     int qtd_coordinates = 8;
//     std::vector<cv::Point> mainCoordinates;
//     std::vector<float> mainCurvatures;
//     float ori;

//     I = cv::imread(argv[1], CV_LOAD_IMAGE_COLOR);
//     cv::resize(I,I,cv::Size(int(640/3),int(480/3)), 0, 0, CV_INTER_LINEAR);
    
//     getLaneParameters(I,color,curve,laneHeight,qtd_coordinates,ori,mainCoordinates,mainCurvatures);
    
//     imshow("Original Image", I);
    
//     cv::flip(I,I,0);
//     polylines(I, mainCoordinates, false, Scalar(0, 0, 255), 2);
//     for(int k = 0; k < qtd_coordinates; k++){
//         circle(I, Point2f(mainCoordinates[k].x,mainCoordinates[k].y), 2, (0,0,255),2);
//     }
//     cv::flip(I,I,0);
    
// 	int i_coord = int(0.5f+mainCoordinates.size()/2)-1; //central coordinate
// 	std::cout << "Z: " << abs(int(640/6) - mainCoordinates[i_coord].x) << endl;
// 	std::cout << "Curvature: " << mainCurvatures[i_coord] << endl;
// 	std::cout << "Orien: " << -(atan2(mainCoordinates[i_coord].x,mainCoordinates[i_coord].y) * 180)/pi << endl;



//     imshow("Image Before", I);
//     waitKey(0);

//     ros::init(argc, argv, "LaneLine");
//     ros::NodeHandle n;
//     ros::Rate loop_rate(10);
      
//     while(ros::ok()) {
//     ros::Time time = ros::Time::now();
//         ros::spinOnce();
//         loop_rate.sleep();
//     }
//     ros::spin();
//     return 0;
// }
