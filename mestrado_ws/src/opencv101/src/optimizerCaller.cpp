#include "ros/ros.h"
#include "std_msgs/String.h"
#include <opencv101/optimizerParams.h>
#include <iostream>
#include "opencv101/optimizerResult.h"
#include <sstream>


extern "C" 
{
	#include "stdio.h"
	#include "stdlib.h"
	#include "string.h"
	#include "optimizer/matrix.h"
	void optimizer(double Params[], double x00[], double u00[], double lij00[], double psiij00[], double betaij00[], double gammaij00[], double result[]);
}

using namespace std;

ros::Publisher pub;
void chatterCallback(const opencv101::optimizerParams& msg)
{

	optimizer(const_cast<double*>(&msg.Parameters[0]), const_cast<double*>(&msg.x[0]), const_cast<double*>(&msg.u[0]), const_cast<double*>(&msg.LFFClij[0]), const_cast<double*>(&msg.LFFCpsiij[0]), const_cast<double*>(&msg.LFFCbetaij[0]), const_cast<double*>(&msg.LFFCgammaij[0]), const_cast<double*>(&msg.optimizerResult[0]));
	
  opencv101::optimizerResult resultado;
  int i = 0;
    
  for (std::vector<double>::const_iterator it = msg.optimizerResult.begin(); it != msg.optimizerResult.end(); ++it) {
      double valor;
      valor = (*it);
      resultado.optimizerResult.push_back(valor);
      i++;
  }
    
  //resultado.optimizerResult.push_back([msg.optimizerResult[0]])
  //resultado.optimizerResult.push_back([msg.optimizerResult[1]])
  //resultado.optimizerResult.push_back([msg.optimizerResult[2]])
  pub.publish(resultado);
  
}

int main(int argc, char **argv)
{

	  
  ros::init(argc, argv, "optimizerCaller");

  /**
   * NodeHandle is the main access point to communications with the ROS system.
   * The first NodeHandle constructed will fully initialize this node, and the last
   * NodeHandle destructed will close down the node.
   */
  ros::NodeHandle n;

  /**
   * The subscribe() call is how you tell ROS that you want to receive messages
   * on a given topic.  This invokes a call to the ROS
   * master node, which keeps a registry of who is publishing and who
   * is subscribing.  Messages are passed to a callback function, here
   * called chatterCallback.  subscribe() returns a Subscriber object that you
   * must hold on to until you want to unsubscribe.  When all copies of the Subscriber
   * object go out of scope, this callback will automatically be unsubscribed from
   * this topic.
   *
   * The second parameter to the subscribe() function is the size of the message
   * queue.  If messages are arriving faster than they are being processed, this
   * is the number of messages that will be buffered up before beginning to throw
   * away the oldest ones.
   */
  ros::Subscriber sub = n.subscribe("/optimizerPub", 1, chatterCallback);
  //pub = n.advertise<opencv101::optimizerResult>("resultadosOtimizador", 10);
  pub = n.advertise<opencv101::optimizerResult>("/resultadosOtimizador", 10);
  /**
   * ros::spin() will enter a loop, pumping callbacks.  With this version, all
   * callbacks will be called from within this thread (the main one).  ros::spin()
   * will exit when Ctrl-C is pressed, or the node is shutdown by the master.
   */
  ros::spin();

  return 0;
}
