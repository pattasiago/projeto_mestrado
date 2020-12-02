#!/usr/bin/env python
# encoding: utf-8

import rospy
import sys
import math
import tf
import numpy as np
import os 
import numpy.ctypeslib as ctl
import ctypes
from opencv101.msg import *
from nav_msgs.msg import Odometry
from geometry_msgs.msg import Twist, PointStamped


class Robot():
  VFCMode = 0
  qsip = 0
  formationK1 = 0
  formationK3 = 0 
  formationQ = 0
  formationP31 = 0
  nmpcVarsWs = 0
  nmpcVarsWu = 0
  nmpcVarsP = 0
  nmpcVarsAlpha = 0
  nmpcVarsBeta = 0

   
  def __init__(self):
    rospy.Subscriber("/desvio_da_curvatura_throttle", desvioParams, self.callback)
    rospy.Subscriber("/husky_velocity_controller/odom", Odometry, self.odometryCb) 
    rospy.Subscriber("/resultadosOtimizador", optimizerResult, self.optimizerCb)
    rospy.Subscriber("/closest_point", PointStamped, self.calcdist)
    self.optimizerPub = rospy.Publisher("/optimizerPub", optimizerParams, queue_size = 10)
    self.cmdVelPub = rospy.Publisher("/cmd_vel", Twist, queue_size=1)  
    self.parametrosVelocidade = Twist()
    self.vRef = 0.2 
    self.wRef = 0.0 
    self.horizontePredicao = 3
    self.horizonteControle = 3
    self.horizonte = 1.7
    self.Ts = 0.2
    
    self.Q = 0.1
    self.R = 0.001
    self.id = 0
    self.useVisualControl = True
    
    self.otimizerResult = [0, 0, 0]
    
    self.LFFClij = [0, 0, 0]
    self.LFFCpsiij = [0, 0, 0]
    self.LFFCbetaij = [0, 0, 0]
    self.LFFCgammaij = [0, 0, 0] 
    
    self.xe = 0
    self.u_opt = 0
    self.errorModelDep = 0
    self.x = [0, 0, 0, 0, 0, 0, 0]
    self.u = [0, 0, 0, 0, 0, 0, 0]
    self.de = [0, 0, 0]
    self.thetae = [0, 0, 0]

    self.throttleCounter = 0

    self.dist = 0
    
  def preditor(self, data):    
    for i in range(self.horizontePredicao-1):      
      self.de[i+1] = self.de[0] + self.Ts * (self.horizonte * self.wRef + (self.vRef + self.de[i] * self.wRef) * math.tan(self.thetae[i]))
      self.thetae[i+1] = self.thetae[i] + self.Ts * self.u[i+1]   

      for i in range(self.horizontePredicao):
        self.u[i] = self.wRef - self.circunferenceCurvature * (self.vRef + self.wRef * self.de[i]/math.cos(self.thetae[i]))
        if i % 2 == 0:
          self.x[i] = self.de[i]
        else:
          self.x[i] = self.thetae[i]
    
    parameters = [self.horizontePredicao, self.horizonteControle, self.Q, self.R, self.Ts, self.circunferenceCurvature, Robot.qsip, self.id, self.vRef, self.v, self.wRef, self.formationK1, self.formationK3, self.formationQ, self.formationP31, Robot.nmpcVarsWs, Robot.nmpcVarsWu, self.wRef, self.horizonte, self.horizonte, self.useVisualControl, Robot.nmpcVarsP, Robot.nmpcVarsAlpha, Robot.VFCMode, self.xe, self.u_opt, self.thetae[0], Robot.nmpcVarsBeta, self.errorModelDep, Robot.nmpcVarsBeta, self.v] 
    
    h = std_msgs.msg.Header()
    h.stamp = rospy.Time.now()
    
    self.optimizerPub.publish(h, parameters, self.x, self.u, self.LFFClij, self.LFFCpsiij, self.LFFCbetaij, self.LFFCgammaij, self.otimizerResult);  

  def callback(self, data):
    self.circunferenceCurvature = data.curvature
    self.de[0] = data.de0
    self.thetae[0] = data.thetae0
    self.preditor(data)
    self.cmdVelPub.publish(self.parametrosVelocidade)
   
  def odometryCb(self, data):
    self.v = data.twist.twist.linear.x
    self.w = data.twist.twist.angular.z

  def optimizerCb(self, data):
    jMin = data.optimizerResult[0]
    self.u_opt = data.optimizerResult[1]
    self.wRef = (self.u_opt*math.cos(self.thetae[0])+self.circunferenceCurvature*self.vRef)/((math.cos(self.thetae[0])-self.circunferenceCurvature*self.de[0]))
    self.parametrosVelocidade.linear.x = self.vRef
    self.parametrosVelocidade.angular.z = self.wRef
    self.cmdVelPub.publish(self.parametrosVelocidade)

  def calcdist(self, data):
    self.dist=math.sqrt(math.pow(data.point.x,2)+math.pow(data.point.y,2))
    rospy.logerr("Distance: {}".format(self.dist))
      
def main(args):
  rospy.init_node('preditor', anonymous=True)
  robot = Robot()
  try:
    rospy.spin()
  except KeyboardInterrupt:
    print("Shutting down")
    cv2.destroyAllWindows()


if __name__ == '__main__':
    main(sys.argv)

