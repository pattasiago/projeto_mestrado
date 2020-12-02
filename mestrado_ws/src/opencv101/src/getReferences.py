#!/usr/bin/env python
# encoding: utf-8

import rospy
import math
#import numpy as np
from std_msgs.msg import String

# para importar os tipos de mensagens criados por mim
from opencv101.msg import *

class referenceGetter(object):

  def __init__(self):
      rospy.Subscriber("/camera_params", imageParams, self.callback)
      # mudar nome dos nós /parametros_de_controle/desvio_da_curvatura
      # definir tipo da saída
     # self.publisher = rospy.Publisher("/desvio_da_curvatura", desvioParams, queue_size = 10)
      self.publisher = rospy.Publisher("/desvio_da_curvatura", desvioParams, queue_size = 10)

  def callback(self, data):
      # rospy.loginfo(rospy.get_caller_id() + "Posição relativa: %s", data.poseori)
      kde = 1100*2
      ktheta = 0.3978/3.75
      if data.poseori[0]>=320:
        de0 = math.sqrt((data.poseori[0]-320)*(data.poseori[0]-320)+(240-240)*(240-240))/kde
        # quem é ktheta?
        thetae0 = ktheta*(0.69 - math.atan2(240,data.poseori[0]))
      else:
        de0 = -math.sqrt((data.poseori[0]-320)*(data.poseori[0]-320)+(240-240)*(240-240))/kde
        thetae0 = ktheta*(0.69-math.atan2(240,data.poseori[0]))
      
      h = std_msgs.msg.Header()
      h.stamp = rospy.Time.now()

      # publisher.publish(de0, tetae, c)
      print("Curvature: {}".format(data.circParams[2]))
      rospy.logerr("De0: {}, thetae0: {}, curvature: {}".format(de0, thetae0, data.circParams[2]))
      self.publisher.publish(h, de0, thetae0, data.circParams[2])
      
                

def main(args):
  ref = referenceGetter()
  rospy.init_node('get_references', anonymous=True)
  try:
    rospy.spin()
  except KeyboardInterrupt:
    print("Shutting down")
    cv2.destroyAllWindows()


if __name__ == '__main__':
    main(sys.argv)
