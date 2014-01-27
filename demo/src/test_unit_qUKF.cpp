#include <cv.h>
#include <highgui.h>
#include <lcpp_gmphd/gmphd_filter.h>


int main() {

  float pose = 0.0;

  // Instanciate the motion filter
  motion_estimator = new MotionEstimation(&pose, 1.0, 1.0, 1.0).

  // Deal with the OpenCV window..
  float angle = CV_PI/2-0.03;
  int width   = 800;
  int height  = 800;

  IplImage * image = cvCreateImage(cvSize(width,height),8,3);

  for(;;angle += 0.01) {
	  cvZero(image);
	  //int nb_target_all = 1+rand()%(max_target-1);
	  int nb_target_all = 1;

	  vector<float> measurements;
	  measurements.resize(dim_measure*nb_target_all);
	  int nb_target = rand()%nb_target_all;
	  for (int i=0; i< nb_target>>1; ++i) {
		  int x = (width>>1)  + 300*cos(angle) + (rand()%2==1?-1:1)*(rand()%10);
		  int y = (height>>1) + 300*sin(angle) + (rand()%2==1?-1:1)*(rand()%10);
		  measurements[dim_measure*i]   = x;
		  measurements[dim_measure*i+1] = y;
		  cvDrawCircle(image,cvPoint(x,y),2,cvScalar(0,0,255),2);
	  }
	  for (int i=nb_target>>1; i< nb_target_all; ++i) {
		  int x = (width>>1)  + 300*cos(angle+CV_PI) + (rand()%2==1?-1:1)*(rand()%10);
		  int y = (height>>1) + 300*sin(angle+CV_PI) + (rand()%2==1?-1:1)*(rand()%10);
		  measurements[dim_measure*i]   = x;
		  measurements[dim_measure*i+1] = y;
		  cvDrawCircle(image,cvPoint(x,y),2,cvScalar(0,0,255),2);
	  }

	  gmphd.setNewMeasurements(measurements);
	  gmphd.propagate();
	  vector<float> weight;
	  vector<float> position;
	  vector<float> speed;
	  gmphd.getTrackedTargets(gmphd_threshold,position,speed,weight);


	  birth_model.clear();
	  for (int i=0; i< position.size(); i+=2) {
// 			GaussianModel birth_gaussian;
// 			birth_gaussian.cov       = MatrixXf::Identity(dim_state,dim_state);
// 			birth_gaussian.mean      = MatrixXf::Zero(dim_state,1);
// 			birth_gaussian.mean(0,0) = position[i];
// 			birth_gaussian.mean(1,0) = position[i+1];
// 			birth_model.push_back(birth_gaussian);
		  cvDrawCircle(image,cvPoint(position[i],position[i+1]),5,cvScalar(0,255,0),2);
	  }
	  //gmphd.setBirthModel(birth_model);

	  cvShowImage("image",image);
	  printf("-----------------------------------------------------------------\n");
	  if( cvWaitKey(5) == 27)
		  break;
  }
  cvReleaseImage(&image);
  return 1;
}
