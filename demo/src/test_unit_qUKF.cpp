#include <cv.h>
#include <highgui.h>
#include <motion_filtering.h>

int main() {

  // Deal with the OpenCV window..
  float angle = CV_PI/2-0.03;
  unsigned int width   = 800;
  unsigned int height  = 800;

  int n_targets = 5;

  IplImage * image = cvCreateImage(cvSize(width,height),8,3);

  // Instanciate the motion filters
  float poses[3] = {0,0,0};

  MotionEstimation **motion_estimators;
  motion_estimators = new MotionEstimation*[n_targets];

  for (unsigned int i=0; i<n_targets; ++i) {
      motion_estimators[i] = new MotionEstimation(&poses[i], 10.0, 1.0, 1.0);
    }

  // Track the circling targets
  float measurements[3];
  float filtered_state[3];

  for(;;angle += 0.01) {
      cvZero(image);

      // Create a new measurement for every target, and do the update
      for (unsigned int i=0; i< n_targets; ++i) {
          // Propagate the previous state
          motion_estimators[i]->predict();

          // Update the state with a new noisy measurement :
          measurements[0] = (width>>1)  + 300*cos(angle) + (rand()%2==1?-1:1)*(rand()%50);
          measurements[1] = (height>>1) + 300*sin(angle) + (rand()%2==1?-1:1)*(rand()%50);

          motion_estimators[i]->update(measurements);

          // Get the filtered state :
          motion_estimators[i]->getLatestSate(filtered_state);

          // Draw both the noisy input and the filtered state :
          cvDrawCircle(image,cvPoint(measurements[0],measurements[1]),2,cvScalar(0,0,255),2);
          cvDrawCircle(image,cvPoint(filtered_state[0],filtered_state[1]),2,cvScalar(0,255,0),2);
        }


      // Show this stuff
      cvShowImage("image",image);
      printf("-----------------------------------------------------------------\n");
      int k = cvWaitKey(33);

      if ((k == 27) || (k == 1048603))
        break;
      else if (k != -1)
        printf("Key pressed : %d\n", k);
    }


  // Close everything and leave
  cvReleaseImage(&image);
  for (int i=0; i<n_targets; ++i) {
      delete motion_estimators[i];
    }

  delete []motion_estimators;

  return 1;
}
