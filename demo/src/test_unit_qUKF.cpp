#include <cv.h>
#include <highgui.h>
#include <motion_filtering.h>
#include <circ_buffer.h>
#include <utility>
#include <memory>

using namespace cv;
using namespace std;

void draw(IplImage *image,
          float const  *measurements,
          float const *predicted_state,
          CircBuffer<Point2f> const & filtered_states) {

    cvDrawCircle(image, cvPoint(measurements[0],measurements[1]),2,cvScalar(0,0,255),2);

    Point2f start = filtered_states[0];

    for (unsigned int i=1; i<filtered_states.size(); ++i)
    {
        auto const & pose = filtered_states[i];
        cvDrawLine(image, cvPoint(start.x, start.y), cvPoint(pose.x,pose.y), cvScalar(0,255,0),2);
        start = pose;
    }

    // Show the predicted motion vector (?)
    auto const & lastVal = filtered_states.front();

    cvDrawLine( image,
                cvPoint(lastVal.x, lastVal.y),
                cvPoint(predicted_state[0],predicted_state[1]), cvScalar(255,0,0),2);
}

bool updateDisplay(IplImage *image, int waitTime)
{
    cvShowImage("image",image);

    int k = cvWaitKey(waitTime);

    if ((k == 27) || (k == 1048603))
    {
        return false;
    }
    else if (k != -1)
    {
        printf("Key pressed : %d\n Nothing to do", k);
    }
    return true;
}

int main() {

    // Deal with the OpenCV window..
    float angle = 0.f, angularStep = 0.01f;
    unsigned int width = 800, height  = 800;

    int const N_TARGETS = 5;
    int const MAX_POSES = 30;
    int const WAIT_TIME_MS = 40;

    // Allocate all the buffers to store historical data
    vector<CircBuffer<Point2f> *> vec_poses(N_TARGETS);
    for (int i=0; i< N_TARGETS; ++i)
    {
        vec_poses[i] = new CircBuffer<Point2f>(MAX_POSES, true);
    }

    IplImage * image = cvCreateImage(cvSize(width,height),8,3);

    // Instanciate the motion filters
    float poses[2] = {0, 0};
    vector<unique_ptr<MotionEstimation>> motion_estimators(N_TARGETS);

    for (auto & motion_estimator :  motion_estimators)
    {
        motion_estimator.reset( new MotionEstimation(&poses[0], 1e3, 0.1, 0.05));
    }

    // Track the circling targets
    Eigen::Vector2f measurements, filtered_state, predicted_state, previous_state;
    int t;

    for(;;angle += angularStep) {
        cvZero(image);
        t = 0;

        // Create a new measurement for every target, and do the update
        for (auto & motion_estimator : motion_estimators)
        {
            motion_estimator->getStatePost(previous_state);

            // Propagate the previous state
            motion_estimator->predict();

            // Get the predicted state :
            motion_estimator->getStatePost(predicted_state);

            // Update the state with a new noisy measurement :
            measurements = { (width>>1)  + 300*cos(angle) + (rand()%2==1?-1:1)*(rand()%30),
                             (height>>1) + 300*sin(angle) + (rand()%2==1?-1:1)*(rand()%30) };

            // Get the filtered state :
            motion_estimator->update(measurements);
            motion_estimator->getStatePost(filtered_state);
            vec_poses[t]->add(Point2f( filtered_state(0), filtered_state(1)));

            // Draw both the noisy input and the filtered state :
            draw(image, measurements.data(), predicted_state.data(), *vec_poses[t++]);
        }

        if(!updateDisplay(image, WAIT_TIME_MS))
        {
            break;
        }
    }

    cvReleaseImage(&image);
    return 1;
}
