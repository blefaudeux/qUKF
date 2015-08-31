#include <cv.h>
#include <highgui.h>
#include <motion_filtering.h>
#include <circ_buffer.h>
#include <utility>

void draw(IplImage *image,
          const float *measurements,
          const float *predicted_state,
          CircBuffer<std::pair<float, float>> const & filtered_states) {

    cvDrawCircle(image,cvPoint(measurements[0],measurements[1]),2,cvScalar(0,0,255),2);

    vector<std::pair<float, float>> const & buffer = filtered_states.buffer();
    for (auto const & pose : buffer)
    {
        cvDrawCircle(image,cvPoint(pose.first,pose.second),2,cvScalar(0,255,0),2);
    }

    // Show the predicted motion vector (?)
    cvDrawLine(image, cvPoint(buffer[0].first, buffer[1].second), cvPoint(predicted_state[0],predicted_state[1]), cvScalar(255,0,0),2);
}


int main() {

    // Deal with the OpenCV window..
    float angle = CV_PI/2-0.03;
    unsigned int width   = 800, height  = 800;
    int const n_targets = 5;
    int const MAX_POSES = 10;

    vector<CircBuffer<std::pair<float, float>> *> vec_poses(n_targets);
    for (int i=0; i< n_targets; ++i)
    {
        vec_poses[i] = new CircBuffer<std::pair<float, float>>(MAX_POSES);
    }

    IplImage * image = cvCreateImage(cvSize(width,height),8,3);

    // Instanciate the motion filters
    float poses[3] = {0,0,0};

    vector<MotionEstimation *> motion_estimators(n_targets);

    int i=0;
    for (auto & motion_estimator :  motion_estimators)
    {
        motion_estimator = new MotionEstimation(&poses[i++], 1e3, 1.0, 0.05);
    }

    // Track the circling targets
    float measurements[6] = {0,0,0,0,0,0};
    float filtered_state[6], predicted_state[6], previous_state[6];
    int t;

    for(;;angle += 0.01) {
        cvZero(image);
        t= 0;

        // Create a new measurement for every target, and do the update
        for (auto & motion_estimator : motion_estimators)
        {
            motion_estimator->getLatestState(previous_state);

            // Propagate the previous state
            motion_estimator->predict();

            // Get the predicted state :
            motion_estimator->getPropagatedState(predicted_state);

            // Update the state with a new noisy measurement :
            measurements[0] = (width>>1)  + 300*cos(angle) + (rand()%2==1?-1:1)*(rand()%50);
            measurements[1] = (height>>1) + 300*sin(angle) + (rand()%2==1?-1:1)*(rand()%50);

            // Define a new 'speed' measurement
            measurements[3] = measurements[0] - previous_state[0];
            measurements[4] = measurements[1] - previous_state[1];

            motion_estimator->update(measurements);

            // Get the filtered state :
            motion_estimator->getLatestState(filtered_state);
            vec_poses[t]->add(std::make_pair<float, float>(float(filtered_state[0]), float(filtered_state[1])));

            // Draw both the noisy input and the filtered state :
            draw(image, measurements, predicted_state, *vec_poses[t]);
        }

        // Show this stuff
        cvShowImage("image",image);
        printf("-----------------------------------------------------------------\n");
        int k = cvWaitKey(20);

        if ((k == 27) || (k == 1048603))
            break;
        else if (k != -1)
            printf("Key pressed : %d\n", k);
    }


    // Close everything and leave
    cvReleaseImage(&image);
    for (auto & motion_estimator : motion_estimators) {
        delete motion_estimator;
    }

    return 1;
}
