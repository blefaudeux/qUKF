#include <cv.h>
#include <highgui.h>
#include <lcpp_gmphd/gmphd_filter.h>


int main() {
	int   max_target              = 10;
	int   max_prune               = max_target;
	int   dim_state               = 2;
	int   dim_measure             = 2;
	float sampling                = 1.f;
	float process_noise           = 10;//1.f;
	float detection_probability   = 0.8f;
	float pose_measurement_noise  = 10;//0.5f;
	float background_measurement  = 0.00000001; // (proba fausse détection)²
	float probability_of_survival = 0.99;//0.9f;
	float gmphd_threshold         = 0.1;//0.005f;
	float prune_truncature_thld   = 0.1;//0.01f;
	float prune_merge_thld        = 4;//2.0f;

	vector<float> measurement_noise;
	measurement_noise.push_back(pose_measurement_noise);
	measurement_noise.push_back(pose_measurement_noise);

	GmphdFilter gmphd(max_target,dim_state,dim_measure,false,0);
	gmphd.reset();
	gmphd.setDynamicsModel(sampling,process_noise);
	gmphd.setObservationModel(detection_probability, measurement_noise, background_measurement);
	gmphd.setPruningParameters(prune_truncature_thld,prune_merge_thld,max_prune);
	gmphd.setSurvivalProbability(probability_of_survival);
	vector<SpawningModel> spawn_model; 
	gmphd.setSpawnModel(spawn_model);


	float angle = CV_PI/2-0.03;
	int width  = 800;
	int height = 800;
	vector<GaussianModel> birth_model;
	{
		GaussianModel birth_gaussian(dim_state);// recouvrement des zone dapparition la plus probable.
		birth_gaussian.cov         = 60000*MatrixXf::Identity(dim_state,dim_state);
		//birth_gaussian.mean      = MatrixXf::Zero(dim_state,1);
   		birth_gaussian.mean(0,0) = (width>>1) ;// + 300*cos(angle);
   		birth_gaussian.mean(1,0) = (height>>1);// + 300*sin(angle);
 		birth_gaussian.weight    = 0.1f;//birth_weight[i]; // Uniform_spread ^ dims
		birth_model.push_back(birth_gaussian);

	}
	gmphd.setBirthModel(birth_model);

	IplImage * image = cvCreateImage(cvSize(width,height),8,3);
	for(;;angle += 0.01) {
		cvZero(image);
		int nb_target_all = 1+rand()%(max_target-1);
		//printf("nb_target measure = %d\n\n",nb_target_all);
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