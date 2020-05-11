/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50;  // TODO: Set the number of particles
  std::default_random_engine gen;
  
  // creates normal (Gaussian) distributions for x initial position
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std_theta[2]);
  
  for(int i = 0; i < num_particles; ++i) {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    
    particles.push_back(p);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  
  double v_theta_dot = 0.0;
  double v_dt = velocity * delta_t;
  double theta_dot_dt = yaw_rate * delta_t;
  bool zero_theta_dot = true;
  
  if(yaw_rate > std::numeric_limits<double>::epsilon()) {
    v_theta_dot = velocity/yaw_rate;
    zero_theta_dot = false;
  }
  
  for(auto& p : particles) {
    if(zero_theta_dot) {
      p.x += v_dt * cos(p.theta);
      p.y += v_dt * sin(p.theta);
    } else {
      p.x += v_theta_dot * (sin(p.theta + theta_dot_dt) - sin(p.theta));
      p.y += v_theta_dot * (cos(p.theta) - cos(p.theta + theta_dot_dt));
      p.theta += theta_dot_dt;
    }
    
    normal_distribution<double> dist_x(p.x, std_pos[0]);
    normal_distribution<double> dist_y(p.y, std_pos[1]);
    normal_distribution<double> dist_theta(p.theta, std_pos[2]);
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for(auto& obs : observations) {
    vecton<double> distances;
    for(auto& pred : predicted) {
      distances.pushback(dist(pred.x, pred.y, obs.x, obs.y));
    }
    auto shortest_dist = std::min_element(distances.begin(), distance.end());
    int index = std::distance(distances.begin(), shortest_dist);
    obs.id = predicted[index].id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  // to use dataAssociation, we need a vector of landmarks in the format of LandmarkObs
  vector<LandmarkObs> landmarks;
  for(auto& lm : map_landmarks.landmark_list) {
    LandmarkObs lm_temp;
    lm_temp.id = lm.id_i;
    lm_temp.x = lm.x_f;
    lm_temp.y = lm.y_f;
    landmarks.pushback(lm_temp);
  }

  double sigma_x = std_landmark[0];
  double sigma_y = std_landmark[1];
  double sigma_x2 = 2 * sigma_x * sigma_x;
  double sigma_y2 = 2 * sigma_y * sigma_y;
  double normalizer = 1.0 / (2 * M_PI * sigma_x * sigma_y);

  for(auto& p : particles) {
    // transform observations into map coordinates
    double cos_theta = cos(p.theta);
    double sin_theta = sin(p.theta);
    vector<LandmarkObs> observations_map;
    for(auto &obs : observations) {
      double xm = p.x + (cos_theta * obs.x) - (sin_theta * obs.y);
      double ym = p.y + (sin_theta * obs.x) + (cos_theta * obs.y);
      LandmarkObs obs_temp;
      obs_temp.id = 0;
      obs_temp.x = xm;
      obs_temp.y = xy;
      observations_map.push_back(obs_temp);
      p.sense_x.push_back(xm);
      p.sense_y.push_back(ym);
    }
    // get associations
    dataAssociation(landmarks, observations_map);
    double prob = 1.0
    for(auto& o : observations_map) {
      p.associations.push_back(o.id)
      auto lm = map_landmarks.landmark_list[o.id-1];
      prob *= normalizer * exp(-( pow(o.x-lm.x_f, 2)/sigma_x2
                               + pow(o.y-lm.y_f, 2)/sigma_y2));
    }
    // update weight with pdf
    p.weight = prob;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}