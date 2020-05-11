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

  if(is_initialized) {
    return;
  }

  num_particles = 20;  // TODO: Set the number of particles
  
  // random number generator
  std::default_random_engine gen;
  // creates normal (Gaussian) distributions around (x, y, theta) initial position
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  for(int i = 0; i < num_particles; ++i) {
    Particle p;
    p.id = i+1;         // number particles from 1, like landmarks
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

  // std::cout << "v: " << velocity << "  yaw_rate: " << yaw_rate << std::endl;

  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);
  
  double v_theta_dot = 0.0;
  double v_dt = velocity * delta_t;
  double theta_dot_dt = yaw_rate * delta_t;
  
  
  // for(auto& p : particles) {
  for(int i = 0; i < num_particles; ++i) {
    Particle& p = particles[i];
    // if(fabs(yaw_rate) > std::numeric_limits<double>::epsilon()) {
    if(fabs(yaw_rate) < 0.00001) {
      p.x += v_dt * cos(p.theta);
      p.y += v_dt * sin(p.theta);
    } else {
      p.x += v_theta_dot * (sin(p.theta + theta_dot_dt) - sin(p.theta));
      p.y += v_theta_dot * (cos(p.theta) - cos(p.theta + theta_dot_dt));
      p.theta += theta_dot_dt;
    }
    
    // add random Gaussian noise
    p.x += dist_x(gen);
    p.y += dist_y(gen);
    p.theta += dist_theta(gen);
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

  // for(auto& obs : observations) {
  for(int i = 0; i < observations.size(); ++i) {
    LandmarkObs& obs = observations[i];
    // compute distance from observation to each landmark
    vector<double> distances;
    // for(auto& pred : predicted) {
    for(int j = 0; j < predicted.size(); ++j) {
      distances.push_back(dist(obs.x, obs.y, predicted[j].x, predicted[j].y));
    }
    // get an iterator pointing to the element with the shortest distance
    auto iter = std::min_element(distances.begin(), distances.end());
    // get the index of the element
    int index = std::distance(distances.begin(), iter);
    // set the id to the id of the nearest landmark
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

  std::cout << "---> updateWeights()" << std::endl;
  std::cout << "sensor_range: " << sensor_range << std::endl;
  std::cout << "n observations: " << observations.size() << std::endl;
  for(auto o : observations) {
    std::cout << "(x, y): (" << o.x << ", " << o.y << ")" << std::endl;
  }

  // compute some intermediate results once
  double sigma_x = std_landmark[0];
  double sigma_y = std_landmark[1];
  double sigma_2x2 = 2 * sigma_x * sigma_x;
  double sigma_2y2 = 2 * sigma_y * sigma_y;
  double normalizer = 1.0 / (2 * M_PI * sigma_x * sigma_y);

  weights.clear();

  // for(auto& p : particles) {
  for(int i = 0; i < num_particles; ++i) {
    Particle& p = particles[i];
    // find all the landmarks within sensor range
    vector<LandmarkObs> landmarks;
    for(auto lm : map_landmarks.landmark_list) {
      if(dist(p.x, p.y, lm.x_f, lm.y_f) <= sensor_range) {
        landmarks.push_back(LandmarkObs{lm.id_i, lm.x_f, lm.y_f});
      }
    }

    // transform observations into map coordinates
    double cos_theta = cos(p.theta);
    double sin_theta = sin(p.theta);
    vector<LandmarkObs> obs_transform;
    for(auto o : observations) {
      double xm = p.x + (cos_theta * o.x) - (sin_theta * o.y);
      double ym = p.y + (sin_theta * o.x) + (cos_theta * o.y);
      obs_transform.push_back(LandmarkObs{o.id, xm, ym});
    }
    // get associations
    dataAssociation(landmarks, obs_transform);
    
    double prob = 1.0;
    double pdf;
    for(auto o : obs_transform) {
      // get the associated landmark
      auto lm = map_landmarks.landmark_list[o.id-1];
      // compute pdf for the observation and accumulate product
      pdf = normalizer * exp(-(pow(o.x-lm.x_f, 2)/sigma_2x2
                               + pow(o.y-lm.y_f, 2)/sigma_2y2));
      prob *= pdf;
    }
    // update weight with pdf
    p.weight = prob;
    weights.push_back(p.weight);
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<Particle> new_particles;

  std::default_random_engine gen;
  double beta = 0.0;
  double max_wt;

  max_wt = *std::max_element(weights.begin(), weights.end());

  std::uniform_int_distribution<> random_int(0, num_particles-1);
  std::uniform_real_distribution<> beta_incr(0.0, 2.0 * max_wt);

  int index = random_int(gen);
  for(int i=0; i<num_particles; ++i) {
    beta += beta_incr(gen);
    while(beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }
  particles.clear();
  for(auto& np : new_particles) {
    particles.push_back(np);
  }
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
  particle.associations = associations;
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