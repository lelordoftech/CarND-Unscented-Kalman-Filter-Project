#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "measurement_package.h"
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[])
{
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 3)
  {
    has_valid_args = true;
    //cerr << usage_instructions << endl;
  }
  else if (argc > 3)
  {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args)
  {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name)
{
  if (!in_file.is_open())
  {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open())
  {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[])
{
  check_arguments(argc, argv);

  string in_file_name_ = "../data/obj_pose-laser-radar-synthetic-input.txt"; //argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = "output.txt"; //argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  /**********************************************
   *  Set Measurements                          *
   **********************************************/

  vector<MeasurementPackage> measurement_pack_list;
  vector<VectorXd> gt_pack_list;

  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line))
  {
    string sensor_type;
    MeasurementPackage meas_package;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0)
    {
      // laser measurement
      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }
    else if (sensor_type.compare("R") == 0)
    {
      // radar measurement
      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

      // read ground truth data to compare later
      float x_gt;
      float y_gt;
      float vx_gt;
      float vy_gt;
      iss >> x_gt;
      iss >> y_gt;
      iss >> vx_gt;
      iss >> vy_gt;
      VectorXd gt_values = VectorXd(4);
      gt_values << x_gt, y_gt, vx_gt, vy_gt;
      gt_pack_list.push_back(gt_values);
  }

  /**********************************************
   *  Measurements                              *
   **********************************************/

  // compute the accuracy (RMSE)
  Tools tools;

  // Create a UKF instance
  UKF ukf;
  float start_point = atof(argv[1]);
  float step = atof(argv[2]);
  vector<float> rmse_list;

  for (int test_id = 0; test_id < 1; test_id++)
  {
    ukf.x_ = VectorXd::Zero(5);
    ukf.P_ = MatrixXd::Zero(5, 5);
    ukf.std_a_ = 0.292;
    ukf.std_yawdd_ = 0.261;
    ukf.is_initialized_ = false;
    ukf.n_x_ = 5;
    ukf.n_aug_ = 7;
    ukf.Xsig_pred_ = MatrixXd::Zero(ukf.n_x_, 2 * ukf.n_aug_ + 1);
    ukf.weights_ = VectorXd::Zero(2 * ukf.n_aug_ + 1);
    ukf.P_ << 0.422, 0,    0,     0,     0,
              0,     0.48, 0,     0,     0,
              0,     0,    3.159, 0,     0,
              0,     0,    0,     0.033, 0,
              0,     0,    0,     0,     4.3;

    //ukf.P_(2,2) = start_point;

    // used to compute the RMSE later
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;

    // start filtering from the second frame (the speed is unknown in the first
    // frame)

    size_t number_of_measurements = measurement_pack_list.size();

    for (size_t k = 0; k < number_of_measurements; ++k)
    {
      // Call the UKF-based fusion
      ukf.ProcessMeasurement(measurement_pack_list[k]);

      // convert ukf x vector to cartesian to compare to ground truth
      VectorXd ukf_x_cartesian_ = VectorXd(4);

      float x_estimate_ = ukf.x_(0);
      float y_estimate_ = ukf.x_(1);
      float vx_estimate_ = ukf.x_(2) * cos(ukf.x_(3));
      float vy_estimate_ = ukf.x_(2) * sin(ukf.x_(3));

      ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

      estimations.push_back(ukf_x_cartesian_);
      ground_truth.push_back(gt_pack_list[k]);
    }
    VectorXd rmse = tools.CalculateRMSE(estimations, ground_truth);
    rmse_list.push_back(rmse[1]);
    cout << "RMSE: " << rmse[0] << " " << rmse[1] << " " << rmse[2] << " " << rmse[3] << endl;
    cout << " " << test_id << std::flush;

    start_point += step;
  }
  cout << endl;

  vector<float>::iterator min_index = min_element(rmse_list.begin(), rmse_list.end());
  cout << "min value at " << distance(rmse_list.begin(), min_index) << endl;
  cout << "min value " << *min_index << endl;
  cout << "set value " << atof(argv[1]) + step*distance(rmse_list.begin(), min_index) << endl;

  // close files
  if (out_file_.is_open())
  {
    out_file_.close();
  }

  if (in_file_.is_open())
  {
    in_file_.close();
  }

  cout << "Done!" << endl;
  return 0;
}
