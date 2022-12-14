// // -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// // we only include RcppEigen.h which pulls Rcpp.h in for us
// #include <RcppEigen.h>
// #include <random>
//
// void log_sum_exp(Eigen::Map<Eigen::VectorXd> & probs,
//                  Eigen::Map<Eigen::VectorXd> & log_probs){
//
//   size_t resolution = probs.size();
//
//   double max_log_prob = log_probs.maxCoeff();
//
//   for (int i = 0; i < resolution; i++){
//     probs(i) = exp(log_probs(i) - max_log_prob);
//   }
//
// }
//
//
// // Function to remove 2 columns of mixing matrix
// void remove_matrix_columns(Eigen::Map<Eigen::MatrixXd> & column_reduced_mixing_matrix,
//                            Eigen::MatrixXd A,
//                            int col1, int col2){
//   size_t Q = A.cols();
//   int keep_index = 0;
//   for (int i = 0; i < Q; i++){
//     if (i != col1 & i != col2){
//       column_reduced_mixing_matrix.col(keep_index) = A.col(i);
//       keep_index++;
//     }
//   }
// }
//
//
//
// // calculate log probability of rotation
// double evaluate_rotation_log_probability(Eigen::Matrix2d A, Eigen::Matrix2d B,
//                                          double cos_theta, double sin_theta){
//   double log_prob = 0.0;
//
//   log_prob = A(0, 0) * B(0, 0) * pow(cos_theta, 2) + A(0, 1) * B(0,0) * cos_theta * sin_theta  +
//     A(0, 1) * B(0,0) * sin_theta * cos_theta + A(1,1) * B(0,0) * pow(sin_theta, 2) +
//     A(0,0)*B(1,1)*pow(sin_theta, 2) - A(0,1)*B(1,1)*sin_theta*cos_theta -
//     A(0,1)*B(1,1)*sin_theta*cos_theta + A(1,1)*B(1,1)*pow(cos_theta, 2);
//
//   return(log_prob);
// }
//
//
//
// // Function to sample the 2 x 2 rotation matrix
// void sample_rotation_matrix(Eigen::Map<Eigen::MatrixXd> & rotation_matrix,
//                             Eigen::MatrixXd scpns,
//                             Eigen::Map<Eigen::MatrixXd> atilde,
//                             Eigen::Map<Eigen::MatrixXd> btilde,
//                             Eigen::VectorXd angle_range,
//                             Eigen::VectorXd sin_theta_vec,
//                             Eigen::VectorXd cos_theta_vec,
//                             Eigen::Map<Eigen::VectorXd> & log_probs,
//                             Eigen::Map<Eigen::VectorXd> & probs){
//
//   std::random_device rd;
//   std::mt19937 gen(rd());
//
//   size_t angle_resolution = angle_range.size();
//
//   double cos_theta = 0.0;
//   double sin_theta = 0.0;
//
//   // Calculate the log probability under each rotation
//   for (int i = 0; i < angle_resolution; i++){
//     cos_theta = cos_theta_vec(i);
//     sin_theta = sin_theta_vec(i);
//     log_probs(i) = evaluate_rotation_log_probability(atilde, btilde, cos_theta, sin_theta);
//     log_probs(i + angle_resolution) = log_probs(i);
//     // Take care of sign related differences in the log probabilities
//     log_probs(i)                    += scpns(0,0) * cos_theta + scpns(1,0) * sin_theta + scpns(0,1) * (-1.0) * sin_theta - scpns(1,1) * (-1.0) * cos_theta;
//     log_probs(i + angle_resolution) += scpns(0,0) * cos_theta + scpns(1,0) * sin_theta + scpns(0,1) * (1.0) * sin_theta - scpns(1,1) * (1.0) * cos_theta;
//   }
//
//   // log sum exp to get probabilities
//   log_sum_exp(probs, log_probs);
//   // Draw a random number, use to sample with weights given by above probabilities
//   std::uniform_real_distribution<double> distribution(0.0, 1.0);
//   double samp = distribution(gen);
//   //std::cout << "sample was: " << samp << std::endl;
//   // Associate the samp variable with an index
//   double sum_probs = 0.0;
//   int selection = 0;
//   for (int i = 0; i < 2*angle_resolution; i++){
//     sum_probs += probs(i);
//     if (sum_probs >= samp){
//       selection = i;
//       break;
//     }
//   }
//
//   double final_s = -1.0;
//   double final_angle = 0.0;
//   if (selection >= angle_resolution){
//     final_s = 1.0;
//     final_angle = angle_range(selection - angle_resolution);
//   } else {
//     final_angle = angle_range(selection);
//   }
//
//   // Create the corresponding rotation matrix
//   rotation_matrix(0,0) = cos(final_angle);
//   rotation_matrix(1,0) = sin(final_angle);
//   rotation_matrix(0,1) = final_s * sin(final_angle);
//   rotation_matrix(1,1) = final_s * (-1.0) * cos(final_angle);
//
// }
//
// //
// // TODO - different solution for Q = 2
// // A = matrix(rnorm(100), nrow=10, ncol=10)
// // gibbs_sample_mixing_matrix(A, matrix(0, nrow=10, ncol=10), matrix(0, nrow=10, ncol=10), seq(from = 0.0, to = 1.0, length.out = 100))
// // [[Rcpp::export]]
// Eigen::MatrixXd gibbs_sample_mixing_matrix(Eigen::Map<Eigen::MatrixXd> & A,
//                                      Eigen::MatrixXd H,
//                                      Eigen::MatrixXd B,
//                                      Eigen::VectorXd angle_range,
//                                      Eigen::VectorXd sin_theta,
//                                      Eigen::VectorXd cos_theta,
//                                      Eigen::Map<Eigen::MatrixXd> & column_reduced_mixing_matrix,
//                                      Eigen::Map<Eigen::MatrixXd> & selected_columns_mixing_matrix,
//                                      Eigen::Map<Eigen::MatrixXd> & selected_columns_proj_null_space,
//                                      Eigen::Map<Eigen::Matrix2d> & btilde,
//                                      Eigen::Map<Eigen::Matrix2d> & atilde,
//                                      Eigen::Map<Eigen::Matrix2d> & rotation_matrix,
//                                      Eigen::Map<Eigen::MatrixXd> & null_rotation,
//                                      Eigen::Map<Eigen::VectorXd> & log_probs,
//                                      Eigen::Map<Eigen::VectorXd> & probs,
//                                      Eigen::Map<Eigen::MatrixXd> & null_space) {
//
//   // Dimension of orthonormal matrix
//   size_t Q = A.rows();
//   size_t angle_resolution = angle_range.size();
//
//   // Get the order we will visit columns of mixing matrix
//   std::vector<int> index_order(Q);
//   for (int i = 0; i < Q; i++){
//     index_order[i] = i;
//   }
//   std::random_shuffle ( index_order.begin(), index_order.end() );
//
//   int col1 = 0;
//   int col2 = 0;
//
//   // define intermediate quantities
//   //Eigen::MatrixXd column_reduced_mixing_matrix(Q, Q-2);
//   //Eigen::MatrixXd selected_columns_mixing_matrix(Q, 2);
//   //Eigen::MatrixXd selected_columns_proj_null_space(2, 2);
//   //Eigen::Matrix2d btilde(2, 2);
//   //Eigen::Matrix2d atilde(2, 2);
//   //Eigen::Matrix2d rotation_matrix(2, 2);
//   //Eigen::MatrixXd null_rotation(Q, 2);
//   //Eigen::VectorXd log_probs(2*angle_resolution);
//   //Eigen::VectorXd probs(2*angle_resolution);
//
//   for (int i = 0; i < Q; i++){
//     for (int j = i+1; j < Q; j++){
//       col1 = index_order[i];
//       col2 = index_order[j];
//
//       // Extract all columns of the mixing matrix except the select pair
//       remove_matrix_columns(column_reduced_mixing_matrix, A, col1, col2);
//
//       // Get the null space
//       Eigen::FullPivLU<Eigen::MatrixXd> lu(column_reduced_mixing_matrix.transpose());
//       null_space = lu.kernel();
//
//       // Quantities for rotation update
//       btilde(1, 1) = B(col1, col1);
//       btilde(2, 2) = B(col2, col2);
//       atilde = null_space.transpose() * A * null_space;
//
//       selected_columns_mixing_matrix.col(0) = H.col(col1);
//       selected_columns_mixing_matrix.col(1) = H.col(col2);
//       selected_columns_proj_null_space = null_space.transpose() * selected_columns_mixing_matrix;
//
//       // Sample rotation matrix
//       sample_rotation_matrix(rotation_matrix, selected_columns_proj_null_space,
//                              atilde, btilde, angle_range, sin_theta, cos_theta, log_probs, probs);
//
//       // Apply to null space - columns are new columns of mixing matrix
//       null_rotation = null_space * rotation_matrix;
//
//       A.col(col1) = null_rotation.col(0);
//       A.col(col2) = null_rotation.col(1);
//
//     }
//   }
//
//   return A;
// }
//
//
//
