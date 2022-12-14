// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <random>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//



// and the inner product returns a scalar
//
// [[Rcpp::export]]
Eigen::MatrixXd eigen_matrix_prod_XYt(const Eigen::Map<Eigen::MatrixXd> & X,
                                      const Eigen::Map<Eigen::MatrixXd> & Y) {
    Eigen::MatrixXd v =  X * Y.transpose();
    return v;
}



// [[Rcpp::export]]
void utility_update_matrix_across_nodes_inplace(Eigen::Map<Eigen::MatrixXd> & X,
                                                const Rcpp::List & Xadd) {

  Eigen::MatrixXd test = Xadd[0];
  X *= 0.0; // set the input matrix to zero
  for (int i = 0; i < Xadd.length(); i++){
    test = Xadd[i];
    X += test;
  }
}


// [[Rcpp::export]]
void matmul_inplace(Eigen::Map<Eigen::MatrixXd> & mat_to_update,
    const Eigen::Map<Eigen::MatrixXd> & left_side,
    const Eigen::Map<Eigen::MatrixXd> & right_side){

  mat_to_update = left_side * right_side;

}


//
// Y  is V x NQ
// Si is V x NQ
// [[Rcpp::export]]
Eigen::MatrixXd calculate_YiSit(const Eigen::Map<Eigen::MatrixXd> & Y,
                                const Eigen::Map<Eigen::MatrixXd> & Si, int Q) {
    int N = Y.cols() / Q;
    int V = Y.rows() / Q;
    int sInd = -1;
    int eInd = -1;
    Eigen::MatrixXd YSi = Eigen::MatrixXd::Zero(Q*N, Q);
    for (int i = 0; i < N; i++){
        sInd = eInd + 1;
        eInd += Q;
        YSi.block(sInd, 0, Q, Q) = Y.block(0, sInd, V, Q).transpose() * Si.block(0, sInd, V, Q);
    }
    return YSi;
}

// [[Rcpp::export]]
void calculate_YiSit_inplace(   Eigen::Map<Eigen::MatrixXd> & YSi,
                                Eigen::Map<Eigen::MatrixXd> & QxQStore,
                                const Eigen::Map<Eigen::MatrixXd> & Y,
                                const Eigen::Map<Eigen::MatrixXd> & Si,
                                int Q) {

    //std::cout << "in: calculate_YiSit_inplace" << std::endl;

    int N = Y.cols() / Q;
    int V = Y.rows();

    //std::cout << "Size of YSi is: " << YSi.rows() << " by " << YSi.cols() << std::endl;

    for (int i = 0; i < N; i++){
        YSi.block(i*Q, 0, Q, Q).noalias() = Y.block(0, i*Q, V, Q).transpose() * Si.block(0, i*Q, V, Q);
        //QxQStore.noalias() = Y.block(0, sInd, V, Q).transpose() * Si.block(0, sInd, V, Q);
        /*
        for (int q = 0; q < Q; q++){
            for (int qprime = 0; qprime < Q; qprime++){
                YSi(sInd + qprime, q) = QxQStore(qprime, q);
            }
        }
         */

    }

}





// [[Rcpp::export]]
void calculate_AitYi_inplace(   Eigen::Map<Eigen::MatrixXd> & Si,
                                Eigen::Map<Eigen::MatrixXd> & VxQStore,
                                const Eigen::Map<Eigen::MatrixXd> & A,
                                const Eigen::Map<Eigen::MatrixXd> & Y,
                                int Q) {

  //std::cout << "in: calculate_YiSit_inplace" << std::endl;

  int N = Y.cols() / Q;
  int V = Y.rows();

  //std::cout << "Size of YSi is: " << YSi.rows() << " by " << YSi.cols() << std::endl;

  for (int i = 0; i < N; i++){
    VxQStore.noalias() = Y.block(0, i*Q, V, Q) * A.block(i*Q, 0, Q, Q);
     for (int q = 0; q < Q; q++){
       for (int v = 0; v < V; v++){
         Si(v, i*Q + q) = VxQStore(v, q);
       }
     }

  }

}



//
// [[Rcpp::export]]
void calculate_AitYi(Eigen::Map<Eigen::MatrixXd> & Si,
                     const Eigen::Map<Eigen::MatrixXd> & A,
                     const Eigen::Map<Eigen::MatrixXd> & Y,
                     int Q) {

    //std::cout << "in: calculate_AitYi" << std::endl;


    int N = A.rows() / Q;
    int V = Y.rows();
    int sInd = -1;
    int eInd = -1;

    /*
    std::cout << "size of Si =  " << Si.rows() << " by " << Si.cols() << std::endl;
    std::cout << "size of A =  " << A.rows() << " by " << A.cols() << std::endl;
    std::cout << "size of Y =  " << Y.rows() << " by " << Y.cols() << std::endl;
    std::cout << "Q =  " << Q << std::endl;
    */

    for (int i = 0; i < N; i++){
       //std::cout << "i = " << i << std::endl;
        sInd = eInd + 1;
        eInd += Q;
       //std::cout << "start: " << sInd << " end: " << eInd << std::endl;
       //std::cout << "not running!" << std::endl;
        Si.block(0, sInd, V, Q).noalias() = Y.block(0, sInd, V, Q) * A.block(sInd, 0, Q, Q);
    }

    //std::cout << "done " << std::endl;
}



//
// [[Rcpp::export]]
    void sample_u(Eigen::Map<Eigen::MatrixXd> & u,
              const Eigen::Map<Eigen::VectorXd> & stick_breaking_weights,
              const Eigen::Map<Eigen::MatrixXi> & cluster_membership) {

       //std::cout << "in: sample_u" << std::endl;


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    int V = u.rows();
    int Q = u.cols();

    for (int q = 0; q < Q; q++){
        for (int v = 0; v < V; v++){
            u(v, q) = distribution(gen) * stick_breaking_weights(cluster_membership(v, q));
        }
    }

}



void likelihood_log_sum_exp(Eigen::VectorXd& probs,
                 Eigen::VectorXd& log_probs){

    size_t resolution = probs.size();

    double max_log_prob = log_probs.maxCoeff();

    if (max_log_prob == (-std::numeric_limits<double>::infinity())){
        std::cout << "Warning - largest probability is underflowing" << std::endl;
    }

    for (int i = 0; i < resolution; i++){
        probs(i) = exp(log_probs(i) - max_log_prob);
    }

    double sumprob = probs.sum();
    for (int i = 0; i < resolution; i++){
      probs(i) =  probs(i) / sumprob;
    }

}


// [[Rcpp::export]]
void likelihood_log_sum_exp1(Eigen::Map<Eigen::VectorXd> & probs,
                            Eigen::Map<Eigen::VectorXd> & log_probs){

  size_t resolution = probs.size();

  double max_log_prob = log_probs.maxCoeff();

  if (max_log_prob == (-std::numeric_limits<double>::infinity())){
    std::cout << "Warning - largest probability is underflowing" << std::endl;
  }

  for (int i = 0; i < resolution; i++){
    probs(i) = exp(log_probs(i) - max_log_prob);
  }

  double sumprob = probs.sum();
  for (int i = 0; i < resolution; i++){
    probs(i) =  probs(i) / sumprob;
  }

}

// [[Rcpp::export]]
void likelihood_log_sum_exp2(Eigen::Map<Eigen::VectorXd> & probs,
                             Eigen::Map<Eigen::VectorXd> & log_probs){

  size_t resolution = probs.size();

  double max_log_prob = log_probs.maxCoeff();

  if (max_log_prob == (-std::numeric_limits<double>::infinity())){
    std::cout << "Warning - largest probability is underflowing" << std::endl;
  }

  for (int i = 0; i < resolution; i++){
    probs(i) = max_log_prob + log(exp(log_probs(i) - max_log_prob));
  }

  double sumprob = probs.sum();
  for (int i = 0; i < resolution; i++){
    probs(i) =  probs(i) / sumprob;
  }



}


//
// [[Rcpp::export]]
void sample_cluster_membership_indicators(Eigen::Map<Eigen::MatrixXi> & cluster_memberships,
                                          Eigen::Map<Eigen::VectorXd> & n_in_cluster,
                                          const Eigen::Map<Eigen::MatrixXd> & S0,
                                          const Eigen::Map<Eigen::MatrixXd> & u,
                                          const Eigen::Map<Eigen::VectorXd> & stick_breaking_weights,
                                          const Eigen::Map<Eigen::VectorXd> & mu_h,
                                          const Eigen::Map<Eigen::VectorXd> & sigma_sq_h,
                                          const Eigen::Map<Eigen::VectorXd> & sigma_sq_q,
                                          int H) {

   //std::cout << "in: sample_cluster_membership_indicators" << std::endl;


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double draw = 0.0;
    double sumprob = 0.0;

    int V = cluster_memberships.rows();
    int Q = cluster_memberships.cols();

    H = n_in_cluster.size();
    Eigen::VectorXd loglik = Eigen::VectorXd::Zero(H);
    Eigen::VectorXd cluster_probs = Eigen::VectorXd::Zero(H);

    double resid = 0.0;

    // Zero out cluster counts
    for (int h = 0; h < n_in_cluster.size(); h++){
        n_in_cluster(h) = 0;
    }

    for (int q = 0; q < Q; q++){
        for (int v = 0; v < V; v++){

            for (int h = 0; h < n_in_cluster.size(); h++){

                if (u(v, q) < stick_breaking_weights(h)){
                    resid = pow(S0(v, q) - mu_h(h), 2);
                    loglik(h) = -0.5 * log(sigma_sq_h(h) * sigma_sq_q(q)) - 0.5 * resid / (sigma_sq_h(h)*sigma_sq_q(q));
                } else {
                    loglik(h) = -std::numeric_limits<double>::infinity();
                }

            } // clusters (h)

            // Log sum exp
            likelihood_log_sum_exp(cluster_probs, loglik);

          /*
            if (q == 1){
              if (v == 1){
                std::cout << "normalized probabilities:" << std::endl;
                for (int h = 0; h < n_in_cluster.size(); h++){
                  if (h < 12){
                    std::cout << "S0 = " << S0(v, q) << " index: " << h << " cluster mean: " << mu_h(h) << " clustervar: " << sigma_sq_h(h) << " prob: " << cluster_probs(h) << std::endl;
                  }
                }
              }
            }
           */

            draw = distribution(gen);
            sumprob = 0.0;
            for (int h = 0; h < n_in_cluster.size(); h++){
                sumprob += cluster_probs(h);
                if (draw < sumprob){
                    cluster_memberships(v, q) = h;
                    n_in_cluster(h) += 1;
                    break;
                } else if (h == n_in_cluster.size()-1){
                    cluster_memberships(v, q) = h;
                    n_in_cluster(h) += 1;
                }


            }


        }
    }

}




// Deal with clusters that no longer have any members
// [[Rcpp::export]]
void generate_cluster_mapping_vector_inplace(Eigen::Map<Eigen::VectorXd> & mapping,
                                     const Eigen::Map<Eigen::VectorXd> & n_in_cluster){

   //std::cout << "in: generate_cluster_mapping_vector_inplace" << std::endl;


    int H = n_in_cluster.size();
    int cluster_count = 0;

    for (int h = 0; h < H; h++){
        mapping(h) = 0;
    }

    for (int h = 0; h < H; h++){
        if (n_in_cluster(h) > 0){
            mapping(cluster_count) = h;
            cluster_count++;
        }
    }


}


// [[Rcpp::export]]
void cleanup_cluster_ordering_inplace(Eigen::Map<Eigen::MatrixXi> & cluster_membership,
                                      Eigen::Map<Eigen::VectorXd> & mu_h,
                                      Eigen::Map<Eigen::VectorXd> & sigma_sq_h,
                                      Eigen::Map<Eigen::VectorXd> & n_in_cluster,
                                      Eigen::Map<Eigen::VectorXd> & mapping,
                                      int actual_cluster_count){

    int H = mu_h.size();
    int cluster_count = 0;

    int V = cluster_membership.rows();
    int Q = cluster_membership.cols();

    for (int h = 0; h < actual_cluster_count; h++){
        if (h != mapping(h)){

            // Fix cluster memberships
            for (int q = 0; q < Q; q++){
                for (int v = 0; v < V; v++){
                    if (cluster_membership(v, q) == mapping(h)){
                        cluster_membership(v, q) = h;
                    }
                }
            }

            // Fix cluster parameter ordering
            mu_h(h)       = mu_h(mapping(h));
            sigma_sq_h(h) = sigma_sq_h(mapping(h));
            n_in_cluster(h) = n_in_cluster(mapping(h));

        }
    }

    // "Zero-out" remaining cluster positions ? TODO -> might not need to
    for (int h = actual_cluster_count; h < H; h++){
        n_in_cluster(h) = 0;
    }

}






// [[Rcpp::export]]
void sample_spatial_maps_inplace(Eigen::Map<Eigen::MatrixXd> & S0,
                                 Eigen::Map<Eigen::MatrixXd> & beta,
                                 Eigen::Map<Eigen::VectorXd> & posterior_mean,
                                 Eigen::Map<Eigen::MatrixXd> & posterior_variance,
                                 Eigen::Map<Eigen::MatrixXd> & prior_precision,
                                 Eigen::Map<Eigen::VectorXd> & draw,
                                 const Eigen::Map<Eigen::MatrixXd> & AtY,
                                 const Eigen::Map<Eigen::MatrixXd> & X,
                                 const Eigen::Map<Eigen::MatrixXd> & XtX,
                                 const Eigen::Map<Eigen::VectorXd> & sigma_sq_q,
                                 const Eigen::Map<Eigen::MatrixXd> & tau_sq,
                                 const Eigen::Map<Eigen::MatrixXd> & lambda_sq,
                                 const Eigen::Map<Eigen::VectorXd> & mu_h_div_sigma_sq_h,
                                 const Eigen::Map<Eigen::VectorXd> & sigma_sq_h,
                                 const Eigen::Map<Eigen::MatrixXi> & cluster_membership,
                                 const Eigen::Map<Eigen::MatrixXd> & I){

   //std::cout << "in: sample_spatial_maps_inplace" << std::endl;


    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

   //std::cout << "dims" << std::endl;
    int Q = S0.cols();
    int N = X.rows();
    int V = S0.rows();
    int P = X.cols() - 1; // X includes a column for the intercept (S0)

    int n_fallback_trigger = 0;

    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;

    Eigen::MatrixXd posterior_precision_safe = Eigen::MatrixXd::Zero(P+1, P+1);

    // update each spatial map voxel
    for (int q = 0; q < Q; q++){

        // setup the prior precision for this component (without the lambda term)
        for (int p = 1; p < P+1; p++){
            prior_precision(p, p) = 1.0 / (tau_sq(q, p-1) * sigma_sq_q(q));
        }

        for (int v = 0; v < V; v++){

            // Finish setup for prior precision
            // first diagonal element is S0 -> comes from DPM cluster membership
            // remaining P diagonals are product of variance terms
            prior_precision(0, 0) = 1.0 / (sigma_sq_h(cluster_membership(v, q)) *  sigma_sq_q(q));

            for (int p = 0; p < P; p++){
                prior_precision(p+1, p+1) = prior_precision(p+1, p+1) / lambda_sq(v, q*P + p);
            }


            // Finish calculation of posterior variance
            // TODO make sure update in place!
            // TODO XtX / sigma outside of loop!
            posterior_variance = prior_precision + XtX / sigma_sq_q(q);


            // Calculate the posterior mean
            posterior_mean =  Eigen::VectorXd::Zero(P+1);
            for (int i = 0; i < N; i++){
                posterior_mean += AtY(v, q + Q*i) * X.row(i);
            }

            /*
            if (v == 14877 & q == 3){
                std::cout << "posterior mean before division" << std::endl;
                std::cout << posterior_mean << std::endl;
            }
             */

            posterior_mean = posterior_mean / sigma_sq_q(q);

            // Add prior mean
            posterior_mean(0) += mu_h_div_sigma_sq_h(cluster_membership(v, q)) / sigma_sq_q(q);

            // Backsolve to get final posterior mean
            posterior_precision_safe = posterior_variance;

            Eigen::FullPivLU<Eigen::Ref<Eigen::MatrixXd> > lu(posterior_variance);

            // draw standard normal random variable
            for (int p = 0; p < P+1; p++){
                draw(p) = distribution(gen);
            }

            if (lu.isInvertible()){
                posterior_mean = lu.solve(posterior_mean);
                // dress it up to MVN draw with posterior mean and variance
                draw = cod.compute(posterior_precision_safe).solve(I).llt().matrixL() * draw;
                posterior_mean += draw;
            } else {
                // This can trigger if product of global and local shrinkage parameter is too small
                // can happen by chance, but should not happen consistently
                n_fallback_trigger++ ;
                for (int p = 0; p < P+1; p++){
                    posterior_mean(p) = posterior_mean(p) / posterior_precision_safe(p, p) +
                        sqrt(1.0 / posterior_precision_safe(p, p)) * draw(p);
                }
            }

            /*
            if (v == 8839 & q == 10){
                std::cout << "posterior mean unnormalized" << std::endl;
                std::cout << posterior_mean << std::endl;
            }

            if (v == 8839 & q == 10){
                std::cout << "uninverted:" << std::endl;
                std::cout << posterior_precision_safe << std::endl;
            }
             */

            /*
            if (v == 14877 & q == 3){
                std::cout << "prior precision at setp 5 should be same" << std::endl;
                std::cout << posterior_variance << std::endl;
                std::cout << "versus" << std::endl;
            }
            */

            /*
            if (v == 8839 & q == 10){
                std::cout << "posterior mean" << std::endl;
                std::cout << posterior_mean << std::endl;
            }
             */

            /*
            if (v == 14877 & q == 3){
                std::cout << "Uninverted matrix" << std::endl;
                std::cout << posterior_variance << std::endl;
                //std::cout << posterior_variance.inverse() << std::endl;
            }


            if  (v == 8839 & q == 10){
                std::cout << "Checking invertable:" << std::endl;
                std::cout << cod.compute(posterior_variance).isInvertible() << std::endl;
                std::cout << "Inverse" << std::endl;
                std::cout << cod.compute(posterior_variance).solve(I) << std::endl;
                //std::cout << posterior_variance.inverse() << std::endl;
            }
             */






             // Update S0 and beta elements

             S0(v, q) = posterior_mean(0);
             for (int p = 1; p < P+1; p++){
                beta(v, q*P + p-1) = posterior_mean(p);
             }



            // Revert change to prior precision
            for (int p = 0; p < P; p++){
                prior_precision(p+1, p+1) = prior_precision(p+1, p+1) * lambda_sq(v, q*P + p);
            }


        } // v


    } // q

    //std::cout << "Triggered Fallback sampler " << n_fallback_trigger << " times" << std::endl;


   //std::cout << "done in spatial maps" << std::endl;

}














// [[Rcpp::export]]
void update_spatial_map_posterior_tracking_inplace(Eigen::Map<Eigen::MatrixXd> & S0_pm,
                                 Eigen::Map<Eigen::MatrixXd> & beta_pm,
                                 Eigen::Map<Eigen::MatrixXd> & Beta_ngt0,
                                 Eigen::Map<Eigen::MatrixXd> & lambda_sq_posterior_mean,
                                 const Eigen::Map<Eigen::MatrixXd> & S0,
                                 const Eigen::Map<Eigen::MatrixXd> & beta,
                                 const Eigen::Map<Eigen::MatrixXd> & lambda_sq,
                                 double n_mcmc){

    int Q = S0.cols();
    int V = S0.rows();
    int P = beta.cols();

    for (int q = 0; q < Q; q++){
        for (int v = 0; v < V; v++){
            S0_pm(v, q) += (1 / n_mcmc) * S0(v, q);
        }
    }

    // note that P here is all columns of Beta, so is really QP
    for (int p = 0; p < P; p++){
        for (int v = 0; v < V; v++){
            beta_pm(v, p) += (1 / n_mcmc) * beta(v, p);
            Beta_ngt0(v, p) += (beta(v, p) > 0.0);
            lambda_sq_posterior_mean(v, p) += (1 / n_mcmc) * lambda_sq(v, p);
        }
    }


}










// [[Rcpp::export]]
void calculate_DPM_hyperparameter_suffstats_inplace(Eigen::Map<Eigen::MatrixXd> & sum_Ssq_qh,
                                                    Eigen::Map<Eigen::MatrixXd> & Sbar_qh,
                                                    Eigen::Map<Eigen::MatrixXd> & n_in_clust_qh,
                                                    const Eigen::Map<Eigen::MatrixXd> & S0,
                                                    const Eigen::Map<Eigen::MatrixXi> & cluster_membership){

   //std::cout << "in: calculate_DPM_hyperparameter_suffstats_inplace" << std::endl;


    int Q = S0.cols();
    int V = S0.rows();
    int maxH = sum_Ssq_qh.rows();

    // Zero out current contents of sufficient statistics
    for (int q = 0; q < Q; q++){
        for (int h = 0; h < maxH; h++){
            sum_Ssq_qh(h, q) = 0.0;
            Sbar_qh(h, q) = 0.0;
            n_in_clust_qh(h, q) = 0.0;
        }
    }

    // Evaluate sufficient statistics
    for (int q = 0; q < Q; q++){
        for (int v = 0; v < V; v++){
            sum_Ssq_qh(cluster_membership(v, q), q) += pow(S0(v, q), 2);
            Sbar_qh(cluster_membership(v, q), q)    += S0(v, q);
            n_in_clust_qh(cluster_membership(v, q), q)    += 1;
        }

        // Now obtain the S - Sbar squared term

        // turn sum into mean for Sbar_qh
        // - this gets done after total is calculated to make sure
        // each worker processes' data is weighted correctly
        //for (int h = 0; h < maxH; h++){
        //    if (n_in_clust_qh(q, h) > 0){
        //        Sbar_qh(h, q) = Sbar_qh(h, q) / n_in_clust_qh(h, q);
        //    }
        //}
    }

}






// [[Rcpp::export]]
void sample_DPM_hyperparameters_inplace(Eigen::Map<Eigen::VectorXd> & mu_h,
                                        Eigen::Map<Eigen::VectorXd> & sigma_sq_h,
                                        const Eigen::Map<Eigen::MatrixXd> & sum_Ssq_qh,
                                        const Eigen::Map<Eigen::MatrixXd> & Sbar_qh,
                                        const Eigen::Map<Eigen::MatrixXd> & DPM_hyper_sse_unscaled,
                                        const Eigen::Map<Eigen::MatrixXd> & n_in_cluster_qh,
                                        const Eigen::Map<Eigen::VectorXd> & n_in_cluster,
                                        const Eigen::Map<Eigen::VectorXd> & sigma_sq_q,
                                        double prior_shape,
                                        double prior_scale,
                                        int maxH){


   //std::cout << "in: sample_DPM_hyperparameters_inplace" << std::endl;


    std::random_device rd;
    std::mt19937 gen(rd());
    std::gamma_distribution<double> gamma(1.0, 1.0); // shape, scale param
    std::normal_distribution<double> norm(0.0, 1.0);

    int Q    = sum_Ssq_qh.cols();

    // intermediate quantities
    double posterior_shape = 0.0;
    double posterior_scale = 0.0;

    double T1 = 0.0;
    double T2 = 0.0;
    double T3 = 0.0;

    for (int h = 0; h < maxH; h++){

        if ( n_in_cluster(h)  > 0){

            T1 = 0.0;
            T2  = 0.0;
            T3   = 0.0;

            posterior_shape = prior_shape + n_in_cluster(h) / 2.0 + 2.0;

            for (int q = 0; q < Q; q++){
                //T1 += sum_Ssq_qh(h, q) / sigma_sq_q(q);
                //T2 += 2 * n_in_cluster_qh(h, q)*Sbar_qh(h, q) / sigma_sq_q(q);
                //T3 += n_in_cluster_qh(h, q) / sigma_sq_q(q);
                T1 += sum_Ssq_qh(h, q) / sigma_sq_q(q);
                T2 += 2 * n_in_cluster_qh(h, q)*Sbar_qh(h, q) / sigma_sq_q(q);
                T3 += n_in_cluster_qh(h, q) / sigma_sq_q(q);
            }

            posterior_scale = (T1 + 2.0 * prior_scale - pow(T2, 2) / (4 * (T3 + 2)) );

            // Sample sigma_sq_h - if X ~ InverseGamma(a, b), then
            // 1 / X ~ Gamma(a, 1/b)
            gamma = std::gamma_distribution<double>(posterior_shape, 2.0 / posterior_scale);
            //gamma = std::gamma_distribution<double>(posterior_shape,posterior_scale);
            sigma_sq_h(h) = 1.0 / gamma(gen);

            /*
            std::cout << "Updating cluster " << h << std::endl;
            std::cout << "T1: " << T1 << std::endl;
            std::cout << "T2: " << T2 << std::endl;
            std::cout << "T3: " << T3 << std::endl;
             */

            mu_h(h) = norm(gen) * sqrt(sigma_sq_h(h) / ((T3 + 2.0) ) ) + (T2 / (2*(T3+2)));
        }


    }




}










// [[Rcpp::export]]
void evaluate_error_inplace(Eigen::Map<Eigen::MatrixXd> & Ei,
                            Eigen::Map<Eigen::MatrixXd> & Si,
                            const Eigen::Map<Eigen::MatrixXd> & AtY,
                            const Eigen::Map<Eigen::MatrixXd> & S0,
                            const Eigen::Map<Eigen::MatrixXd> & beta,
                            const Eigen::Map<Eigen::MatrixXd> & X){

    int V = S0.rows();
    int Q = S0.cols();
    int N = X.rows();
    int P = X.cols() - 1;

    for (int i = 0; i < N; i++){
        for (int q = 0; q < Q; q++){
            for (int v = 0; v < V; v++){
                Si(v, Q*i + q) = S0(v, q);
                // TODO replace this with vector product
                for (int p = 0; p < P; p++){
                    Si(v, Q*i + q) += X(i, p+1) * beta(v, q*P + p);
                }
                Ei(v, Q*i + q) = AtY(v, Q*i + q) - Si(v, Q*i + q);

            }
        }
    }

}











// [[Rcpp::export]]
void evaluate_sse_inplace(Eigen::Map<Eigen::VectorXd> & sse,
                          const Eigen::Map<Eigen::MatrixXd> & Ei){


   //std::cout << "in: evaluate_sse_inplace" << std::endl;

    int V = Ei.rows();
    int Q = sse.size();
    int N = Ei.cols() / Q;

    for (int q  = 0; q < Q; q++){
        sse(q) = 0.0;
        for (int v = 0; v < V; v++){
            for (int i = 0; i < N; i++){
                sse(q) += pow(Ei(v, Q*i + q), 2);
            } // i
        } // v
    } // q


}

// [[Rcpp::export]]
void evaluate_ssY_inplace(Eigen::Map<Eigen::VectorXd> & sse,
                          const Eigen::Map<Eigen::MatrixXd> & AtY){


  //std::cout << "in: evaluate_sse_inplace" << std::endl;

  int V = AtY.rows();
  int Q = sse.size();
  int N = AtY.cols() / Q;

  std::cout << Q << " " << V << " " << N << "" << std::endl;

  for (int q  = 0; q < Q; q++){
    sse(q) = 0.0;
    for (int v = 0; v < V; v++){
      for (int i = 0; i < N; i++){
        sse(q) += pow(AtY(v, Q*i + q), 2);
      } // i
    } // v
  } // q


}





// [[Rcpp::export]]
void evaluate_reg_scale_inplace(Eigen::Map<Eigen::VectorXd> & reg_scale,
                          const Eigen::Map<Eigen::MatrixXd> & S0,
                          const Eigen::Map<Eigen::MatrixXd> & beta,
                          const Eigen::Map<Eigen::VectorXd> & mu_h,
                          const Eigen::Map<Eigen::VectorXd> & sigma_sq_h,
                          const Eigen::Map<Eigen::MatrixXd> & tau_sq,
                          const Eigen::Map<Eigen::MatrixXi> & cluster_membership,
                          const Eigen::Map<Eigen::MatrixXd> & lambda_sq){

   //std::cout << "in: evaluate_reg_scale_inplace" << std::endl;

    int V = S0.rows();
    int Q = S0.cols();
    int P = beta.cols() / Q;

    for (int q  = 0; q < Q; q++){
        reg_scale(q) = 0.0;
        for (int v = 0; v < V; v++){
            reg_scale(q) += pow(S0(v, q) - mu_h(cluster_membership(v, q)), 2) /
                sigma_sq_h(cluster_membership(v, q));
            for (int p = 0; p < P; p++){
                reg_scale(q) += pow(beta(v, q*P + p), 2) / (tau_sq(q, p) * lambda_sq(v, q*P + p));
            }
        } // v
    } // q

}



/*
void evaluate_sigma_scale_inplace(Eigen::Map<Eigen::VectorXd> & reg_scale,
                                const Eigen::Map<Eigen::MatrixXd> & diag_projection,
                                const Eigen::Map<Eigen::VectorXd> & mu_h,
                                const Eigen::Map<Eigen::VectorXd> & sigma_sq_h,
                                const Eigen::Map<Eigen::MatrixXd> & tau_sq,
                                const Eigen::Map<Eigen::MatrixXi> & cluster_membership,
                                const Eigen::Map<Eigen::MatrixXd> & lambda_sq){

    //std::cout << "in: evaluate_reg_scale_inplace" << std::endl;

    int V     = diag_projection.rows();
    int Q     = tau_sq.rows();
    // Pstar is P + 1
    int Pstar = diag_projection.cols() / Q;

    for (int q  = 0; q < Q; q++){
        reg_scale(q) = 0.0;
        for (int v = 0; v < V; v++){
            //reg_scale(q) += pow(S0(v, q) - mu_h(cluster_membership(v, q)), 2) /
            //    sigma_sq_h(cluster_membership(v, q));
            reg_scale(q) += diag_projection(v, q*P)
            for (int p = 0; p < P; p++){
                //reg_scale(q) += pow(beta(v, q*P + p), 2) / (tau_sq(q, p) * lambda_sq(v, q*P + p));
            }
        } // v
    } // q

}
*/












// [[Rcpp::export]]
void sample_local_shrinkage_parameters_inplace(Eigen::Map<Eigen::MatrixXd> & lambda_sq,
                                               Eigen::Map<Eigen::MatrixXd> & nu,
                                               Eigen::Map<Eigen::MatrixXd> & cauchy_mixing,
                                               Eigen::Map<Eigen::MatrixXd> & cauchy_mixing_prior,
                                               const Eigen::Map<Eigen::MatrixXd> & beta,
                                               const Eigen::Map<Eigen::MatrixXd> & tau_sq,
                                               const Eigen::Map<Eigen::VectorXd> & sigma_sq_q){

   //std::cout << "in: sample_local_shrinkage_parameters_inplace" << std::endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::gamma_distribution<double> gamma(1.0, 1.0); // shape, scale param

    int V = beta.rows();
    int Q = sigma_sq_q.size();
    int P = lambda_sq.cols() / Q;

    // temporary
    double denom = 0.0;
    double scale = 0.0;

    for (int q = 0; q < Q; q++){
        for (int p = 0; p < P; p++){

            denom = 2.0 * tau_sq(q, p) * sigma_sq_q(q);
            for (int v = 0; v < V; v++){

                // if X ~ InverseGamma(a, b), with rate b, then
                // then  1 / X ~ Gamma(a, b)
                // c++ gamma function default second arg is the scale
                // Paper defines inverse gamma in terms of a scale parameter beta
                // inverse gamma with scale = 1/beta is equivalent to 1 / gamma(alpha, scale = beta)


                // this is actually the rate, not the scale
                scale = (1.0 / nu(v, P*q + p)) +
                    pow(beta(v, P*q + p), 2) / denom;
                //gamma = std::gamma_distribution<double>(1.0, 1.0 / scale);
                gamma = std::gamma_distribution<double>(1.0, 1.0 / scale);
                lambda_sq(v, P*q + p) = 1.0 / gamma(gen);

                scale = (1.0 / cauchy_mixing(v, P*q + p)) + (1.0 / lambda_sq(v, P*q + p));
                gamma = std::gamma_distribution<double>(1.0, 1.0 / scale);
                nu(v, P*q + p) = 1.0 / gamma(gen);

                scale = (1.0 / cauchy_mixing_prior(v, P*q + p)) + (1.0 / nu(v, P*q + p));
                gamma = std::gamma_distribution<double>(1.0, 1.0 / scale);
                cauchy_mixing(v, P*q + p) = 1.0 / gamma(gen);

                scale = 1.0 + (1.0 / cauchy_mixing(v, P*q + p));
                gamma = std::gamma_distribution<double>(1.0, 1.0 / scale);
                cauchy_mixing_prior(v, P*q + p) = 1.0 / gamma(gen);
            }
        }
    }


}






// [[Rcpp::export]]
void calculate_tausum_inplace(Eigen::Map<Eigen::MatrixXd> & tau_sum,
                              const Eigen::Map<Eigen::MatrixXd> & beta,
                              const Eigen::Map<Eigen::MatrixXd> & lambda_sq,
                              const Eigen::Map<Eigen::VectorXd> & sigma_sq_q){

   //std::cout << "in: calculate_tausum_inplace" << std::endl;

    int V = beta.rows();
    int Q = sigma_sq_q.size();
    int P = lambda_sq.cols() / Q;

    for (int q = 0; q < Q; q++){
        for (int p = 0; p < P; p++){
            tau_sum(q, p) = 0.0;
            for (int v = 0; v < V; v++){
                tau_sum(q, p) += (pow(beta(v, q*P + p), 2) / (2.0 * lambda_sq(v, q*P + p) * sigma_sq_q(q)));
            }
        }
    }


}




// [[Rcpp::export]]
void update_vector_in_place(Eigen::Map<Eigen::VectorXd> & variable,
                     const Eigen::Map<Eigen::VectorXd> & new_level){

  if (variable.size() != new_level.size()){
    std::cout << "ERROR - SIZES DO NOT MATCH" << std::endl;
    return;
  }

  for (int i = 0; i < variable.size(); i++){
    variable(i) = new_level(i);
  }

}

// [[Rcpp::export]]
void update_matrix_in_place(Eigen::Map<Eigen::MatrixXd> & variable,
                            const Eigen::Map<Eigen::MatrixXd> & new_level){

  if (variable.rows() != new_level.rows()){
    std::cout << "ERROR - ROW DIMENSIONS DO NOT MATCH" << std::endl;
    return;
  }

  if (variable.cols() != new_level.cols()){
    std::cout << "ERROR - COLUMN DIMENSIONS DO NOT MATCH" << std::endl;
    return;
  }

  for (int j = 0; j < variable.cols(); j++){
    for (int i = 0; i < variable.rows(); i++){
      variable(i, j) = new_level(i, j);
    }
  }

}


// NOTE indices start at 1 - have to offset
// [[Rcpp::export]]
void update_matrix_in_place_rowspec(Eigen::Map<Eigen::MatrixXd> & variable,
                            const Eigen::Map<Eigen::MatrixXd> & new_level,
                            const Eigen::Map<Eigen::VectorXi> & indices){

  if (variable.rows() != indices.size()){
    std::cout << "ERROR - number of indices provided in indices variable does not match number of rows in variable to update" << std::endl;
    return;
  }

  if (variable.cols() != new_level.cols()){
    std::cout << "ERROR - COLUMN DIMENSIONS DO NOT MATCH" << std::endl;
    return;
  }

  for (int j = 0; j < variable.cols(); j++){
    for (int i = 0; i < variable.rows(); i++){
      variable(i, j) = new_level(indices(i) - 1, j);
    }
  }

}


