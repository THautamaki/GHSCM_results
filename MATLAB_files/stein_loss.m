% Author: Tuomas Hautamäki, University of Oulu
function loss = stein_loss(true_matrix, est_matrix)
  % This function calculates Stein's loss between true and estimated positive-definite matrices,
  % which equals two times the Kullback-Leibler divergence, and is defined as follows:
  %   tr(M_est*M_true^-1) - log|M_est*M_true^-1| - p,
  % where * is matrix product, tr() is trace of a matrix, || is determinant of the matrix and p is
  % dimension of the matrix.
  %
  % Inputs
  %   true_matrix The true p by p covariance or precision matrix.
  %   est_matrix  The estimate of the p by p covariance or precision matrix.
  %
  % Outputs
  %   loss        Scalar, which is calculated Stein's loss between two matrices.
  %
  % References
  %   James, W. and Stein, C. (1961). Estimation with quadratic loss. In Proceedings of the Fourth
  %   Berkeley Symposium on Mathematical Statistics and Probability, volume 1, pages 361--379.
  %   University of California Press.
  p = size(true_matrix, 1);
  P = est_matrix / true_matrix;
  loss = trace(P) - log(det(P)) - p;
end