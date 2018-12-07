#' @title Solve $l_0$-regularized least squres with Single Best Replacement (SBR) algorithm
#'
#' @description This is the main function for solving $l_0$-regularized least squres for a given
#' penalty parameter, based on the Single Best Replacement (SBR) algorithm originally proposed by Soussen et al (2011)
#'
#' @param X A n by p numeric design matrix
#' @param y A n vector of numeric response
#' @param lambda A non-negative numeric penalty parameter
#' @export
#' @importFrom stats dnorm pnorm
#' @importFrom utils capture.output modifyList
#'
l0reg = function (X, y, lambda, control = list()) {
  A <- X
  z <- y
  lambda = 2 * lambda

	M = ncol(A) # size of vector x

	# Optional argument
	control.default = list(
	VERB = 0,
	qinit = c(),
	explore = "all",
	K = M
	)
	control = modifyList(control.default, control)
	VERB = control$VERB
	qinit = control$qinit
	explore = control$explore
	K = control$K
	K = min(K, min(dim(A)))

	# WE ASSUME THAT THE GRAM MATRIX A'A CAN BE COMPUTED AND STORED ONCE FOR ALL
	AA = t(A) %*% A
	Az = t(A) %*% z # < h_i,z>
	h2 = diag(AA) # ||h_i||^2
	z2 = sum(z^2) # ||z||^2
	Jit = rep(NA, length = K + 1)

	# zero-valued initial solution
	q = rep(0, length = M) # sparse vector of size Mx1, initialized to 0
	vecm1 = c() # support = active indices
	lvecm1 = 0       # size of support
	lvecm0 = M       # complement
#	R = c()           # Cholesky factor of the Gram matrix:
                  # R is the upper triangular matrix such that
                  #    AA(active_columns,active_columns) = R'*R
                  # With the notations of [Soussen et al, 2011], R is the
                  # transpose of L_Q (see Appendix C of [Soussen et al, 2011])
	crit = z2        # Jinit = J(0) = ||z||^2
	FLAG = 1         # Termination of SBR if FLAG = 0 (normal) or -1 (anomaly)
	iter = 0
	Jit[1] = crit
  sol_cour = 0    # R'\Az(vecm1);

# USER INITIALIZATION
if (!is.null(qinit)) {
    q[qinit] = 1
    vecm1 = qinit
    lvecm1 = length(vecm1)
    lvecm0 = M - lvecm1
    if (lvecm1) {
        R = chol(AA[vecm1, vecm1])   # factorisation AA = R'R
        sol_cour = solve(t(R), Az[vecm1])
        crit = z2 - t(sol_cour) %*% sol_cour + lambda * lvecm1
        Jit[lvecm1 + 1] = crit
    }
}

if (VERB) {cat("0:", "\t", "\t", crit, "\n")}

# EXPLORATION STRATEGY: ALL COLUMNS OR NOT

# all the columns are explored, starting by
                             # insertion tests

if (explore == 'all') {
	while ((FLAG == 1) & (lvecm1 < K)) {
		# The current value of the cost function J(x) is already stored in crit
        dcritopt = 0;       # J_new - J_old, dcritopt is always the lowest obj value
                            # in current loop
#        mopt = [];          # position of the new column to insert
#        indopt = [];        # location of the index to remove in vecm1
        mouv = 0;           # 1 if insertion, -1 if removal
        vecm0 = which(q == 0); # updated at the beginning of the loop contrary
                            # to vecm1

        #--------------------------------------------------------------
        #-------- I. INSERTION/REMOVAL TESTS     ----------------------
        #--------------------------------------------------------------
        # INSERTION TRIALS h_m (m-th column) for all m
# non empty support
	if (lvecm1) {
        	tmp_vect = solve(t(R), AA[vecm1, vecm0, drop = FALSE])  # size lvecm1 x lvecm0
            F22_vect = h2[vecm0] - colSums(tmp_vect^2);
            dcritnew_vect = lambda - (t(tmp_vect) %*% sol_cour - Az[vecm0])^2 / F22_vect;
        }        else {
        	dcritnew_vect = lambda - Az^2 / h2;
        }   # empty support

	val_min = min(dcritnew_vect);
	ind = which.min(dcritnew_vect);

        if (val_min < dcritopt) {
        	            mouv = 1;
            mopt = vecm0[ind];   # the best insertion
            dcritopt = val_min;  # deltaJ(insert m) = the least deltaJ
        }
        # END INSERTION TESTS

        # In some cases, it is not worth testing removals:
        if ((lambda > 0) & (dcritopt > -lambda) & (lvecm1 > 2)) {
            # REMOVAL TESTS
            # Inversion of R matrix
            R_inv = solve(R);
            # Computation of current amplitudes
            x_loc = R_inv %*% sol_cour;
            dcritnew_vect = x_loc^2 / rowSums(R_inv^2) - lambda;
            	val_min = min(dcritnew_vect);
	ind = which.min(dcritnew_vect);

            if(val_min < dcritopt) {
                mouv = -1;
                indopt = ind;       # index in VECM1 (not in q) of best removal
                dcritopt = val_min;
            }
            # END REMOVAL TESTS
        }

        #--------------------------------------------------------------
        #-------- II. DO INSERTION, REMOVAL OR NOTHING    -------------
        #--------------------------------------------------------------
        if (dcritopt < 0) {       # Modification of the support
            iter = iter + 1;
            crit = crit + dcritopt;

            if ((mouv == 1) & (lvecm1 > 0)) {  # instability test
                tmp = solve(t(R), AA[vecm1, mopt, drop = FALSE]);
                FLAG =  1 - 2 * (h2[mopt] < sum(tmp^2));
            }

            if ((crit < -1e-16) | (FLAG == -1)) {
                cat("ABORT - NUMERICAL INSTABILITY", "\t", "FLAG =", FLAG)
#                disp(sprintf('*** ABORT - NUMERICAL INSTABILITY, FLAG = #d ***\n',FLAG));
                FLAG = -1;
                # Amplitude computation
                x = rep(0, M);
                x[vecm1] = solve(R, sol_cour);
                stop;
            }

            if (mouv == 1) {  # INSERTION
                # Update of R
                if (lvecm1 > 0) {
                    R = cbind(R, tmp);
                    R = rbind(R, c(rep(0, lvecm1), sqrt(h2[mopt] - sum(tmp^2))));
                } else {
                    R = sqrt(h2[mopt]);
                }
                vecm1 = c(vecm1, mopt); # insertion of mopt at the last position
                lvecm0 = lvecm0 - 1;
                lvecm1 = lvecm1 + 1;
                #  note: vecm0 is updated in the beginning of the loop
                if (VERB) {cat(iter, ":", "\t", "+", mopt, "\t", crit, "\n")}
            } else {          # REMOVAL
                # Update of R
              if (indopt < nrow(R)) {
                FF = R[(indopt + 1) : nrow(R), (indopt + 1) : ncol(R)];
                E = R[indopt, (indopt + 1) : ncol(R)];
                R = R[-indopt, ];
                R = R[, -indopt];
                R[indopt : nrow(R), indopt : ncol(R)] = chol(t(FF) %*% FF + E %*% t(E));
              } else {
                R = R[-indopt, ];
                R = R[, -indopt];
              }
                mopt = vecm1[indopt];
                vecm1 = vecm1[-indopt];
                # vecm0 is updated in the beginning of the loop
                lvecm0 = lvecm0 + 1;
                lvecm1 = lvecm1 - 1;
                if (VERB) {cat(iter, ":", "\t", "-", mopt, "\t", crit, "\n")}
            }
            q[mopt] = 1 - q[mopt];   # 1-> 0 or 0 -> 1
            sol_cour = solve(t(R), Az[vecm1]);
            Jit[lvecm1 + 1] = crit;
        } else {
            FLAG = 0;   # Stop of SBR (no more decrease)
        }
	}
}	else if (explore == 'remove') { # the removals are first explored.
                                    # If one of them yields a decrease of
                                    # the cost, insertions are not tested
		    while ((FLAG == 1) & (lvecm1 < K)) {
        # The current value of the cost function J(x) is already stored in crit
        dcritopt = 0;       # J_new - J_old
#        mopt = [];          # position of the new column to insert
#        indopt = [];        # location of the index to remove in vecm1
        mouv = 0;           # 1 if insertion, -1 if removal
        vecm0 = which(q == 0); # updated at the beginning of the loop contrary
                            # to vecm1

        #--------------------------------------------------------------
        #-------- I. INSERTION/REMOVAL TESTS     ----------------------
        #--------------------------------------------------------------
        if ((lambda > 0) & (lvecm1 > 2)) {
            # REMOVAL TESTS
            # Inversion of R matrix
            R_inv = solve(R);
            # Computation of current amplitudes
            x_loc = R_inv %*% sol_cour;
            dcritnew_vect = x_loc^2 / rowSums(R_inv^2) - lambda;

            val_min = min(dcritnew_vect);
			ind = which.min(dcritnew_vect);

            if (val_min < dcritopt) {
                mouv = -1;
                indopt = ind;  # index in VECM1 (not in q) of best removal
                dcritopt = val_min;
            }
            # END REMOVAL TESTS
        }

        # If mouv=-1 the best removal is selected, else insertions are explored
        if (mouv == 0) {
            # INSERTION TRIALS h_m (m-th column) for all m
            if (lvecm1) {           # support of at least 2 elements
                tmp_vect = solve(t(R), AA[vecm1, vecm0, drop = FALSE]);  # size lvecm1 x lvecm0
                F22_vect = h2[vecm0] - colSums(tmp_vect^2);
                dcritnew_vect = lambda - (t(tmp_vect) %*% sol_cour - Az[vecm0])^2 / F22_vect;
            } else {   # empty support
                dcritnew_vect = lambda - Az^2 / h2;
            }

                       val_min = min(dcritnew_vect);
			ind = which.min(dcritnew_vect);

            if(val_min < dcritopt) {
                mouv = 1;
                mopt = vecm0[ind];  # the best insertion
                dcritopt = val_min; # deltaJ(insert m) = the least deltaJ
            }
            # END INSERTION TESTS
        }

        #--------------------------------------------------------------
        #-------- II. DO INSERTION, REMOVAL OR NOTHING    -------------
        #--------------------------------------------------------------
        if (dcritopt < 0) {       # Modification of the support
            iter = iter + 1;
            crit = crit + dcritopt;
            if ((mouv == 1) & (lvecm1 > 0)) {  # instability test
                tmp = solve(t(R), AA[vecm1, mopt, drop = FALSE]);
                FLAG = 1 - 2 * (h2[mopt] < sum(tmp^2));
            }

            if ((crit < -1e-16) | (FLAG == -1)) {
                cat("ABORT - NUMERICAL INSTABILITY", "\t", "FLAG =", FLAG);
                FLAG = -1;
                # Amplitude computation
                x = rep(0, M);
                x[vecm1] = solve(R, sol_cour);
                stop;
            }

            if (mouv == 1) {  # INSERTION
                # Update of R
                if (lvecm1 > 0) {
                    R = cbind(R, tmp);
                    R = rbind(R, c(rep(0, lvecm1), sqrt(h2[mopt] - sum(tmp^2))));
                } else {
                    R = sqrt(h2[mopt]);
                }
                vecm1 = c(vecm1, mopt); # insertion of mopt at the last position
                lvecm0 = lvecm0 - 1;
                lvecm1 = lvecm1 + 1;
                #  note: vecm0 is updated in the beginning of the loop
                if (VERB) {cat(iter, ":", "\t", "+", mopt, "\t", crit, "\n")}
            } else {          # REMOVAL
                # Update of R
              if (indopt < nrow(R)) {
                FF = R[(indopt + 1) : nrow(R), (indopt + 1) : ncol(R)];
                E = R[indopt, (indopt + 1) : ncol(R)];
                R = R[-indopt, ];
                R = R[, -indopt];
                R[indopt : nrow(R), indopt : ncol(R)] = chol(t(FF) %*% FF + E %*% t(E));
              } else {
                R = R[-indopt, ];
                R = R[, -indopt];
              }
              mopt = vecm1[indopt];
                vecm1 = vecm1[-indopt];
                # vecm0 is updated in the beginning of the loop
                lvecm0 = lvecm0 + 1;
                lvecm1 = lvecm1 - 1;
                if (VERB) {cat(iter, ":", "\t", "-", mopt, "\t", crit, "\n")}
            }
            q[mopt] = 1 - q[mopt];   # 1-> 0 or 0 -> 1
            sol_cour = solve(t(R), Az[vecm1]);
            Jit[lvecm1 + 1] = crit;
        } else {
            FLAG = 0;   # Stop of SBR (no more decrease)
        }
    } # WHILE
}

# Final amplitude computation
#
x = rep(0, M);
if (!all(sol_cour == 0)) {
  x[vecm1] = solve(R, sol_cour);
}
x.sparse = cbind(which(x != 0), x[which(x != 0)])
colnames(x.sparse) = c("position", "value")

obj.nonzero <- cbind(
  num.nonzero = which(!is.na(Jit)) - 1,
  opt.obj.value = Jit[!is.na(Jit)]
)

result <- list(nonzero.coef = x.sparse, obj.val = obj.nonzero, status = FLAG,
               coef = x, X = X, y = y, lambda = lambda / 2, explore = explore, K = K
               )
class(result) <- 'l0reg'

return(result)
}

#' @export

summary.l0reg <- function (result, ...) {
  result[1 : 3]
}

#' @export

print.l0reg <- function (result, ...) {
  print(summary.l0reg(result))
}
