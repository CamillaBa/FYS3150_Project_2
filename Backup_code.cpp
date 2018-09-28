//
//
//// update similarity transformation matrix S
//S = setzero(S, n, n);
//for (int i = 0; i < n; i++) { S[i][i] = 1; }
//S[l][l] = c; S[k][k] = c; S[k][l] = -s; S[l][k] = s;
//condition = false;
//





//sp_mat A(n, n); sp_mat A_nodiag(n, n); sp_mat S(n, n);
//A.diag() += 2;
//A.diag(1) += -1;
//A.diag(-1) += -1;


//
//int k;
//int l;
//double tau;
//uvec sub;
//bool condition = true;
//double t; //tan(theta)
//double c; //cos(theta)
//double s; //sin(theta)
//double root;

//while (condition) {
//	// check if max(a_ij) > epsilon
//	A_nodiag = A;
//	A_nodiag.diag() = zeros<vec>(n); //remove diagonal from A
//	A_nodiag = square(A_nodiag);  //square matrix entries (in order to seek the largest)
//	condition = (A_nodiag.max() > epsilon); //update condition (loop stops if condition = false)

//	// calculate values for tau, tan(theta), cos(theta), and sin(theta)
//	k = A_nodiag.index_max(); // obtain linear index for max value
//	sub = ind2sub(size(A_nodiag), k); // convert linear index to sub = (k,l)
//	k = sub(0); l = sub(1); // extract k and l from sub
//	tau = (A(l, l) - A(k, k)) / (2.0 * A(k, l));
//	root = sqrt(1+tau*tau);
//	t = fmin(-tau-root,-tau+root);
//	c = 1.0 /sqrt(1+t*t);
//	s = t * c;

//	// update similarity transformation matrix S
//	S = zeros<sp_mat>(n, n);
//	S.diag() += 1;
//	S(l, l) = c; S(k, k) = c; S(k, l) = -s; S(l, k) = s;

//	// update A
//	A = S.t()*A*S;
//}

//for (int i = 0; i < n; i++) {
//	for (int j = 0; j < n; j++) {
//		std::cout << std::setw(8)<< A(i, j) << " ";
//	}
//	std::cout << std::endl;
//}
