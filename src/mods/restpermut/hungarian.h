#ifndef __HUNGARIAN_H
#define __HUNGARIAN_H

#include <valarray>

// Description of this solution to the linear ssignment problem
//
// http://www.public.iastate.edu/~ddoty/HungarianAlgorithm.html
//


// Hungarian Algorithm for Linear assignement - the original 
// algorithm due to Kuhn and Munkres. Its pretty slow ( O(n^4) )
// 
// Kuhn, H. W.: The Hungarian method for the assignment problem. 
// Naval Research Logistics Quarterly 2, 83-97 (1955).
int LinearAssignmentHungarian(
	std::valarray<int> &costm, 
	std::vector<int> &row2col, 
	std::vector<int> &col2row
);


/// A much faster algorithm for linear assignement: 
///  Reference: Jonker, R. and Volgenant, A., 1987. A shortest augmenting
///  path algorithm for dense and sparse linear assignment problems. 
///  Computing 38, pp. 325–340
int LinearAssignmentJVC(
	std::valarray<int> &costm,
	std::vector<int> &row2col,
	std::vector<int> &col2row
);

#endif



