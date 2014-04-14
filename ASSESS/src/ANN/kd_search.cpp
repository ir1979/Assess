//----------------------------------------------------------------------
// File:			kd_search.cpp
// Programmer:		Sunil Arya and David Mount
// Description:		Standard kd-tree search
// Last modified:	01/04/05 (Version 1.0)
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------
// History:
//	Revision 0.1  03/04/98
//		Initial release
//	Revision 1.0  04/01/05
//		Changed names LO, HI to ANN_LO, ANN_HI
//----------------------------------------------------------------------

#include "kd_search.h"					// kd-search declarations
#include <iostream>

//----------------------------------------------------------------------
//	Approximate nearest neighbor searching by kd-tree search
//		The kd-tree is searched for an approximate nearest neighbor.
//		The point is returned through one of the arguments, and the
//		distance returned is the squared distance to this point.
//
//		The method used for searching the kd-tree is an approximate
//		adaptation of the search algorithm described by Friedman,
//		Bentley, and Finkel, ``An algorithm for finding best matches
//		in logarithmic expected time,'' ACM Transactions on Mathematical
//		Software, 3(3):209-226, 1977).
//
//		The algorithm operates recursively.  When first encountering a
//		node of the kd-tree we first visit the child which is closest to
//		the query point.  On return, we decide whether we want to visit
//		the other child.  If the box containing the other child exceeds
//		1/(1+eps) times the current best distance, then we skip it (since
//		any point found in this child cannot be closer to the query point
//		by more than this factor.)  Otherwise, we visit it recursively.
//		The distance between a box and the query point is computed exactly
//		(not approximated as is often done in kd-tree), using incremental
//		distance updates, as described by Arya and Mount in ``Algorithms
//		for fast vector quantization,'' Proc.  of DCC '93: Data Compression
//		Conference, eds. J. A. Storer and M. Cohn, IEEE Press, 1993,
//		381-390.
//
//		The main entry points is annkSearch() which sets things up and
//		then call the recursive routine ann_search().  This is a recursive
//		routine which performs the processing for one node in the kd-tree.
//		There are two versions of this virtual procedure, one for splitting
//		nodes and one for leaves.  When a splitting node is visited, we
//		determine which child to visit first (the closer one), and visit
//		the other child on return.  When a leaf is visited, we compute
//		the distances to the points in the buckets, and update information
//		on the closest points.
//
//		Some trickery is used to incrementally update the distance from
//		a kd-tree rectangle to the query point.  This comes about from
//		the fact that which each successive split, only one component
//		(along the dimension that is split) of the squared distance to
//		the child rectangle is different from the squared distance to
//		the parent rectangle.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//		To keep argument lists short, a number of global variables
//		are maintained which are common to all the recursive calls.
//		These are given below.
//----------------------------------------------------------------------

int				ANNkdDim;				// dimension of space
ANNpoint		ANNkdQ;					// query point
double			ANNkdMaxErr;			// max tolerable squared error
ANNpointArray	ANNkdPts;				// the points
ANNmin_k		*ANNkdPointMK;			// set of k closest points

//----------------------------------------------------------------------
//	annkSearch - search for the k nearest neighbors
//----------------------------------------------------------------------

void ANNkd_tree::annkSearch(
	ANNpoint			q,				// the query point
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps)			// the error bound
{

	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;
	ANNptsVisited = 0;					// initialize count of points visited

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

		ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	// search starting at the root
	root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));

	for (int i = 0; i < k; i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	delete ANNkdPointMK;				// deallocate closest point set
}


//----------------------------------------------------------------------
//	kd_split::ann_search - search a splitting node
//----------------------------------------------------------------------

void ANNkd_split::ann_search(ANNdist box_dist)
{
	// check dist calc term condition
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

	if (!active)
		return;

	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		child[ANN_LO]->ann_search(box_dist);// visit closer child first

		ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
		// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
			ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

		// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_HI]->ann_search(box_dist);

	}
	else {								// right of cutting plane
		child[ANN_HI]->ann_search(box_dist);// visit closer child first

		ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
		// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
			ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

		// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_LO]->ann_search(box_dist);

	}
	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited
}

//----------------------------------------------------------------------
//	kd_leaf::ann_search - search points in a leaf node
//		Note: The unreadability of this code is the result of
//		some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_search(ANNdist box_dist)
{
	if (!active)
		return;

	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		pp = ANNkdPts[bkt[i]];			// first coord of next data point
		qq = ANNkdQ;					// first coord of query point
		dist = 0;

		for(d = 0; d < ANNkdDim; d++) {
			ANN_COORD(1)				// one more coordinate hit
				ANN_FLOP(4)					// increment floating ops

				t = *(qq++) - *(pp++);		// compute length and adv coordinate
			// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
				break;
			}
		}

		if (d >= ANNkdDim &&					// among the k best?
			(ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
				// add it to the list
				ANNkdPointMK->insert(dist, bkt[i]);
				min_dist = ANNkdPointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited
}


//----------------------------------------------------------------------
//	annkSearchContinuous - search for the k nearest neighbors 
//	                       with minimum distance from the query point
//----------------------------------------------------------------------

int ANNkd_tree::annkSearchContinuous(
	ANNpoint			q,				// the query point
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps,			// the error bound
	double				minDist)		// minimum distance of a neighbor	
{

	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

		ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited

	// search starting at the root
	//root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));
	root->ann_searchContinuous(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), minDist);


	for (int i = 0; i < min(k,ANNkdPointMK->n); i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	int found = min(k,ANNkdPointMK->n);
	delete ANNkdPointMK;				// deallocate closest point set
	return found;

}


//----------------------------------------------------------------------
//	kd_split::ann_searchContinuous - search a splitting node
//----------------------------------------------------------------------

void ANNkd_split::ann_searchContinuous(ANNdist box_dist, double minDist)
{
	// check dist calc term condition
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		child[ANN_LO]->ann_searchContinuous(box_dist, minDist);// visit closer child first

		ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
		// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
			ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

		// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_HI]->ann_searchContinuous(box_dist, minDist);

	}
	else {								// right of cutting plane
		child[ANN_HI]->ann_searchContinuous(box_dist, minDist);// visit closer child first

		ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
		// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
			ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

		// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_LO]->ann_searchContinuous(box_dist, minDist);

	}
	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited
}

//----------------------------------------------------------------------
//	kd_leaf::ann_searchContinuous - search points in a leaf node
//		Note: The unreadability of this code is the result of
//		some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_searchContinuous(ANNdist box_dist, double minDist)
{
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		pp = ANNkdPts[bkt[i]];			// first coord of next data point
		qq = ANNkdQ;					// first coord of query point
		dist = 0;

		for(d = 0; d < ANNkdDim; d++) {
			ANN_COORD(1)				// one more coordinate hit
				ANN_FLOP(4)					// increment floating ops

				t = *(qq++) - *(pp++);		// compute length and adv coordinate
			// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
				break;
			}
		}

		if (d >= ANNkdDim &&					 // among the k best?
			(ANN_ALLOW_SELF_MATCH || dist!=0) &&  // and no self-match problem
			dist>ANN_POW(minDist)) {							 // add it to the list
				ANNkdPointMK->insert(dist, bkt[i]);
				min_dist = ANNkdPointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited
}


//----------------------------------------------------------------------
//	kd_leaf::ann_searchContinuous - search points in a leaf node
//		Note: The unreadability of this code is the result of
//		some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_searchContinuous(ANNdist box_dist, checkActiveClass1 *checkActive)
{
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		if (!checkActive->isActiveFast(bkt[i]))
			continue;
		pp = ANNkdPts[bkt[i]];			// first coord of next data point
		qq = ANNkdQ;					// first coord of query point
		dist = 0;

		for(d = 0; d < ANNkdDim; d++) {
			ANN_COORD(1)				// one more coordinate hit
				ANN_FLOP(4)					// increment floating ops

				t = *(qq++) - *(pp++);		// compute length and adv coordinate
			// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
				break;
			}
		}

		if (d >= ANNkdDim &&					 // among the k best?
			(ANN_ALLOW_SELF_MATCH || dist!=0)   // and no self-match problem
			) {							 // add it to the list
				ANNkdPointMK->insert(dist, bkt[i]);
				min_dist = ANNkdPointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited
}

void ANNkd_leaf::ann_searchContinuous(ANNdist box_dist, checkActiveClass2 *checkActive)
{
	if (!active) {
		//ANNptsVisited += n_pts;				// increment number of points visited
		return;
	}
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		if (bkt[i]<0 || !checkActive->isActiveFast(bkt[i]))
			continue;
		pp = ANNkdPts[bkt[i]];			// first coord of next data point
		qq = ANNkdQ;					// first coord of query point
		dist = 0;

		for(d = 0; d < ANNkdDim; d++) {
			ANN_COORD(1)				// one more coordinate hit
				ANN_FLOP(4)					// increment floating ops

				t = *(qq++) - *(pp++);		// compute length and adv coordinate
			// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
				break;
			}
		}

		if (d >= ANNkdDim &&					 // among the k best?
			//(ANN_ALLOW_SELF_MATCH || dist!=0)   // and no self-match problem
				//&& checkActive->isActiveLazy(bkt[i],dist)
					(checkActive->isActiveLazy(dist) || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						//(dist>=checkActive->lastDist || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						) {							 // add it to the list
							ANNkdPointMK->insert(dist, bkt[i]);
							min_dist = ANNkdPointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
	ANN_PTS(n_pts)						// increment points visited
	ANNptsVisited += n_pts;				// increment number of points visited
}


void ANNkd_leaf::ann_searchContinuous(ANNdist box_dist, checkActiveClass3 *checkActive)
{
#ifdef DEBUG
	std::cout << std::hex << "I'm in the leaf (" << (int) this % 0x10000 << ")\n" << std::dec;
#endif
	if (!active) {
		//ANNptsVisited += n_pts;				// increment number of points visited
		return;
	}
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		if (bkt[i]<0 || !checkActive->isActiveFast(bkt[i]))
			continue;
		//if (ANNkdPointMK->n>0 && bkt[i]<checkActive->curNode)
		//	continue;
		pp = ANNkdPts[bkt[i]];			// first coord of next data point
		qq = ANNkdQ;					// first coord of query point
		dist = 0;

		for(d = 0; d < ANNkdDim; d++) {
			ANN_COORD(1)				// one more coordinate hit
				ANN_FLOP(4)					// increment floating ops

				t = *(qq++) - *(pp++);		// compute length and adv coordinate
			// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
				break;
			}
		}

		if (d >= ANNkdDim &&					 // among the k best?
			//(ANN_ALLOW_SELF_MATCH || dist!=0)   // and no self-match problem
				//&& checkActive->isActiveLazy(bkt[i],dist)
					(checkActive->isActiveLazy(dist) || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						//(dist>=checkActive->lastDist || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						) {							 // add it to the list
							ANNkdPointMK->insert(dist, bkt[i]);
							min_dist = ANNkdPointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited
}

void ANNkd_leaf::ann_searchContinuous(ANNdist box_dist, checkActiveClass3 *checkActive, compressed_matrix<double> *distances)
{
	if (!active) {
		//ANNptsVisited += n_pts;				// increment number of points visited
		return;
	}
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		if (bkt[i]<0 || !checkActive->isActiveFast(bkt[i]))
			continue;

		//if (ANNkdPointMK->n>0 && bkt[i]<checkActive->curNode)
		//	continue;
		if ((*distances)(min(bkt[i],checkActive->curNode),max(bkt[i],checkActive->curNode))==0) {
			pp = ANNkdPts[bkt[i]];			// first coord of next data point
			qq = ANNkdQ;					// first coord of query point
			dist = 0;

			for(d = 0; d < ANNkdDim; d++) {
				ANN_COORD(1)				// one more coordinate hit
					ANN_FLOP(4)					// increment floating ops

					t = *(qq++) - *(pp++);		// compute length and adv coordinate
				// exceeds dist to k-th smallest?
				if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
					break;
				}
			}
			if (d>= ANNkdDim)
				(*distances)(min(bkt[i],checkActive->curNode),max(bkt[i],checkActive->curNode)) = dist;
		} else {
			dist = (*distances)(min(bkt[i],checkActive->curNode),max(bkt[i],checkActive->curNode));
			d = ANNkdDim;
		}

		if (d >= ANNkdDim &&					 // among the k best?
			//(ANN_ALLOW_SELF_MATCH || dist!=0)   // and no self-match problem
				//&& checkActive->isActiveLazy(bkt[i],dist)
					(checkActive->isActiveLazy(dist) || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						//(dist>=checkActive->lastDist || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						) {							 // add it to the list
							ANNkdPointMK->insert(dist, POSITIVE(bkt[i]));
							min_dist = ANNkdPointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited
}



void ANNkd_leaf::ann_searchContinuous_BottomUp(ANNdist box_dist, checkActiveClass2 *checkActive, ANNkd_ptr caller, int nodeIdx)//, ANNkd_ptr *searchStart)
{
	if (!active) {
		//ANNptsVisited += n_pts;				// increment number of points visited
		return;
	}
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";

	if (n_pts==1) {
		/*if (checkActive->redirectIdx[bkt[0]]!=checkActive->curNode)			// can be removed
		cout << "Error, you have started from incorrect leaf node!\a";
		*/
		ANNptsVisited++;				// increment number of points visited
		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_searchContinuous_BottomUp(0, checkActive, this, nodeIdx);//, searchStart);
		return;
	} else {
		for (int i = 0; i < n_pts; i++) {	// check points in bucket

			if (bkt[i]<0 || !checkActive->isActiveFast(bkt[i]))
				continue;
			pp = ANNkdPts[bkt[i]];			// first coord of next data point
			qq = ANNkdQ;					// first coord of query point
			dist = 0;

			for(d = 0; d < ANNkdDim; d++) {
				ANN_COORD(1)				// one more coordinate hit
					ANN_FLOP(4)					// increment floating ops

					t = *(qq++) - *(pp++);		// compute length and adv coordinate
				// exceeds dist to k-th smallest?
				if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
					break;
				}
			}

			if (d >= ANNkdDim &&					 // among the k best?
				//(ANN_ALLOW_SELF_MATCH || dist!=0) &&  // and no self-match problem
					//&& checkActive->isActiveLazy(bkt[i],dist)
						(checkActive->isActiveLazy(dist) || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
							//(dist>=checkActive->lastDist || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
							) {							 // add it to the list
								ANNkdPointMK->insert(dist, POSITIVE(bkt[i]));
								min_dist = ANNkdPointMK->max_key();
			}
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited

	if (parent!=NULL && parent!=KD_TRIVIAL)
		parent->ann_searchContinuous_BottomUp(0, checkActive,this, nodeIdx);//, searchStart);
}

void ANNkd_leaf::ann_searchContinuous_BottomUp(ANNdist box_dist, checkActiveClass3 *checkActive, ANNkd_ptr caller, int nodeIdx)//, ANNkd_ptr *searchStart)
{
#ifdef DEBUG
	std::cout << std::hex << "I'm in the leaf (" << (int)this % 0x10000 << ")\n" << std::dec;
#endif
	if (!active) {
		//ANNptsVisited += n_pts;				// increment number of points visited
		return;
	}
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";

	if (n_pts==1) {
		/*if (checkActive->redirectIdx[bkt[0]]!=checkActive->curNode)			// can be removed
		cout << "Error, you have started from incorrect leaf node!\a";
		*/
		ANNptsVisited++;				// increment number of points visited
		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_searchContinuous_BottomUp(0, checkActive, this, nodeIdx);//, searchStart);
		return;
	} else {
		for (int i = 0; i < n_pts; i++) {	// check points in bucket
			if (bkt[i]<0 || !checkActive->isActiveFast(bkt[i]))
				continue;
			pp = ANNkdPts[bkt[i]];			// first coord of next data point
			qq = ANNkdQ;					// first coord of query point

			dist = 0;
			for(d = 0; d < ANNkdDim; d++) {
				ANN_COORD(1)				// one more coordinate hit
					ANN_FLOP(4)					// increment floating ops

					t = *(qq++) - *(pp++);		// compute length and adv coordinate
				// exceeds dist to k-th smallest?
				if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
					break;
				}
			}

			if (d >= ANNkdDim &&					 // among the k best?
				//(ANN_ALLOW_SELF_MATCH || dist!=0) &&  // and no self-match problem
					//checkActive->isActiveLazy(bkt[i],dist)
						(checkActive->isActiveLazy(dist) || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						) {							 // add it to the list
							ANNkdPointMK->insert(dist, POSITIVE(bkt[i]));
							min_dist = ANNkdPointMK->max_key();
			}
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited

	if (parent!=NULL && parent!=KD_TRIVIAL)
		parent->ann_searchContinuous_BottomUp(0, checkActive,this, nodeIdx);//, searchStart);
}

void ANNkd_leaf::ann_search_BottomUp(ANNdist box_dist, ANNkd_ptr caller, int nodeIdx)//, ANNkd_ptr *searchStart)
{
#ifdef DEBUG
	std::cout << std::hex << "I'm in the leaf (" << (int)this % 0x10000 << ")\n" << std::dec;
#endif
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";

	if (n_pts==1) {
		/*if (checkActive->redirectIdx[bkt[0]]!=checkActive->curNode)			// can be removed
		cout << "Error, you have started from incorrect leaf node!\a";
		*/
		ANNptsVisited++;				// increment number of points visited
		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_search_BottomUp(0, this, nodeIdx);//, searchStart);
		return;
	} else {
		for (int i = 0; i < n_pts; i++) {	// check points in bucket
			pp = ANNkdPts[bkt[i]];			// first coord of next data point
			qq = ANNkdQ;					// first coord of query point

			dist = 0;
			for(d = 0; d < ANNkdDim; d++) {
				ANN_COORD(1)				// one more coordinate hit
					ANN_FLOP(4)					// increment floating ops

					t = *(qq++) - *(pp++);		// compute length and adv coordinate
				// exceeds dist to k-th smallest?
				if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
					break;
				}
			}

			if (d >= ANNkdDim					 // among the k best?
				//(ANN_ALLOW_SELF_MATCH || dist!=0) &&  // and no self-match problem
					//checkActive->isActiveLazy(bkt[i],dist)
						// && (checkActive->isActiveLazy(dist) || ANNkdMaxErr>1)   //*__  && ANNkdPointMK->n)<=0) 
						) {							 // add it to the list
							ANNkdPointMK->insert(dist, POSITIVE(bkt[i]));
							min_dist = ANNkdPointMK->max_key();
			}
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
		ANN_PTS(n_pts)						// increment points visited
		ANNptsVisited += n_pts;				// increment number of points visited

	if (parent!=NULL && parent!=KD_TRIVIAL)
		parent->ann_search_BottomUp(0,this, nodeIdx);//, searchStart);
}



bool ANNkd_leaf::annDeleteNode(int nodeIdx, ANNdist box_dist)
{
	/*
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;
	*/
	//min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	//cout << "searching a leaf node with " << n_pts << " nodes.\n";
	if (!active)
		return false;

	//cout << " ==>> " << this << "(" << bkt[0] << ")\n ";

	if (n_pts==1) {
		if (bkt[0]==nodeIdx) {
			bkt[0]=~bkt[0];
			active = false;
			return true;
		}
		return false;
	}
	else {
		bool flag = false;
		int activeNodesCount = 0;
		for (int i = 0; i < n_pts; i++) {	// check points in bucket
			if (bkt[i]==nodeIdx) {
				bkt[i]=~bkt[i];			//  negative value denotes inactive node
				flag = true;
			} else if (bkt[i]>0) 
				activeNodesCount++;

			/*
			pp = ANNkdPts[abs(bkt[i])];			// first coord of next data point
			qq = ANNkdQ;					// first coord of query point
			dist = 0;

			for(d = 0; d < ANNkdDim; d++) {

			t = *(qq++) - *(pp++);		// compute length and adv coordinate
			// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
			break;
			}
			}

			if (d >= ANNkdDim)					 // among the k best?
			{
			ANNkdPointMK->insert(dist, abs(bkt[i]));
			min_dist = ANNkdPointMK->max_key();
			}
			*/
		}
		if (activeNodesCount==0) {
			active = false;
		}
		return flag;
	}
}



bool ANNkd_leaf::annDeleteNode_BottomUp(int nodeIdx)
{
	if (!active)
		return false;					// It is not required to make me inactive

	//cout << "I am deleting " << std::hex << (int)this%0x10000 << std::dec << endl;

	bool deleteResult = false;
	
	if (n_pts==1) {
		if (bkt[0]==nodeIdx) {			// can be removed, it must be true always
			bkt[0]=~bkt[0];
			active = false;
			if (parent!=NULL && parent!=KD_TRIVIAL && parent->parent!=NULL && parent->parent!=KD_TRIVIAL) {
				deleteResult = parent->annDeleteNode_BottomUp(nodeIdx);
				if (deleteResult) {
					delete parent;
					parent = 0;
				}
			}
			return deleteResult;		// if I cannot delete the parent, I cannot delete the child
		}
		return false;					// This indicates an error!
	}
	else {
		//bool flag = false;
		int activeNodesCount = 0;
		for (int i = 0; i < n_pts; i++) {	// check points in bucket
			if (bkt[i]==nodeIdx) {
				bkt[i]=~bkt[i];				//  negative value denotes inactive node
				//flag = true;
			} else if (bkt[i]>0) 
				activeNodesCount++;
		}
		if (!activeNodesCount) {
			active = false;
			if (parent!=NULL && parent!=KD_TRIVIAL && parent->parent!=NULL && parent->parent!=KD_TRIVIAL)
				deleteResult = parent->annDeleteNode_BottomUp(nodeIdx);
			if (deleteResult) {
				delete parent;
				parent = 0;
			}
			return deleteResult;
		}
		return false;
	}
	return false;					// this does not occur and is added only for completeness
}




//----------------------------------------------------------------------
//	kd_split::ann_searchContinuous - search a splitting node
//----------------------------------------------------------------------

void ANNkd_split::ann_searchContinuous(ANNdist box_dist, checkActiveClass1 *checkActive)
{
	// check dist calc term condition
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
		// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
			ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

		// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);

	}
	else {								// right of cutting plane
		child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
		// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
			ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

		// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);

	}
	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited
}

void ANNkd_split::ann_searchContinuous(ANNdist box_dist, checkActiveClass2 *checkActive)
{
	//if (!active) {							// if current node is active
	//	cout << "This is inactive.";
	//	return;
	//}

	// check dist calc term condition
	//if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

	/*if (!(child[ANN_LO]->active || child[ANN_HI]->active))
	return;*/
	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		if (child[ANN_LO]->active)
			child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_HI]->active) {
			ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);
		}
	}
	else {								// right of cutting plane
		if (child[ANN_HI]->active)
			child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_LO]->active) {
			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);
		}
	}

	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited
}

void ANNkd_split::ann_searchContinuous(ANNdist box_dist, checkActiveClass3 *checkActive)
{
#ifdef DEBUG
	std::cout << std::hex << "I'm in the split (" << (int) this % 0x10000 << ")\n" << std::dec;
#endif


	//if (!active) {							// if current node is active
	//	cout << "This is inactive.";
	//	return;
	//}

	// check dist calc term condition
	//if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

	/*if (!(child[ANN_LO]->active || child[ANN_HI]->active))
	return;*/
	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		if (child[ANN_LO]->active)
			child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_HI]->active) {
			ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * (((ANNkdPointMK->n)>0)?ANNkdMaxErr:1) < ANNkdPointMK->max_key())
				child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);
		}
	}
	else {								// right of cutting plane
		if (child[ANN_HI]->active)
			child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_LO]->active) {
			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * (((ANNkdPointMK->n)>0)?ANNkdMaxErr:1) < ANNkdPointMK->max_key())
				child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);
		}
	}

	ANN_FLOP(10)						// increment floating ops
	ANN_SPL(1)							// one more splitting node visited
}

void ANNkd_split::ann_searchContinuous(ANNdist box_dist, checkActiveClass3 *checkActive, compressed_matrix<double>*distances)
{
	//if (!active) {							// if current node is active
	//	cout << "This is inactive.";
	//	return;
	//}

	// check dist calc term condition
	//if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

	/*if (!(child[ANN_LO]->active || child[ANN_HI]->active))
	return;*/
	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		if (child[ANN_LO]->active)
			child[ANN_LO]->ann_searchContinuous(box_dist, checkActive,distances);// visit closer child first

		if (child[ANN_HI]->active) {
			ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[ANN_HI]->ann_searchContinuous(box_dist, checkActive,distances);
		}
	}
	else {								// right of cutting plane
		if (child[ANN_HI]->active)
			child[ANN_HI]->ann_searchContinuous(box_dist, checkActive,distances);// visit closer child first

		if (child[ANN_LO]->active) {
			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[ANN_LO]->ann_searchContinuous(box_dist, checkActive,distances);
		}
	}

	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited
}

void ANNkd_split::ann_searchContinuous_BottomUp(ANNdist box_dist, checkActiveClass2 *checkActive, ANNkd_ptr caller, int nodeIdx)//, ANNkd_ptr *searchStart)
{
	//if (!active) {							// if current node is active
	//	cout << "This is inactive.";
	//	return;
	//}

	// check dist calc term condition
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited && ANNkdPointMK->n > 0) return;

	if (caller==NULL) {
		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_searchContinuous_BottomUp(0, checkActive,this, nodeIdx);//, searchStart);
		return;
	}


	/*if (!(child[ANN_LO]->active || child[ANN_HI]->active))
	return;*/
	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (caller==child[ANN_LO]) {					// left of cutting plane is searched before
		//if (child[ANN_LO]->active)	// searched before
		//	child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_HI]->active) {
			ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);
		}
	}
	else {								// right of cutting plane
		//if (child[ANN_HI]->active)
		//	child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_LO]->active) {
			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
				child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);
		}
	}

	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited

		if ((ANNkdPointMK->n) == 0) {
#ifdef DEBUG
			//cout << "search start changed for " << std::dec << nodeIdx << " from " << std::hex << (int)searchStart[nodeIdx]%0x10000 << std::dec << " to " << std::hex << (int)this % 0x10000 << endl;
			//cout << "search start changed for " << std::dec << nodeIdx << " from " << std::hex << (int)searchStart[nodeIdx]%0x10000 << std::dec << " to " << std::hex << (int)this % 0x10000 << endl;
#endif
			//searchStart[nodeIdx] = this;
		}

		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_searchContinuous_BottomUp(0, checkActive,this, nodeIdx);//, searchStart);
}

void ANNkd_split::ann_searchContinuous_BottomUp(ANNdist box_dist, checkActiveClass3 *checkActive, ANNkd_ptr caller, int nodeIdx)//, ANNkd_ptr *searchStart)
{
	//if (!active) {							// if current node is active
	//	cout << "This is inactive.";
	//	return;
	//}

	// check dist calc term condition
#ifdef DEBUG
	std::cout << std::hex << "I'm in the split (" << (int)this % 0x10000 << ")\n" << std::dec;
#endif
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited && ANNkdPointMK->n > 0) return;

	if (caller==NULL) {
		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_searchContinuous_BottomUp(0, checkActive,this, nodeIdx);//, searchStart);
		return;
	}


	/*if (!(child[ANN_LO]->active || child[ANN_HI]->active))
	return;*/
	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (caller==child[ANN_LO]) {					// left of cutting plane is searched before
		//if (child[ANN_LO]->active)	// searched before
		//	child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_HI]->active) {
			ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * (((ANNkdPointMK->n)>0)?ANNkdMaxErr:1) < ANNkdPointMK->max_key())
				child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);
		}
	}
	else {								// right of cutting plane
		//if (child[ANN_HI]->active)
		//	child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_LO]->active) {
			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * (((ANNkdPointMK->n)>0)?ANNkdMaxErr:1) < ANNkdPointMK->max_key())
				child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);
		}
	}

	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited

		if ((ANNkdPointMK->n) == 0) {
#ifdef DEBUG
			//cout << "search start changed for " << std::dec << nodeIdx << " from " << std::hex << (int)searchStart[nodeIdx]%0x10000 << " to " << std::hex << (int)this % 0x10000 << std::dec << endl;
#endif
			//searchStart[nodeIdx] = this;
		}

		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_searchContinuous_BottomUp(0, checkActive,this, nodeIdx);//, searchStart);
}

void ANNkd_split::ann_search_BottomUp(ANNdist box_dist, ANNkd_ptr caller, int nodeIdx)//, ANNkd_ptr *searchStart)
{
	//if (!active) {							// if current node is active
	//	cout << "This is inactive.";
	//	return;
	//}

	// check dist calc term condition
#ifdef DEBUG
	std::cout << std::hex << "I'm in the split (" << (int)this % 0x10000 << ")\n" << std::dec;
#endif
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited && ANNkdPointMK->n > 0) return;

	if (caller==NULL) {
		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_search_BottomUp(0,this, nodeIdx);//, searchStart);
		return;
	}


	/*if (!(child[ANN_LO]->active || child[ANN_HI]->active))
	return;*/
	// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (caller==child[ANN_LO]) {					// left of cutting plane is searched before
		//if (child[ANN_LO]->active)	// searched before
		//	child[ANN_LO]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_HI]->active) {
			ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * (((ANNkdPointMK->n)>0)?ANNkdMaxErr:1) < ANNkdPointMK->max_key())
				child[ANN_HI]->ann_search(box_dist);
		}
	}
	else {								// right of cutting plane
		//if (child[ANN_HI]->active)
		//	child[ANN_HI]->ann_searchContinuous(box_dist, checkActive);// visit closer child first

		if (child[ANN_LO]->active) {
			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
			if (box_diff < 0)				// within bounds - ignore
				box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			if (box_dist * (((ANNkdPointMK->n)>0)?ANNkdMaxErr:1) < ANNkdPointMK->max_key())
				child[ANN_LO]->ann_search(box_dist);
		}
	}

	ANN_FLOP(10)						// increment floating ops
		ANN_SPL(1)							// one more splitting node visited

		if ((ANNkdPointMK->n) == 0) {
#ifdef DEBUG
			//cout << "search start changed for " << std::dec << nodeIdx << " from " << std::hex << (int)searchStart[nodeIdx]%0x10000 << " to " << std::hex << (int)this % 0x10000 << std::dec << endl;
#endif
			//searchStart[nodeIdx] = this;
		}

		if (parent!=NULL && parent!=KD_TRIVIAL)
			parent->ann_search_BottomUp(0, this, nodeIdx);//, searchStart);
}

bool ANNkd_split::annDeleteNode(int nodeIdx, ANNdist box_dist)
{
	/*
	if (!active) {							// if current node is active
	return false;
	}

	if (!(child[ANN_LO]->active || child[ANN_HI]->active)) {
	active = false;
	return false;
	}
	*/									// distance to cutting plane


	if (!active) {
		return false;
		//cout << "This is not active. dont call it." << nodeIdx << " ";
	}

	//cout << " -> " << this;

	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	bool flag = false;
	if (cut_diff < 0) {					// left of cutting plane
		flag = child[ANN_LO]->active && child[ANN_LO]->annDeleteNode(nodeIdx, box_dist); // visit closer child first

		/*if (!flag && child[ANN_HI]->active)
		{
		ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];

		if (box_diff < 0)				// within bounds - ignore
		box_diff = 0;
		// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
		ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

		// visit further child if close enough
		//if (box_dist == 0) // ANNkdPointMK->max_key())
		flag = child[ANN_HI]->annDeleteNode(nodeIdx, box_dist);

		}*/
	}
	else {								// right of cutting plane
		flag = child[ANN_HI]->active && child[ANN_HI]->annDeleteNode(nodeIdx, box_dist); // visit closer child first

		if (!flag && cut_diff==0 && ANNkdQ[cut_dim] <= cd_bnds[ANN_HI] && child[ANN_LO]->active) 
		{
			/*
			ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
			if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
			// distance to further box
			box_dist = (ANNdist) ANN_SUM(box_dist,
			ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

			// visit further child if close enough
			//if (box_dist == 0) // * ANNkdMaxErr < ANNkdPointMK->max_key())
			*/
			flag = child[ANN_LO]->annDeleteNode(nodeIdx, box_dist);
		}
	}
	if (flag)
		active = child[ANN_LO]->active || child[ANN_HI]->active;
	return flag;
}


bool ANNkd_split::annDeleteNode_BottomUp(int nodeIdx)
{
	
	//if (!active) {							// if current node is active
	//return false;
	//}

	//if (!(child[ANN_LO]->active || child[ANN_HI]->active)) {
	//active = false;
	//return false;
	//}
	//									// distance to cutting plane


	if (!active) {
		//_/return false;
		return true;
		//cout << "This is not active. dont call it." << nodeIdx << " ";
	}

	//cout << "I am deleting " << std::hex << (int)this%0x10000 << std::dec << endl;
	
	bool deleteResult = false;

	if (!(child[ANN_LO]->active || child[ANN_HI]->active)) {
		active = false;
		if (parent!=NULL && parent!=KD_TRIVIAL) 
			deleteResult = parent->annDeleteNode_BottomUp(nodeIdx);
			//_/return (parent->annDeleteNode_BottomUp(nodeIdx));
		//_/return false;
		if (deleteResult) {
			delete parent;
			parent = 0;
		}
		return true;
	} 
	
	

	ANNkd_split *curNode = this;
	ANNkd_split *splitParent = (ANNkd_split *) parent;

	if (parent==NULL || parent==KD_TRIVIAL)// || splitParent->parent==NULL || splitParent->parent==KD_TRIVIAL)
		return false;

	if (child[ANN_LO]->active) {
		if (curNode==splitParent->child[ANN_LO]) {
			splitParent->child[ANN_LO] = child[ANN_LO];
			//child[ANN_LO]->parent = splitParent;
		} else {
			splitParent->child[ANN_HI] = child[ANN_LO];
			//child[ANN_LO]->parent = splitParent;
		}
		
		//delete this;
		child[ANN_LO]->parent = parent;
		child[ANN_LO] = child[ANN_HI] = 0;
		//_/return false;
		return true;
	} else if (child[ANN_HI]->active) {
		if (curNode==splitParent->child[ANN_LO]) {
			splitParent->child[ANN_LO] = child[ANN_HI];
			//child[ANN_HI]->parent = splitParent;
		} else {
			splitParent->child[ANN_HI] = child[ANN_HI];
			//child[ANN_HI]->parent = splitParent;
		}
		//delete this;
		child[ANN_HI]->parent = parent;
		child[ANN_LO] = child[ANN_HI] = 0;
		//_/return false;
		return true;
	}

	return false; // only for completeness, does not occur
}

/*
bool ANNkd_split::annDeleteNode_BottomUp(int nodeIdx)
{
	/*
	if (!active) {							// if current node is active
	return false;
	}

	if (!(child[ANN_LO]->active || child[ANN_HI]->active)) {
	active = false;
	return false;
	}
	*/									// distance to cutting plane


	//if (!active) {
	//	return false;
	//	//cout << "This is not active. dont call it." << nodeIdx << " ";
	//}

	//if (!(child[ANN_LO]->active || child[ANN_HI]->active)) {
	//	active = false;
	//	if (parent!=NULL && parent!=KD_TRIVIAL)
	//		return (parent->annDeleteNode_BottomUp(nodeIdx));
	//	return false;
	//} 
	//
	//if (parent==NULL || parent==KD_TRIVIAL)
	//	return false;

	//ANNkd_split *curNode = this;
	//ANNkd_split *splitParent = (ANNkd_split *) parent;

	//bool childDirection;

	//if (curNode==splitParent->child[ANN_LO])
	//	childDirection = ANN_LO;
	//else
	//	childDirection = ANN_HI;

	//do {
	//	if (child[ANN_LO]->active) {
	//		if (childDirection==ANN_LO) {
	//			splitParent->child[ANN_LO] = child[ANN_LO];
	//		} else {
	//			splitParent->child[ANN_HI] = child[ANN_LO];
	//		}
	//		//child[ANN_LO]->parent = splitParent;
	//	} else if (child[ANN_HI]->active) {
	//		if (childDirection==ANN_LO) {
	//			splitParent->child[ANN_LO] = child[ANN_HI];
	//		} else {
	//			splitParent->child[ANN_HI] = child[ANN_HI];
	//		}
	//		//child[ANN_HI]->parent = splitParent;
	//	}

	//	if ((splitParent->child[ANN_LO]->active==true && splitParent->child[ANN_HI]->active==true) || splitParent->parent==NULL)
	//		break;
	//	cout << "+";
	//	curNode = splitParent;
	//	splitParent = (ANNkd_split *)splitParent->parent;
	//	if (curNode==splitParent->child[ANN_LO])
	//		childDirection = ANN_LO;
	//	else
	//		childDirection = ANN_HI;
	//} while (true);

	//if (child[ANN_LO]->active) 
	//	child[ANN_LO]->parent = splitParent;
	//else
	//	child[ANN_HI]->parent = splitParent;


	
	/*if (child[ANN_LO]->active) {
		if (childDirection==ANN_LO) {
			splitParent->child[ANN_LO] = child[ANN_LO];
		} else {
			splitParent->child[ANN_HI] = child[ANN_LO];
		}
		child[ANN_LO]->parent = splitParent;
		return false;
	} else if (child[ANN_HI]->active) {
		if (childDirection==ANN_LO) {
			splitParent->child[ANN_LO] = child[ANN_HI];
		} else {
			splitParent->child[ANN_HI] = child[ANN_HI];
		}
		child[ANN_HI]->parent = splitParent;
		return false;
	}
	
	return false;
}
*/

int ANNkd_tree::annkSearchContinuous(
	ANNpoint			q,				// the query point
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps,			// the error bound
	checkActiveClass1	*checkActive)	// check for active (acceptable) node
{

	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

		ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited

	// search starting at the root
	//root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));
	root->ann_searchContinuous(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), checkActive);


	for (int i = 0; i < min(k,ANNkdPointMK->n); i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	int found = min(k,ANNkdPointMK->n);
	delete ANNkdPointMK;				// deallocate closest point set
	return found;

}

int ANNkd_tree::annkSearchContinuous(
	ANNpoint			q,				// the query point
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps,			// the error bound
	checkActiveClass2	*checkActive)	// check for active (acceptable) node
{
	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

		ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited

	// search starting at the root
	//root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));
	root->ann_searchContinuous(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), checkActive);


	for (int i = 0; i < min(k,ANNkdPointMK->n); i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	int found = min(k,ANNkdPointMK->n);
	delete ANNkdPointMK;				// deallocate closest point set
	return found;

}


int ANNkd_tree::annkSearchContinuous(
	ANNpoint			q,				// the query point
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps,			// the error bound
	checkActiveClass3	*checkActive)	// check for active (acceptable) node
{
	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

		ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited

	// search starting at the root
	//root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));
	root->ann_searchContinuous(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), checkActive);


	for (int i = 0; i < min(k,ANNkdPointMK->n); i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	int found = min(k,ANNkdPointMK->n);
	delete ANNkdPointMK;				// deallocate closest point set
	return found;

}


int ANNkd_tree::annkSearchContinuous_BottomUp(
	int					nodeIdx,		// the index of the query point (between all nodes in the kdtree)
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps,			// the error bound
	checkActiveClass2	*checkActive)	// check for active (acceptable) node
{
	ANNpoint q = pts[nodeIdx];
	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

		ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited

	// search starting at the root
	//root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));
	//root->ann_searchContinuous(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), checkActive);

#ifdef DEBUG
	//std::cout << "** Search for " << std::dec << nodeIdx << " is started from " << std::hex << (int)searchStart[nodeIdx] % 0x10000 << endl;
	//Print(ANNtrue, cout);
#endif

	//searchStart[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx, searchStart);
	//searchStart[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx, searchStart);
	
	//leafAddr[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx, searchStart);
	leafAddr[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx);


	for (int i = 0; i < min(k,ANNkdPointMK->n); i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	int found = min(k,ANNkdPointMK->n);
	delete ANNkdPointMK;				// deallocate closest point set
	return found;

}



int ANNkd_tree::annkSearchContinuous_BottomUp(
	int					nodeIdx,		// the index of the query point (between all nodes in the kdtree)
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps,			// the error bound
	checkActiveClass3	*checkActive)	// check for active (acceptable) node
{
	ANNpoint q = pts[nodeIdx];
	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

	ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited

	// search starting at the root
	//root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));
	//root->ann_searchContinuous(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), checkActive);

#ifdef DEBUG
	//std::cout << "** Search for " << std::dec << nodeIdx << " is started from " << std::hex << (int)searchStart[nodeIdx] % 0x10000 << endl;
	//Print(ANNtrue, cout);
#endif

	//searchStart[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx, searchStart);
	
	//leafAddr[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx, searchStart);
	leafAddr[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx);


	for (int i = 0; i < min(k,ANNkdPointMK->n); i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	int found = min(k,ANNkdPointMK->n);
	delete ANNkdPointMK;				// deallocate closest point set
	return found;

}

int ANNkd_tree::annkSearch_BottomUp(
	int					nodeIdx,		// the index of the query point (between all nodes in the kdtree)
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	double				eps)			// the error bound
{
	ANNpoint q = pts[nodeIdx];
	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)							// increment floating op count

	ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited

	// search starting at the root
	//root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));
	//root->ann_searchContinuous(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), checkActive);

#ifdef DEBUG
	//std::cout << "** Search for " << std::dec << nodeIdx << " is started from " << std::hex << (int)searchStart[nodeIdx] % 0x10000 << endl;
	//Print(ANNtrue, cout);
#endif

	//searchStart[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx, searchStart);
	
	//leafAddr[nodeIdx]->ann_searchContinuous_BottomUp(0,checkActive,0, nodeIdx, searchStart);
	leafAddr[nodeIdx]->ann_search_BottomUp(0,0, nodeIdx);


	for (int i = 0; i < min(k,ANNkdPointMK->n); i++) {		// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
	int found = min(k,ANNkdPointMK->n);
	delete ANNkdPointMK;				// deallocate closest point set
	return found;

}
 

bool ANNkd_tree::annDeleteNode(
	int					nodeIdx)				// the index of the node to be deleted
{


	ANNkdDim = dim;						// copy arguments to static equivs
	ANNkdPts = pts;
	ANNpoint q = ANNkdPts[nodeIdx];
	ANNkdQ = q;
	//int k = 1;

	/*
	if (nodeIdx > n_pts) {					// index is incorrect
	annError("The index of the node to be deleted is more than data points", ANNabort);
	}


	ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited
	*/

	//return root->annDeleteNode(nodeIdx,annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));

	/*std::cout << " ** " << nodeIdx << " ** ";
	if (nodeIdx==5)
	Print(ANNtrue,cout);*/

	return root->annDeleteNode(nodeIdx,0);

}


bool ANNkd_tree::annDeleteNode_BottomUp(
	int					nodeIdx)				// the index of the node to be deleted starting from the leaf node
{
	//ANNkdDim = dim;						// copy arguments to static equivs
	//ANNkdPts = pts;
	//ANNpoint q = ANNkdPts[nodeIdx];
	//ANNkdQ = q;

	//int k = 1;


	/*
	if (nodeIdx > n_pts) {					// index is incorrect
	annError("The index of the node to be deleted is more than data points", ANNabort);
	}


	ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
	ANNptsVisited = 0;					// initialize count of points visited
	*/

	//return root->annDeleteNode(nodeIdx,annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));

#ifdef DEBUG
	std::cout << " ** Deleting " << nodeIdx << " ** ";
#endif
	

	//cout << "I am deleting " << std::hex << (int)this%0x10000 << std::dec << endl;

	//return root->annDeleteNode(nodeIdx);
	bool deleteResult = false;
	
	deleteResult = leafAddr[nodeIdx]->annDeleteNode_BottomUp(nodeIdx);

	if (deleteResult) {
		delete leafAddr[nodeIdx];
		leafAddr[nodeIdx] = 0;
	}


	return deleteResult;

}


