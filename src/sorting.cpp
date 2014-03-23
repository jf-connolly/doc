/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  sort.h - Descending heap sort algorithm (with tie breaker)                           *
 *   																		  *
 *  ----------------------------------------------------------------------    *
 *  This program is free software; you can redistribute it and/or modify it   *
 *  under the terms of the GNU General Public License as published by the     *
 *  Free Software Foundation; either version 2 of the License, or (at your    *
 *  option) any later version.                                   			  *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful, but       *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of                *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         *
 *  General Public License for more details.                                  *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License along   *
 *  with this program; if not, write to the Free Software Foundation, Inc.,   *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.                 *
 ******************************************************************************/

/*------------------------------- Header files -------------------------------*/
//-- Standard
#include <cstdlib>
#include <stdio.h>

//-- Project specific
#include "sorting.h"

/*============================================================================*
 *  Function name			:	sorting::sort
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Default constructor
 *============================================================================*/
sorting::sorting()
{
	size       = 0;
	fitness    = (float*)NULL;
	tieBreaker = (float*)NULL;
	order      = (int*)  NULL;
}
/*============================================================================*
 *  Function name			:	sorting::initialization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Initialization function
 *----------------------------------------------------------------------------*
 *  Argument				:   int sz	- Size of the
 *============================================================================*/
void sorting::initialization(int sz)
{
	size       = sz;
	fitness    = new float[size]; //- 1st objective (most important)
	tieBreaker = new float[size]; //- 2nd objective (tie breaker)
	order      = new int  [size];
}
/*============================================================================*
 *  Function name			:	sorting::~sort
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Destructor
 *============================================================================*/
sorting::~sorting()
{
	size = 0;
	delete [] fitness;
	delete [] tieBreaker;
	delete [] order;
}
/*============================================================================*
 *  Function name			:	sorting::heapSort
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Sort using a heap data structure
 *  							(brought to you by wikipedia!)
 *============================================================================*/
void sorting::heapSort()
{
	int end = size -1;

	//-- Order initialization
	for(int i=0; i<size; i++)   { order[i]=i; }

	//-- Initial heap
	doHeap();

	//-- Sorting the heap
	while (end > 0)
	{	//-- Swap the root (largest value) with the last element of the heap
		swap(&order[0], &order[end]);

		//-- Put the heap back in max-heap order
		end -= 1;
		siftDown(0,end);
	}
}
/*============================================================================*
 *  Function name			:	sorting::doHeap
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Does a heap structure for the heap sort
 *============================================================================*/
void sorting::doHeap()
{
	//-- Index in a of the last parent node
	int start = (size-2) /2;

	//-- Sift down the node at index start to the proper place such that all
	//   nodes below the start index are in heap order
	while ( start >= 0 )
	{	siftDown(start, size-1);
		start -= 1;
	}
}
/*============================================================================*
 *  Function name			:	sorting::siftDown
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Goes into the heap
 *----------------------------------------------------------------------------*
 *  Argument				:   int start - Start of the heap
 * 								int end   - How far down the heap to sift
 *============================================================================*/
void sorting::siftDown(int start, int end)
{
	int root = start;

	//-- While the root has at least one child
	while ( root *2 +1 <= end )
	{
		int child = root *2 +1;
		//-- If the child has a sibling ...
		if(child +1 <= end)
		{	//-- ... and the child's value is less than its sibling's
			//   (with a tie breaker), point to the sibling instead
			if( fitness[order[child]] < fitness[order[child+1]])
			{	child = child +1;   }
			else if( (fitness   [order[child]] == fitness   [order[child+1]]) &&
			         (tieBreaker[order[child]] <  tieBreaker[order[child+1]]) )
			{	child = child +1;   }
		}

		//-- Out of max-heap order (with tie breaker)
		if( fitness[order[root]] < fitness[order[child]]   ||
			( fitness   [order[root]] == fitness   [order[child]] &&
			  tieBreaker[order[root]] <  tieBreaker[order[child]] )  )
		{	swap(&order[root], &order[child]);
			root = child;
		}
		else
		{	return;   }
	}
}
/*============================================================================*
 *  Function name			:	knn::swap
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Swap two values
 *----------------------------------------------------------------------------*
 *  Argument				:   int *a, *b  - Values to swap
 *============================================================================*/
void sorting::swap( int *a, int *b )   { int c = *a;   *a = *b;   *b = c; }
