/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  artmap.h - ARTMAP network (ony fuzzy ARTMAP for now).                      *
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
//-- (none)

//-- Project specific
#include "configuration.h"
#include "dbase.h"
#include "result.h"

#ifndef ARTMAP_H
#define ARTMAP_H

using namespace std;

/*-------------------------------- Constants ---------------------------------*/
//-- General
#define defaultConvTime		10

//-- Default values : Fuzzy ARTMAP
#define defaultAlpha	0.001
#define defaultBeta		1.0
#define defaultEpsilon	0.001
#define defaultRho		0

#define opTrain  1
#define opTest   2

//-- Macro
#define argMax(val_a,val_b, a,b)   ( ((val_a>val_b)||(val_a==val_b))?a:b )
#define argMin(val_a,val_b, a,b)   ( ((val_a<val_b)||(val_a==val_b))?a:b )

class artmap {
public:
	/*------------------------------ Variables -------------------------------*/
	//-- General informations
	int nClasses;			//- Nb. of classes
	int nFeatures;			//- Nb. of features
	int sizeF2max;			//- Max. nb. of allocated neurons in memory
	int nEpochMax;			//- Maximal nb. of training epochs

	int  nEpochs;			//- Nb. of completed training epochs
	int  nPatternsLearned;	//- Nb. of patterns learned
	int  typeOp;			//- Train or test
	int *order;				//- Patterns presentation order
	int  sizeOrder;			//- Size of the presentation order vector
	int	 error;				//- Error code

	//-- Parametres
	float alpha;			//- Choice parameter
	float beta;				//- Learning parameter
	float epsilon;			//- Match tracking parameter
	float rho;				//- Baseline vigilance

	//-- F0 layer
	int	sizeF0;				//- F0 layer size (nb. of features)

	//-- F1 layer
	int	sizeF1;				//- F1 layer size (2 x nb. of features)

	//-- Weights F1 -> F2
	float *W;				//- Synaptic weights : F1 -> F2
	float *normW;			//- Facteur |Wij|
	float *normWAlpha;		//- Facteur 1/(alpha+|Wij|)

	//-- F2 layer
	int	   sizeF2;			//- Nb. of committed neurons
	int   *sizeClasses;		//- Nb. of committed neurons for each class
	float *fMatch;			//- Match functions
	float *fChoice;			//- Choice functions
	int   *activeMt;		//- Boolean activation vector (for match tr.)
	int   *activeVl;		//- Boolean activation vector (for vigilance)

	//-- Weights  F2 -> Fab
	int *Wab;           	//- Mapfield : F2 -> Fab

	/*------------------------------ Functions -------------------------------*/
	//-- Constructors & destructor
	artmap();
	~artmap();

	//-- Functions - initialization
	void expand(int sizeAdded);
	void initialization(int nCls, int nFtr, int szF2max, int nEpoMax);
	void reinit();

	//-- Function - utilities
	void copy(artmap* nnOriginal);
	void save(char* nameFile, int currRep, int currBlk);
    void load(char* nameFile, configuration* cfg, int currRep, int currBlock);
	void setHp(float a, float b, float e, float r);

	//-- Function - training algorithm and changing pattern presentation order
	void train(dbase* db);
	int  prediction(dbase* db, int n);
	void orderRandomization(int nPatterns);
};
#endif
/*----------------------------------------------------------------------------*/
