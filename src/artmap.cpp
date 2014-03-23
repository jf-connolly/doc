/******************************************************************************
 *  Copyright(C) 2009 by Jean-Francois Connolly - jfconnolly@livia.etsmtl.ca  *
 *																			  *
 *  artmap.c - ARTMAP neural network (only fuzzy ARTMAP for now).             *
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
#include <fstream>
#include <math.h>
#include <string.h>

//-- Project specific
#include "artmap.h"

/*============================================================================*
 *  Function name			:	artmap::artmap
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Default constructor
 *============================================================================*/
artmap::artmap()
{
	error = 100;

	//-- General
	nClasses    	 = 0;
	nFeatures   	 = 0;
	nEpochMax   	 = 0;
	typeOp			 = opTrain;
	nEpochs	    	 = 0;
	nPatternsLearned = 0;

	//-- Hyperparameters
	alpha   = defaultAlpha;
	beta    = defaultBeta;
	epsilon	= defaultEpsilon;
	rho		= defaultRho;

	//-- Layer sizes
	sizeF0    = 0;
	sizeF1    = 0;
	sizeF2    = 0;
	sizeF2max = 0;

	//-- Weights F1 -> F2
	W          = (float*)NULL;
	normW      = (float*)NULL;
	normWAlpha = (float*)NULL;

	//-- F2 layer
	sizeF2      = 0;
	sizeClasses = (int*)NULL;
	fMatch      = (float*)NULL;
	fChoice     = (float*)NULL;
	activeMt    = (int*)NULL;
	activeVl    = (int*)NULL;

	//-- Weights F2 -> Fab
	Wab = (int*)NULL;
}
/*============================================================================*
 *  Function name			:	artmap::initialization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Initialization of an artmap neural network
 *----------------------------------------------------------------------------*
 *  Argument				:   int nCls        - Number of classes
 *                              int nFtr        - Number of features
 *                              int sizeF2Init	- Initial number of uncommited
 * 												  F2 nodes
 *                              int typeTr      - Type of training
 *============================================================================*/
void artmap::initialization(int nCls, int nFtr, int szF2max, int nEpoMax)
{
	//-- General
	nClasses  = nCls;
	nFeatures = nFtr;
	nEpochMax = nEpoMax;

	//-- Patterns presentation order
	sizeOrder = szF2max;
	order     = new int[sizeOrder];
	if(order == NULL)   { error += 6; return; }

	for(int i=0; i<sizeOrder; i++)   { order[i]=i; }

	//-- Layer sizes
	sizeF0    = nFeatures;
	sizeF1    = 2*nFeatures;
	sizeF2max = 0;

	//-- Class sizes
	sizeClasses = new int[nClasses];
	for(int k=0; k<nClasses; k++)   { sizeClasses[k] = 0; };
	if(sizeClasses == NULL)   { error += 6; return; }

	//-- Expand the network (adds uncommited nodes)
	expand(szF2max);

	if(error) { return; }
}
/*============================================================================*
 *  Function name			:	artmap::expand
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Adds F2 neurons and synaptic weights
 *----------------------------------------------------------------------------*
 *  Argument				:   int sizeAdded   - Number of F2 nodes added
 *============================================================================*/
void artmap::expand(int sizeAdded)
{
	int sizeF2Old  = sizeF2max;		  //-- Old F2 size
	sizeF2max     += sizeAdded;   //-- New F2 layer size

	/*----------------- Memory allocation (1D): Info & norms -----------------*/
	//-- Match function , choice function & norms
	fMatch = new float[sizeF2max];
	if(fMatch == NULL)     { error += 6; return; }

	fChoice = new float[sizeF2max];
	if(fChoice == NULL)    { error += 6; return; }

	normW = new float[sizeF2max];
	if(normW == NULL)      { error += 6; return; }

	normWAlpha = new float[sizeF2max];
	if(normWAlpha == NULL) { error += 6; return; }

	//-- Active node boolean vectors (vigilance & match tracking)
	activeMt = new int[sizeF2max];
	if(activeMt == NULL)   { error += 6; return; }

	activeVl = new int[sizeF2max];
	if(activeVl == NULL)   { error += 6; return; }

	/*------------------- Memory allocation (2D): Weights --------------------*/
	//-- Wij (element variable: nF2_Size)
	W = new float[sizeF2max*sizeF1];
	if(W == NULL)          { error += 6; return; }

	//-- Mapfield
	Wab = new int[sizeF2max*nClasses];
	if(Wab == NULL)        { error += 6; return; }

	/*---------------------------- Initialization ----------------------------*/
	float normAlpha = 1/(sizeF1+alpha);

	for(int j=sizeF2Old; j<sizeF2max; j++)
	{
		//-- Match & choice function; |Wij| & 1/(alpha+|Wij|)
		fMatch    [j] = 0;
		fChoice   [j] = 0;
		normW     [j] = sizeF1;
		normWAlpha[j] = normAlpha;

		//-- F1 -> F2 weights (initialized at 1)
		for(int i=0; i<sizeF1; i++)   { *(W + j*sizeF1 +i) = 1; }

		//-- F2 -> Fab weights (initialized at 0)
		for(int k=0; k<nClasses; k++) { *(Wab + j*nClasses +k) = 0; }
	}
}
/*============================================================================*
 *  Function name			:	artmap::reinit
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	NN weights & norms re-initialization
 *============================================================================*/
void artmap::reinit()
{
	/*------------------------ Weights initialization ------------------------*/
	float normAlpha = 1/(sizeF1+alpha);

	//-- Only the commited nodes are reinitialized
	for(int j=0; j<sizeF2; j++)
	{
		//-- Match & choice function; |Wij| & 1/(alpha+|Wij|)
		fMatch    [j] = 0;
		fChoice   [j] = 0;
		normW     [j] = sizeF1;
		normWAlpha[j] = normAlpha;

		//-- F1 -> F2 weights (initialized at 1)
		for(int i=0; i<sizeF1; i++)   { *(W   + j*sizeF1   +i) = 1; }

		//-- F2 -> Fab weights (initialized at 0)
		for(int k=0; k<nClasses; k++) { *(Wab + j*nClasses +k) = 0; }
	}

	//-- F2 layer
	for(int k=0; k<nClasses; k++)   { sizeClasses[k] = 0; }

	/*----------------------- Constants initialization -----------------------*/
	//-- Default values
	error   		 = 100;
	nEpochs 		 = 0;
	nPatternsLearned = 0;

	//-- Hyperparameters
	alpha   = defaultAlpha;
	beta    = defaultBeta;
	epsilon	= defaultEpsilon;
	rho		= defaultRho;

	//-- F2 layer
	sizeF2 = 0;
}
/*============================================================================*
 *  Function name			:	artmap::artmap
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Destructor
 *----------------------------------------------------------------------------*
 *  Argument				:   int nCls        - Number of classes
 *                              int nFtr        - Number of features
 *                              int sizeF2Init	- Initial number of uncommited
 * 												  F2 nodes
 *                              int typeTr      - Type of training
 *============================================================================*/
artmap::~artmap()
{
	//-- F1 -> F2 weights & norms
	delete [] W;
	delete [] normW;
	delete [] normWAlpha;
	delete [] activeMt;
	delete [] activeVl;

	//-- F2 layer
	delete [] fMatch;
	delete [] fChoice;
	delete [] sizeClasses;

	//-- F2 -> Fab weights
	delete [] Wab;

	/*------------------------------ Constantes ------------------------------*/
	//-- General
	error   		 = 100;
	nEpochs 		 = 0;
	nPatternsLearned = 0;

	//-- Hyperparametres
	alpha	= 0;
	beta	= 0;
	epsilon	= 0;
	rho		= 0;

	//-- Couche F0
	sizeF0 = 0;

	//-- Couche F1
	sizeF1 = 0;

	//-- Couche F2
	sizeF2 = 0;
}
/*============================================================================*
 *  Function name			:	artmap::copy
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Copies the artmap neural network
 *----------------------------------------------------------------------------*
 *  Argument				:   artmap *nnOriginal  - Original neural network
 *============================================================================*/
void artmap::copy(artmap* nnOriginal)
{
	size_t sizeData;

	//-- Neural networks existence
	if(nnOriginal == NULL)   { error += 2; return; }

	//-- General
	nEpochs 		 = nnOriginal->nEpochs;
	nPatternsLearned = nnOriginal->nPatternsLearned;

	//-- Hyperparametres
	alpha   = nnOriginal->alpha;
	beta    = nnOriginal->beta;
	epsilon = nnOriginal->epsilon;
	rho     = nnOriginal->rho;

	//-- Sizes
	nFeatures = nnOriginal->sizeF0;
	sizeF1    = nnOriginal->sizeF1;
	sizeF2    = nnOriginal->sizeF2;

	sizeData =  sizeof(int) * (unsigned)(nClasses);
	memcpy(sizeClasses, nnOriginal->sizeClasses, sizeData);

	/*-------------------- Copy - Layers, weight & norms ---------------------*/
	//-- Wji
	sizeData =  sizeof(float) * (unsigned)(sizeF2max * sizeF1);
	memcpy(W, nnOriginal->W, sizeData);

	//-- |Wji| & 1 / (|Wji|+alpha)   
	sizeData =  sizeof(float) * (unsigned)(sizeF2max);
	memcpy(normW,      nnOriginal->normW,      sizeData);
	memcpy(normWAlpha, nnOriginal->normWAlpha, sizeData);

	//-- Mapfield
	sizeData = sizeof(int) * (unsigned)(sizeF2max * nClasses);
	memcpy(Wab, nnOriginal->Wab, sizeData);

	//-- Check
	if(W          == NULL) { error += 6; return; }
	if(normW      == NULL) { error += 6; return; }
	if(normWAlpha == NULL) { error += 6; return; }
	if(Wab        == NULL) { error += 6; return; }
}
/*============================================================================*
 *  Function name			:	artmap::save
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Copies the artmap neural network
 *----------------------------------------------------------------------------*
 *  Argument				:   artmap *nn       - Neural network saved
 * 								char nameFile    - File to save
 *============================================================================*/
void artmap::save(char* nameFile, int currRep, int currBlk)
{
	//-- File opening
	ofstream f(nameFile, ios::app | ios::binary);
	if(!f)   { printf("Error! Can't open the file %s \n", nameFile);  exit(1); }

	//-- General informations
	f.write((char*)&currRep,     sizeof(int)                       );
	f.write((char*)&currBlk,     sizeof(int)                       );
	f.write((char*)&nClasses,    sizeof(int)                       );
	f.write((char*)&nFeatures,   sizeof(int)                       );
	f.write((char*)&sizeF2,      sizeof(int)                       );
	f.write((char*) sizeClasses, sizeof(int) * (unsigned)(nClasses));

	//-- Hyperparameters
	f.write((char*)&alpha,   sizeof(float));
	f.write((char*)&beta,    sizeof(float));
	f.write((char*)&epsilon, sizeof(float));
	f.write((char*)&rho,     sizeof(float));

	//-- Wij (F2 node after F2 node)
	f.write((char*)W,          sizeof(float) * (unsigned)(sizeF2 * sizeF1));
	f.write((char*)normW,      sizeof(float) * (unsigned)(sizeF2)         );
	f.write((char*)normWAlpha, sizeof(float) * (unsigned)(sizeF2)         );

	//-- Map field
	f.write((char*)Wab, sizeof(int) * (unsigned)(sizeF2 * nClasses));

	//-- File closing
	f.close();
}
/*============================================================================*
 *  Function name			:	artmap::load
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Copies the artmap neural network
 *----------------------------------------------------------------------------*
 *  Argument				:   artmap *nn      - Neural network saved
 * 								char nameFile   - Loaded file
 *  Returns					:   int position	- Get pointer's position
 *============================================================================*/
void artmap::load(char* nameFile, configuration* cfg, int currRep,int currBlk)
{
	int savedRep, savedBlk;
	streampos currentPos, offset;

	//-- File opening
	ifstream f(nameFile, ios::binary);
	if(!f)   { printf("Error! Can't open the file %s \n", nameFile);  exit(1); }

	//-- General informations
	int nBlocks       = cfg->nBlocks;
	int nReplications = cfg->nReplications;
	int nRead         = nReplications*nBlocks;

	/*----------------- Neural networks that are passed over -----------------*/
	for(int addr=0; addr<nRead; addr++)
	{
		//-- General informations
		f.read( (char*)&savedRep,  sizeof(int) );
		f.read( (char*)&savedBlk,  sizeof(int) );
		f.read( (char*)&nClasses,  sizeof(int) );
		f.read( (char*)&nFeatures, sizeof(int) );
		f.read( (char*)&sizeF2,    sizeof(int) );

		if( savedBlk!=currBlk || savedRep!=currRep )
		{
			//-- ... jump to the next neural network depending of the info.
			offset = (streampos)( sizeof(float) *( 4 + sizeF2 * (sizeF1+2)   ) +
							      sizeof(int)   *( sizeF2*nClasses +nClasses ) );

			currentPos = f.tellg() + offset;
			f.seekg(currentPos);
		}
		/*--------------------- The right neural network ---------------------*/
		else
		{	//-- Class sizes
			f.read( (char*)sizeClasses, sizeof(int)*(unsigned)(nClasses) );

			//-- Hyperparameters
			f.read( (char*)&alpha,   sizeof(float) );
			f.read( (char*)&beta,    sizeof(float) );
			f.read( (char*)&epsilon, sizeof(float) );
			f.read( (char*)&rho,     sizeof(float) );

			//-- Wij (F2 node after F2 node)
			f.read( (char*)W,          sizeof(float)*(unsigned)(sizeF2*sizeF1) );
			f.read( (char*)normW,      sizeof(float)*(unsigned)(sizeF2) );
			f.read( (char*)normWAlpha, sizeof(float)*(unsigned)(sizeF2) );

			//-- Map field
			f.read( (char*)Wab, sizeof(int)*(unsigned)(sizeF2 * nClasses) );

			break;
		}
	}

	//-- File closing
	f.close();
}
/*============================================================================*
 *  Function name			:	artmap::setHp
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Sets hyperparameters
 *----------------------------------------------------------------------------*
 *  Argument				:   float a, b, e,& r   - Hyperparameters
 *============================================================================*/
void artmap::setHp (float a, float b, float e,	float r)
{
	//-- Validation
	if(a<=0)		{ error += 1; return; }
	if(b<0  || b>1)	{ error += 1; return; }
	if(e<-1 || e>1)	{ error += 1; return; }
	if(r<0  || r>1)	{ error += 1; return; }

	//-- Adjusment
	alpha	= a;
	beta	= b;
	epsilon	= e;
	rho		= r;

	//-- Norm reevaluation 1/(alpha+|wj|)
	for(int j=0; j<sizeF2; j++)  { normWAlpha[j] = 1.0 / (normW[j] + alpha); }
}
/*============================================================================*
 *  Function name			:	artmap::train
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Fuzzy artmap learning algorithm
 *----------------------------------------------------------------------------*
 * Argument                 :   dbase  *db     - Data base
 *============================================================================*/
void artmap::train(dbase* db)
{
	int 	addr;                  //-- Address

	int		newNode;		       //-- Is there a new node
	int		jStar;			       //-- Winning F2 node
	int		reset;				   //-- reset signal
	int		mTracking;			   //-- Match Tracking signal
	float	dRho;				   //-- Internal vigilance parameter
	float 	Wold, dMin;

	float  *data;                  //-- Current pattern
	int		label;				   //-- Current label

	int    nFailed, nDeactivated;  //-- Number of failed / deactivated F2 nodes
	float  tj, tChoosen;		   //-- Tj in the learning algorithm

	//-- Data base existence
	if(db == NULL)                   { error += 2; return; }
	if(db->nClasses  != nClasses )   { error += 3; return; }
	if(db->nFeatures != nFeatures)   { error += 3; return; }

	/*------------------------------------------------------------------------*/
	/*----------------- Main loop (for each training samples) ----------------*/
	/*------------------------------------------------------------------------*/
	//-- Number of training patterns and randomization
	int nPatterns = db->nPatterns;
	orderRandomization(nPatterns);

//	printf("hp: (%1.2f, %1.2f, %1.2f, %1.2f)\n", alpha, beta, epsilon, rho);

	for(int p=0; p<nPatterns; p++)
	{
		//-- Load current sample
		data  = db->getPatternAddr(order[p]);
		label = db->labels[ order[p] ];

		//-- (Re) initialiser rho
		dRho = rho;
 		if(dRho > 0.999) dRho = 0.999;
		
		//************** For MoBo only **************//
		if(dRho > 0.9) dRho = 0.9;		
		//************** For MoBo only **************//

//		printf("Pattern %d, sizeF2: %d \n", p, sizeF2);

		/*--------------------- Prototype selection (1/2) --------------------*/
		//-- Choice and match function evaluation, & prototype selection
//		printf("T(j):");
		jStar = 0;
		for(int j=0; j<sizeF2; j++)
		{
			//-- L^1 distance between sample A and node j = |min(A,Wij)|
			float normAW = 0.0;

			//-- BOTTLE NECK!!!
			for(int i=0; i<sizeF1; i++)
			{	normAW += fmin(data[i], W[ j*sizeF1 +i ]);	}

			fMatch [j]  = normAW / (float)nFeatures; // |min(A,Wij)|/M
			fChoice[j] = normAW * normWAlpha[j];     // |min(A,Wij)|/(a+|Wij|)

			//-- Prototype selection (jStar of all nodes)
			jStar = argMax(fChoice[jStar], fChoice[j], jStar, j);

//			printf(" (%1.2f, %1.2f)", fMatch[j], fChoice[j]);
		}
//		printf("\nWinner: %d, winner's class %d, true class: %d", jStar, );

		/*------------- Activate all neurons and initialize loop -------------*/
		reset        = 1;
		nFailed      = 0;
		nDeactivated = 0;

		for(int j=0; j<sizeF2; j++)
		{	activeMt[j] = 1;   activeVl[j] = 1;   }

		/*------------ Prototype selection (2/2) & match tracking ------------*/
		while(reset)
		{
			reset     = 0;		//-- If the vigilance test is passed and no
								//   match tracking -> to learning
			mTracking = 0;

			//-- If no node passes the tests (at the beginning for nF2 = 0)
			if (nFailed == sizeF2)   { reset = 0;   break; }

			/*------------------------ Vigilance test ------------------------*/
			//-- If failed: reset the node.
			if (fMatch[jStar] < dRho)   { reset = 1; }
			else//-- If pass, check if the prediction is correct.
			{
				//-- If the prediction is incorrect: reset & match tracking
				//-- (We look at pnMapfield[nChoice][nLabel].)
				if (Wab[ jStar*nClasses + label ] == 0)
				{   reset = 1;   mTracking = 1;   }
				else//-- If all tests are passed, we scream our joy & learn
				{	reset = 0;   }
			}

			/*----------------------------- Reset ----------------------------*/
			if (reset)
			{	//-- Reset the node (for vigilance),
				activeVl[jStar] = 0;

				//-- count number of failure,
				nFailed ++;

				if (mTracking)
				{	//-- reset permanently (for this pattern),
					activeMt[jStar] = 0;

					//-- adjust rho (if match tracking),
					dRho = fMatch[jStar] + epsilon;
					if(dRho > 0.99)  { dRho = 0.99; }//-- if rho > 0.999
					if(dRho < 0)     { dRho = 0;   }//-- if rho < 0

					//-- reinitialize vigilance activation vector
					for(int j=0; j<sizeF2; j++)
					{	if (activeMt[j]) { activeVl[j] = 1; }   }

					//-- count the number of deactivation,
					nDeactivated ++;
					nFailed = nDeactivated;
				}

				//-- and choose the next jStar Tj.
				for(int j=0; j<sizeF2; j++)
				{	tChoosen = fChoice[jStar] * activeVl[jStar];
					tj       = fChoice[j]     * activeVl[j];
					jStar    = argMax(tChoosen, tj, jStar, j);
				}
			}
		}

		/*------------ If no neuron satisfied the vigilance test, ------------*/
		/*-------------------- a commited neuron is added --------------------*/
		newNode = 0;
		if (nFailed == sizeF2)
		{	//-- It's new
			newNode = 1;

			//-- The node after the last is chosen
			jStar = sizeF2;

			//-- Learn class pnMapfield[nChoice][nLabel] = 1;
			Wab[ jStar*nClasses + label ] = 1;

			//-- Adjust the number of commited node
			sizeF2             += 1;
			sizeClasses[label] += 1;

			//-- Add uncommited neurons (NOT implemented!!!)
			if (sizeF2 == sizeF2max)
			{	printf("Nb de Jean-Paul : %d\n", sizeF2);
				error += 66; exit(0);
			}
		}

		/*----------------------------- Learning -----------------------------*/
		normW[jStar] = 0.0;

		//-- For new nodes only
		if(newNode)
		{
			for(int i=0; i<sizeF1; i++)
			{
				//-- The new box is a copy of the pattern
				addr          = jStar*sizeF1 +i;
				W[addr]       = data[i];
				normW[jStar] += W[addr];
			}
		}
		//-- Usual training
		else
		{
			for(int i=0; i<sizeF1; i++)
			{
				addr  = jStar*sizeF1 +i;

				//-- Wij nn -> pdWji[nChoice][i]
				Wold = W[addr];
				dMin = fmin(Wold, data[i]);

				//-- w_new = beta*min(I,w_old) + (1-beta)*w_old
				W[addr] = (beta)*dMin + (1 - beta)*Wold;

				//-- |Wij|[nChoice][i]
				normW[jStar]+= W[addr];
			}
		}

		//- Evaluate 1/(Alpha+|Wij|)
		normWAlpha[jStar] = 1.0 / (alpha + normW[jStar]);

//		printf("Pattern %d, sizeF2: %d \n", p, sizeF2);
//		getchar();
	}

	nEpochs++;
}
/*============================================================================*
 *  Function name			:	artmap::prediction
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Makes a prediction
 *----------------------------------------------------------------------------*
 *  Arguments :    dbase *db      - Data base
 *                 int    n       - index of the tested pattern
 *  Returns :      int    kStar   - Prediction (see notation)
 *============================================================================*/
int artmap::prediction(dbase* db, int p)
{
	/*------------------------- Parameters validation ------------------------*/
	//-- Data base
	if(db == NULL)                   { error += 2; return -1; }
	if(db->nClasses  != nClasses )   { error += 3; return -1; }
	if(db->nFeatures != nFeatures)   { error += 3; return -1; }

	//-- n
	if( p < 0 || p > db->nPatterns ) { error += 4; return -1; }

	/*--------------------- Prototype selection (1/2) --------------------*/
	//-- Load current sample
	float *data  = db->getPatternAddr(p);   //- Current pattern

	//-- Choice and match function evaluation, & prototype selection
	int jStar = 0;   //-- jStar F2 node
	for(int j=0; j<sizeF2; j++)
	{
		//-- L^1 distance between sample A and node j = |min(A,Wij)|
		float normAW = 0.0;

		//-- |min(A,Wij)|/(a+|Wij|)  -  BOTTLE NECK!!!
		for(int i=0; i<sizeF1; i++)
		{	normAW += fmin(data[i], W[ j*sizeF1 +i ]);	}

		fChoice[j] = normAW * normWAlpha[j];

		//-- Prototype selection (jStar of all nodes)
		jStar = argMax(fChoice[jStar], fChoice[j], jStar, j);
	}

	/*------------------------ Prediction evaluation -------------------------*/
	//-- Making the prediction
	int	kStar = 0;			       //- Predicted class
	for(int k=0; k<nClasses; k++)
	{
		kStar = argMax( Wab[ jStar*nClasses + kStar ],  Wab[ jStar*nClasses + k ],
						kStar,                          k);
	}

	return kStar;
}
/*============================================================================*
 *  Function name			:	artmap::orderRandomization
 *  Originally written by	:	Jean-Francois Connolly
 *  Description				:	Changes the presentation order
 *============================================================================*/
void artmap::orderRandomization(int nPatterns)
{
	int	   nRnd, nTmp;
	float  dNormalized;

	//-- Reinitialization
	for(int p=0; p<nPatterns; p++)   { order[p] = p; }

	//-- Randomization
	for(int p=nPatterns-1; p>0; p--)
	{
		//-- Real in [0,1[
		dNormalized = ((float)rand() / ((float)(RAND_MAX)+(float)(1)) );

		//-- Integer in {0,1, ... , i-1}
		nRnd = (int)(dNormalized * p);

		//-- Permutation
		nTmp        = order[nRnd];
		order[nRnd] = order[p];
		order[p]    = nTmp;
	}
}
