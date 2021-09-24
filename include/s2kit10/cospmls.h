/***************************************************************************
  **************************************************************************
  
                           S2kit 1.0

          A lite version of Spherical Harmonic Transform Kit

   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
   Copyright 2004 Peter Kostelec, Dan Rockmore

   This file is part of S2kit.

   S2kit is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   S2kit is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with S2kit; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


#ifndef _COSPMLS_H
#define _COSPMLS_H

#if defined __cplusplus
  /// Public APIs will have C linkage if __cplusplus is present.
extern "C" {
#endif

int TableSize( int ,
		      int ) ;

int Spharmonic_TableSize( int ) ;

int Reduced_SpharmonicTableSize( int ,
					int ) ;

int Reduced_Naive_TableSize( int ,
				    int ) ;

int NewTableOffset( int ,
			   int ) ;

void CosPmlTableGen( int ,
			    int ,
			    double * ,
			    double * ) ;

int RowSize( int ,
		    int ) ;

int Transpose_RowSize( int ,
			      int ,
			      int ) ;

void Transpose_CosPmlTableGen( int ,
				      int ,
				      double * ,
				      double * ) ;

double **Spharmonic_Pml_Table( int ,
				      double * ,
				      double * ) ;

double **Transpose_Spharmonic_Pml_Table( double ** ,
						int ,
						double * ,
						double * ) ;

double **SemiNaive_Naive_Pml_Table( int ,
					   int ,
					   double * ,
					   double * ) ;

double **Transpose_SemiNaive_Naive_Pml_Table( double ** , 
						     int ,
						     int ,
						     double * ,
						     double * ) ;

#if defined __cplusplus
} /// End of C linkage section.
#endif

#endif

