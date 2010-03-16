/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke                                  *
 *   mihrke@uni-goettingen.de                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**\file som.h
 \brief \ref status_experimental Self-Organizing Maps

 \warning For now this is only experimental but plans are to complete this file.


\section somalg SOM Algorithm
Self-Organizing Feature Maps are extensively described in 

Kohonen, T. 1995. Self-Organizing Maps. Springer.

We have a set of n neurons \f$ W = \{ \vec{w}_1,\ldots, \vec{w}_n\} \f$ that
are of dimension m: \f$ \vec{w}_i \in \mathcal{R}^m \f$. The neurons
are arranged on a bidirectional, connected graph, such that there
exists a distance between two neurons \f$ d(\vec{w}_i,\vec{w}_j) \f$.

Using input data \f$ \vec{x}_1,\ldots, \vec{x}_n\f$, of dimension m, 
we find the best-matching unit (BMU) 
\f[
w_{\mbox{bmu}} = \mbox{argmin}_{w\in W} || x - w ||^2.
\]
(we can use any other metric as well - many metrices are implemented
in LibEEGTools).

Then, we update each neuron by the learning rule
\f[
w_i(t+1) = w_i(t) + \alpha(t)h(t)[ x - m_i(t) ]
\f]
where \f$\alpha\f$ is a time-decaying function, h(t) is a 
neighbourhood function centered on the BMU (e.g. gaussian) which 
narrows down with increasing time.

In this implementation, we divide the training phase into an initial 
and a fine-tuning phase. During the initial phase, $\alpha$ takes on rather 
large values and the neighbourhood is wide such that an initial 
arrangement may take place. Later, both parameters are narrowed down to enable
fine-tuning.

All of the parameters can be conveniently manipulated in the implementation.

\section somimpl SOM Implementation 

The implementation given here is quite general but not very efficient. 
Topology of the network of neurons is implemented in terms of distance
matrices between two neurons (in terms of their ''physical'' distance D
on the graph) \f$ \mathbf{M} = m_{ij} = D(i,j)\f$.

The user can either supply this matrix directly or generate such a matrix
using our convenience function som_generate_connectivity_matrix(). 


\section somusage Usage

First, you need to initialize a Som-struct using som_init(). You need to 
supply a value of SOMConnectivityType. If it is not CUSTOM, the function will 
allocate a connectivity matrix.

\code
Som *S;
S = som_init( 2, 1000, 1e6, ONED_LINEAR );
\endcode

\section sompython Python Usage


 \todo generalize metric (need to solve Eq. 3.40 (from Kohonen, 1995) for each metric)
 */
#ifndef SOM_H
# define SOM_H

#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  /**\ingroup som_help 
	\{ */
  Som* som_init( int dimension, int n, int nruns, SOMConnectivityType connectivity );
  void som_free( Som *s );
  void som_print( FILE *f, Som *s );
  /** \} */

  /**\ingroup som_initialize
	\{ */ 
  void som_initialize_random( Som *s, double min, double max );
  void som_initialize_random_samples( Som *s, double **X, int dim, int nsamples );
  double** som_generate_connectivity_matrix( SOMConnectivityType type, double **m, int n );
  /** \} */

  /**\ingroup som_neighbourhood
	\{ */
  double som_neighbourhood_gaussian( int x, int bmu, struct som_struct *s, int t);
  /** \} */

  double som_time_decay_linear( int t, int nruns, int initial_runs );

  /**\ingroup som_train
	\{ */
  void som_train_from_data( Som *s, double **X, int dim, int nsamples );
  /** \} */
  
#ifdef __cplusplus
}
#endif

#endif /* SOM_H */
