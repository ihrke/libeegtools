/***************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke   *
 *   mihrke@uni-goettingen.de   *
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

/**\file tools.h
 * Other Tools.
 * \defgroup tools Some Other Tools
 *\{
 *\}
 */
#ifndef TOOLS_H
#define TOOLS_H
#include "mathadd.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \addtogroup tools
 *\{
 */
    void eeg_remove_baseline( EEGdata *eeg, const double *times, double win_from, double win_to );
    void eegtrials_remove_baseline( EEGdata_trials *eeg, double win_from, double win_to );
/** \} */

#ifdef __cplusplus
}
#endif

#endif
