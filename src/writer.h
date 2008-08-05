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

/**\file writer.h
 * \todo implement writer for EEGlab set-files
 */
#ifndef WRITER_H
# define WRITER_H
#include "definitions.h"
#include "helper.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
/* -------------------------------------------------------------------- */
/* functions */
/* -------------------------------------------------------------------- */

  void write_double_matrix_ascii(FILE *out, const double **d, int xdim, int ydim);
  void write_double_matrix_ascii_file(const char *fname, const double **d, int xdim, int ydim);

#ifdef __cplusplus
}
#endif

#endif
