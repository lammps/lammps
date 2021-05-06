/*----------------------------------------------------------------------
xfI  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "reaxc_types.h"

#include <stddef.h>

#include "init_md.h"
#include "list.h"
#include "random.h"
#include "reset_tools.h"
#include "tool_box.h"
#include "vector.h"





/************************ initialize workspace ************************/
/* Initialize Taper params */
void Init_Host_Taper( control_params * const control,  storage * const workspace )
{
    real d1, d7;
    const real swa = control->nonb_low;
    const real swb = control->nonb_cut;
    real swa2, swa3;
    real swb2, swb3;

    if ( FABS( swa ) > 0.01 )
    {
        fprintf( stderr, "[WARNING] non-zero lower Taper-radius cutoff in force field parameters\n" );
    }

    if ( swb < 0.0 )
    {
        fprintf( stderr, "[ERROR] negative upper Taper-radius cutoff in force field parameters\n" );
        MPI_Abort( MPI_COMM_WORLD,  INVALID_INPUT );
    }
    else if ( swb < 5.0 )
    {
        fprintf( stderr, "[WARNING] very low Taper-radius cutoff in force field parameters (%f)\n", swb );
    }

    d1 = swb - swa;
    d7 = POW( d1, 7.0 );
    swa2 = SQR( swa );
    swa3 = CUBE( swa );
    swb2 = SQR( swb );
    swb3 = CUBE( swb );

    workspace->Tap[7] =  20.0 / d7;
    workspace->Tap[6] = -70.0 * (swa + swb) / d7;
    workspace->Tap[5] =  84.0 * (swa2 + 3.0 * swa * swb + swb2) / d7;
    workspace->Tap[4] = -35.0 * (swa3 + 9.0 * swa2 * swb + 9.0 * swa * swb2 + swb3 ) / d7;
    workspace->Tap[3] = 140.0 * (swa3 * swb + 3.0 * swa2 * swb2 + swa * swb3 ) / d7;
    workspace->Tap[2] = -210.0 * (swa3 * swb2 + swa2 * swb3) / d7;
    workspace->Tap[1] = 140.0 * swa3 * swb3 / d7;
    workspace->Tap[0] = (-35.0 * swa3 * swb2 * swb2 + 21.0 * swa2 * swb3 * swb2
            + 7.0 * swa * swb3 * swb3 + swb3 * swb3 * swb ) / d7;
}








