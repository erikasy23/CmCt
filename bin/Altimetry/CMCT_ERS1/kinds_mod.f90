MODULE Kinds_MOD

!  NAME: kinds_mod
!
!  PURPOSE: definition of the symbolic kind names used in TYPE statements
!
!  COMMENTS: All parameters defined in this module are default-type integers.
!            The GLAS naming convention suffix (e.g., _i4) is omitted from
!            these parameters in order to simplify type declarations in other
!            software components.
!
!            GLAS parameter names for integer kinds:
!               i1b = 1-byte integer
!               i2b = 2-byte integer
!               i4b = 4-byte integer
!
!            GLAS parameter names for logical kinds:
!               l1b = 1-byte logical
!               l4b = 4-byte logical
!
!            GLAS parameter names for real kinds:
!               r4b = 4-byte (single precision) real
!               r8b = 8-byte (double precision) real
!
!            GLAS parameter names for complex kinds:
!               c4b = 4-byte (single precision) complex
!               c8b = 8-byte (double precision) complex
!
!            Integer ranges:                minimum value       maximum value
!
!                  1-byte (i1b)                      -126                +127
!                  2-byte (i2b)                    -32768              +32767
!                  4-byte (i4b)               -2147483648         +2147483647
!
!            Real/complex magnitudes:  smallest magnitude   largest magnitude
!
!                  4-byte (r4b,c4b)         1.175494E-38        3.402823E+38
!                  8-byte (r8b,c8b)         2.225074D-308       1.797693D+308
!
!  TO COMPILE: f90 -c kinds_mod.f90
!
!  HISTORY: 11/09/98 J. McMillan/Raytheon ITSS - initial version
!           02/19/99 K. Barbieri/Raytheon ITSS - changed name kindsmod to
!                                                   kinds_mod
!           02/22/99 LA. Roberts/Raytheon ITSS - added kinds l1b, l4b
!
!-------------------------------- end of prolog --------------------------------

   implicit none
   private

!---------------------- parameter names for integer kinds ----------------------

   integer, public, parameter :: i1b = selected_int_kind(1)
   integer, public, parameter :: i2b = selected_int_kind(4)
   integer, public, parameter :: i4b = selected_int_kind(9)
   integer, public, parameter :: i8b = selected_int_kind(15)

!---------------------- parameter names for logical kinds ----------------------

   integer, public, parameter :: l1b = KIND(.TRUE._1)
   integer, public, parameter :: l4b = KIND(.TRUE._4)

!---------------------- parameter names for real kinds -------------------------

   integer, public, parameter :: r4b = selected_real_kind( 6, 37)
   integer, public, parameter :: r8b = selected_real_kind(15,307)

!---------------------- parameter names for complex kinds ----------------------

   integer, public, parameter :: c4b = kind((1.0_r4b,1.0_r4b))
   integer, public, parameter :: c8b = kind((1.0_r8b,1.0_r8b))

END MODULE Kinds_MOD

