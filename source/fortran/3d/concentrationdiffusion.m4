c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c 
define(NDIM,3)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
c
c Coefficient [h(phi)*d_solid+(1-h(phi))*d_liquid]
c
      subroutine concentration_pfmdiffusion(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diff0, diff1, diff2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
      double precision gas_constant_R
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision vphi, hphi, invT
      double precision q0_liquid_invR, q0_solidA_invR
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0-1,ic1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diff0(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1-1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diff1(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1,ic2-1) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diff2(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid +
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      return
      end
c
c same as function concentrationdiffusion0, without accumulating
c component into single D
c
      subroutine concentration_pfmdiffusion_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diffL0, diffL1, diffL2, 
     &   diffA0, diffA1, diffA2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
      double precision gas_constant_R
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision vphi, hphi, invT
      double precision q0_liquid_invR, q0_solidA_invR
      double precision diff_liquid, diff_solidA
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
c
      do ic2 = ifirst2-ngdiff, ilast2+ngdiff
         do ic1 = ifirst1-ngdiff, ilast1+ngdiff
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0-1,ic1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diffL0(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid
               diffA0(ic0,ic1,ic2) =
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      do ic2 = ifirst2-ngdiff, ilast2+ngdiff
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0-ngdiff, ilast0+ngdiff

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1-1,ic2) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diffL1(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid 
               diffA1(ic0,ic1,ic2) =
     &             hphi * diff_solidA
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1-ngdiff, ilast1+ngdiff
            do ic0 = ifirst0-ngdiff, ilast0+ngdiff

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               invT = 2.0d0 /
     &            ( temp(ic0,ic1,ic2-1) + temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )

               diffL2(ic0,ic1,ic2) =
     &            ( 1.0d0 - hphi ) * diff_liquid
               diffA2(ic0,ic1,ic2) =
     &            hphi * diff_solidA

            end do
         end do
      end do
c
      return
      end
c
      subroutine concentration_pfmdiffusion_of_temperature_threephases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, nphi, nphia, ngphi,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2,
     &   diffB0, diffB1, diffB2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   d_solidB, q0_solidB,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphi, nphia, ngphi, ngdiff, ngtemp, three_phases
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solidA, d_solidB
      double precision q0_liquid, q0_solidA, q0_solidB
      double precision gas_constant_R
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi),nphi)
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffB0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffB1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffB2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2, ip
      double precision vphi, phiL, phiA, phiB, invT
      double precision q0_liquid_invR, q0_solidA_invR, q0_solidB_invR
      double precision diff_liquid, diff_solidA, diff_solidB
      double precision interp_func
      double precision average_func
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
      q0_solidB_invR = q0_solidB / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = 0.d0
               do ip = 1, nphia
                  vphi = vphi + average_func(
     &               phi(ic0-1,ic1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = nphia+1, nphi-1
                  vphi = vphi + average_func(
     &               phi(ic0-1,ic1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2,nphi), phi(ic0,ic1,ic2,nphi),
     &            avg_type )
               phiL = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0-1,ic1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
               diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

               diffL0(ic0,ic1,ic2) = phiL * diff_liquid
               diffA0(ic0,ic1,ic2) = phiA * diff_solidA
               diffB0(ic0,ic1,ic2) = phiB * diff_solidB

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphia
                  vphi = vphi + average_func(
     &               phi(ic0,ic1-1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = nphia+1, nphi-1
                  vphi = vphi + average_func(
     &                phi(ic0,ic1-1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &                avg_type )
               enddo
               phiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2,nphi), phi(ic0,ic1,ic2,nphi),
     &            avg_type )
               phiL = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1-1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
               diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

               diffL1(ic0,ic1,ic2) = phiL * diff_liquid
               diffA1(ic0,ic1,ic2) = phiA * diff_solidA
               diffB1(ic0,ic1,ic2) = phiB * diff_solidB

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphia
                  vphi = vphi + average_func(
     &               phi(ic0,ic1,ic2-1,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phiA = interp_func( vphi, interp_type )

               vphi = 0.d0
               do ip = nphia+1, nphi-1
                  vphi = average_func(
     &               phi(ic0,ic1,ic2-1,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phiB = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1,nphi), phi(ic0,ic1,ic2,nphi),
     &            avg_type )
               phiL = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1,ic2-1)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solidA = d_solidA * exp( -q0_solidA_invR * invT )
               diff_solidB = d_solidB * exp( -q0_solidB_invR * invT )

               diffL2(ic0,ic1,ic2) = phiL * diff_liquid
               diffA2(ic0,ic1,ic2) = phiA * diff_solidA
               diffB2(ic0,ic1,ic2) = phiB * diff_solidB

            end do
         end do
      end do

      return
      end
c
c Coefficient \tilde D from Beckermann, Diepers, Steinbach, Karma, Tong, 1999
c
      subroutine concentrationdiffusion_beckermann(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, ngphi,
     &   diff0, diff1, diff2, ngdiff,
     &   partition_coeff, ngk,
     &   d_liquid,
     &   d_solid,
     &   interp_type,
     &   avg_type )
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngphi, ngdiff, ngk
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid
c
c variables in 3d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi))
      double precision partition_coeff(CELL3d(ifirst,ilast,ngk))
      double precision diff0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diff1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diff2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************     
c
      integer ic0, ic1, ic2
      double precision vphi, hphi, k
      double precision diff_liquid, diff_solid
      double precision interp_func
      double precision average_func
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2), phi(ic0,ic1,ic2), avg_type )
               
               hphi = interp_func( vphi, interp_type )
               
               k = partition_coeff(ic0,ic1,ic2)

               diff0(ic0,ic1,ic2) = d_solid +
     &            (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                      /(1.0d0 - hphi+k*hphi)

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               k = partition_coeff(ic0,ic1,ic2)

               diff1(ic0,ic1,ic2) = d_solid +
     &            (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                      /(1.0d0 - hphi+k*hphi)

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1), phi(ic0,ic1,ic2), avg_type )

               hphi = interp_func( vphi, interp_type )

               k = partition_coeff(ic0,ic1,ic2)

               diff2(ic0,ic1,ic2) = d_solid +
     &            (d_liquid-d_solid)*( 1.0d0 - hphi )
     &                      /(1.0d0 - hphi+k*hphi)

            end do
         end do
      end do
c
      return
      end

      subroutine concentration_diffcoeff_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solidA, q0_solidA,
     &   gas_constant_R)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer ngdiff, ngtemp
      double precision d_liquid, d_solidA
      double precision q0_liquid, q0_solidA
      double precision gas_constant_R
c variables in 3d cell indexed
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
c
c***********************************************************************
c***********************************************************************
c
      integer ic0, ic1, ic2
      double precision invT
      double precision q0_liquid_invR, q0_solidA_invR
c
      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solidA_invR = q0_solidA / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1

               invT = 2.0d0 /
     &                ( temp(ic0-1,ic1, ic2) + temp(ic0,ic1,ic2) )

               diffL0(ic0,ic1,ic2) = d_liquid *
     &                               exp( -q0_liquid_invR * invT )
               diffA0(ic0,ic1,ic2) = d_solidA *
     &                               exp( -q0_solidA_invR * invT )

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               invT = 2.0d0 /
     &               ( temp(ic0,ic1-1,ic2) + temp(ic0,ic1,ic2) )

               diffL1(ic0,ic1,ic2) = d_liquid *
     &                               exp( -q0_liquid_invR * invT )
               diffA1(ic0,ic1,ic2) = d_solidA *
     &                               exp( -q0_solidA_invR * invT )

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               invT = 2.0d0 /
     &                ( temp(ic0,ic1,ic2-1) + temp(ic0,ic1,ic2) )

               diffL2(ic0,ic1,ic2) = d_liquid *
     &                               exp( -q0_liquid_invR * invT )
               diffA2(ic0,ic1,ic2) = d_solidA *
     &                               exp( -q0_solidA_invR * invT )

            end do
         end do
      end do
c
      return
      end

      subroutine concentration_pfmdiffusion_of_temperature_multiphases(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phi, nphi, ngphi,
     &   diffL0, diffL1, diffL2,
     &   diffA0, diffA1, diffA2, ngdiff,
     &   temp, ngtemp,
     &   d_liquid, q0_liquid,
     &   d_solid, q0_solid,
     &   gas_constant_R,
     &   interp_type,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphi, ngphi, ngdiff, ngtemp
      character*(*) avg_type, interp_type
      double precision d_liquid, d_solid
      double precision q0_liquid, q0_solid
      double precision gas_constant_R
c
c variables in 2d cell indexed
      double precision phi(CELL3d(ifirst,ilast,ngphi),nphi)
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffL0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffL1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffL2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ic2, ip
      double precision vphi, phil, phis, invT
      double precision q0_liquid_invR, q0_solid_invR
      double precision diff_liquid, diff_solid
      double precision interp_func
      double precision average_func

      q0_liquid_invR = q0_liquid / gas_constant_R
      q0_solid_invR  = q0_solid / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               vphi = 0.d0
c assuming the first nphi-1 order parameters are solid phase
               do ip = 1, nphi-1
                  vphi = vphi + average_func(
     &               phi(ic0-1,ic1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phis = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0-1,ic1,ic2,nphi), phi(ic0,ic1,ic2,nphi),
     &            avg_type )
               phil = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0-1,ic1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid  = d_solid * exp( -q0_solid_invR * invT )

               diffL0(ic0,ic1,ic2) = phil * diff_liquid
               diffA0(ic0,ic1,ic2) = phis * diff_solid
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0
               vphi = 0.d0
               do ip = 1, nphi-1
                  vphi = vphi + average_func(
     &               phi(ic0,ic1-1,ic2,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phis = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1-1,ic2,nphi), phi(ic0,ic1-1,ic2,nphi),
     &            avg_type )
               phil = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1-1,ic2)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid  = d_solid * exp( -q0_solid_invR * invT )

               diffL1(ic0,ic1,ic2) = phil * diff_liquid
               diffA1(ic0,ic1,ic2) = phis * diff_solid
            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0
               vphi = 0.d0
               do ip = 1, nphi-1
                  vphi = vphi + average_func(
     &               phi(ic0,ic1,ic2-1,ip), phi(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               phis = interp_func( vphi, interp_type )

               vphi = average_func(
     &            phi(ic0,ic1,ic2-1,nphi), phi(ic0,ic1-1,ic2,nphi),
     &            avg_type )
               phil = interp_func( vphi, interp_type )

               invT = 2.0d0 / ( temp(ic0,ic1,ic2-1)+temp(ic0,ic1,ic2) )

               diff_liquid = d_liquid * exp( -q0_liquid_invR * invT )
               diff_solid  = d_solid * exp( -q0_solid_invR * invT )

               diffL2(ic0,ic1,ic2) = phil * diff_liquid
               diffA2(ic0,ic1,ic2) = phis * diff_solid
            end do
         end do
      end do

      return
      end
c
c add interface diffusion to A and B diffusion
c
      subroutine ab_diffusion_of_temperature(
     &   ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2,
     &   phia, nphia, phib, nphib, ngphi,
     &   diffA0, diffA1, diffA2,
     &   diffB0, diffB1, diffB2, ngdiff,
     &   temp, ngtemp,
     &   d0, q0,
     &   gas_constant_R,
     &   avg_type)
c***********************************************************************
      implicit none
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2
      integer nphia, nphib
      integer ngphi, ngdiff, ngtemp
      character*(*) avg_type
      double precision d0, q0
      double precision gas_constant_R
c
c
c variables in 2d cell indexed
      double precision phia(CELL3d(ifirst,ilast,ngphi),nphia)
      double precision phib(CELL3d(ifirst,ilast,ngphi),nphib)
      double precision temp(CELL3d(ifirst,ilast,ngtemp))
      double precision diffA0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffA1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffA2(SIDE3d2(ifirst,ilast,ngdiff))
      double precision diffB0(SIDE3d0(ifirst,ilast,ngdiff))
      double precision diffB1(SIDE3d1(ifirst,ilast,ngdiff))
      double precision diffB2(SIDE3d2(ifirst,ilast,ngdiff))
c
      integer ic0, ic1, ic2, ip, jp
      double precision vphi, pa, pb, invT
      double precision q0_invR
      double precision dAB
      double precision average_func
c
      q0_invR = q0 / gas_constant_R
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0+1
               vphi = 0.d0
               do ip = 1, nphia
                  vphi = vphi + average_func(
     &               phia(ic0-1,ic1,ic2,ip), phia(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               pa = vphi

               do ip = 1, nphib
                  vphi = vphi + average_func(
     &               phib(ic0-1,ic1,ic2,ip), phib(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               pb = vphi

               invT = 2.0d0 / ( temp(ic0-1,ic1,ic2)+temp(ic0,ic1,ic2) )

               dAB = 16.d0*pa*pa*pb*pb*d0*exp(-q0_invR*invT)

               diffA0(ic0,ic1,ic2) = diffA0(ic0,ic1,ic2) + dAB
               diffB0(ic0,ic1,ic2) = diffB0(ic0,ic1,ic2) + dAB
            end do
         end do
      end do
c
c
      do ic2 = ifirst2, ilast2
         do ic1 = ifirst1, ilast1+1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphia
                  vphi = vphi + average_func(
     &               phia(ic0,ic1-1,ic2,ip), phia(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               pa = vphi

               vphi = 0.d0
               do ip = 1, nphib
                  vphi = vphi + average_func(
     &                phib(ic0,ic1-1,ic2,ip), phib(ic0,ic1,ic2,ip),
     &                avg_type )
               enddo
               pb = vphi

               invT = 2.0d0 / ( temp(ic0,ic1-1,ic2)+temp(ic0,ic1,ic2) )

               dAB = 16.d0*pa*pa*pb*pb*d0*exp(-q0_invR*invT)

               diffA1(ic0,ic1,ic2) = diffA1(ic0,ic1,ic2) + dAB
               diffB1(ic0,ic1,ic2) = diffB1(ic0,ic1,ic2) + dAB

            end do
         end do
      end do
c
      do ic2 = ifirst2, ilast2+1
         do ic1 = ifirst1, ilast1
            do ic0 = ifirst0, ilast0

               vphi = 0.d0
               do ip = 1, nphia
                  vphi = vphi + average_func(
     &               phia(ic0,ic1,ic2-1,ip), phia(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               pa = vphi

               vphi = 0.d0
               do ip = 1, nphib
                  vphi = average_func(
     &               phib(ic0,ic1,ic2-1,ip), phib(ic0,ic1,ic2,ip),
     &               avg_type )
               enddo
               pb = vphi

               invT = 2.0d0 / ( temp(ic0,ic1,ic2-1)+temp(ic0,ic1,ic2) )

               dAB = 16.d0*pa*pa*pb*pb*d0*exp(-q0_invR*invT)

               diffA2(ic0,ic1,ic2) = diffA2(ic0,ic1,ic2) + dAB
               diffB2(ic0,ic1,ic2) = diffB2(ic0,ic1,ic2) + dAB

            end do
         end do
      end do

      return
      end
