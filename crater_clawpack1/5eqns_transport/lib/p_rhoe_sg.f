c
c ------------------------------------------------------------
c
      double precision function p_rhoe_sg(rhoe,geos0,bn0)
      implicit double precision (a-h,o-z)
c
c     # compute pressure as function of 
c     # density and internal energy
c     # stiffened gas EOS
c
      p_rhoe_sg = (geos0-1.d0)*rhoe-geos0*bn0
      return
c
      end
