c
c ------------------------------------------------------------
c
      double precision function p_rhoe_ideal(rhoe,geos0)
      implicit double precision (a-h,o-z)
c
c     # compute phasic pressure with both density and
c     # specific internal energy known apriori
c     # ideal gas EOS
c
      p_rhoe_ideal = (geos0-1.d0)*rhoe
      return
c
      end
