c
c ----------------------------------------------------------------
c
      program driver
      implicit double precision (a-h,o-z)
      parameter(maxmx = 500)
      parameter(maxmy = 1000)
      parameter(mwork = 9736792)
      parameter(mbc = 2)
      parameter(meqn = 6)
      parameter(mwaves = 3)
      parameter(maux = 7)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      dimension mthlim(mwaves)
      dimension work(mwork)
      real*4 cpusec,cputm(2)
c
      call cpu_time(cpustart)
c
      call claw2ez(maxmx,maxmy,meqn,mwaves,mbc,maux,mwork,
     &     mthlim,q,work,aux)
c
      call cpu_time(cpuend)
      cpusec = cpuend-cpustart
c
      write(23,*) cpusec
c
      write(10,*)
     &    'Program has been running for', cpusec, 'seconds.'
      stop
c
      end
