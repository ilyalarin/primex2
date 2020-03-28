/*
      gfortran -mcmodel=medium yinit_c12.F -o yinit_c12

      INPUT:
      tables/eche.dat, 
      tables/pi0_coul_c12.dat, tables/pi0_stro_c12.dat,
      tables/inc_c12.dat, 
      tables/tmt_c12.dat

      OUTPUT:
      tables/dfp_c12.dat
*/
      program yinit_c12
      implicit none
*
*     parameters
*
      double precision PI, TWOPI, DEGRAD
      PARAMETER (PI     = 3.141592653589793238d0)
      PARAMETER (TWOPI  = 2.0d0*PI)
      PARAMETER (DEGRAD = PI/180.d0)
*
      integer nech, maxtht, maxthr, nterms
      parameter(nech = 180, maxtht = 5500, maxthr = 550, nterms = 5)
*
*     common variables and blocks
*
      integer icmodel
      double precision    ksi, alphae
      common/cparam_com/  ksi, alphae, icmodel
*
      real    fc_a1, fc_a2, fs_b1, fs_b2, fs_c1, fs_c2
      common/ff_com/fc_a1(nech,maxtht), fc_a2(nech,maxtht),
     &              fs_b1(nech,maxtht), fs_b2(nech,maxtht),
     &              fs_c1(nech,maxtht), fs_c2(nech,maxtht)
*
*     local variables
*
      integer iargc, i, j, k, i1, ic, jc, ii,
     &        ipwo, ibinning, nbinr, ntetdv, nsub,
     &        imin(maxtht,nech), imax(maxtht,nech)
*
      real    ffd1, ffd2, ffd3, ffd4, fk1, fk2,
     &        etc(nech), rmat(maxtht,maxthr,nech),
     &        dfprob(nterms,maxthr,nech)
*
      double precision dth, th, dom, dsdo(nterms), dsdoic(nech,maxtht),
     &        ebeam
      character*(4) cha
*
*     preset default values
*
      ksi     = 0.25d0
      alphae  = 1.2d0 ! 1.6d0
      icmodel = 1
*
*     read arguments
*
      if(iargc().lt.2.or.iargc().gt.5) then
        print*, 'should be 2 - 5 arguments'
        call exit(1)
      elseif(iargc().eq.3) then
        call getarg(3,cha)
        if(cha(1:1).eq.'c') then        ! Cornell
          icmodel = 1
        elseif(cha(1:1).eq.'s') then    ! Sibirtsev
          icmodel = 2
        elseif(cha(1:1).eq.'l') then    ! Laget
          icmodel = 3
        else
          print*, '3rd argument should be c, s or l'
          call exit(1)
        endif
      elseif(iargc().eq.4) then
        call getarg(4,cha)
        read(cha,*) ksi
      elseif(iargc().eq.5) then
        call getarg(5,cha)
        read(cha,*) alphae
      endif
*
      call getarg(1,cha)
      read(cha(1:1),'(i1)') ibinning
      if(ibinning.lt.1.or.ibinning.gt.3) then
        print*, 'Bad bining option ', cha(1:1)
        call exit(1)
      endif
      call getarg(2,cha)
      read(cha(1:2),'(i2)') ipwo
      if(ipwo.lt.-1.or.ipwo.gt.1) then
        print*, 'Bad pwo option ', cha
        call exit(1)
      endif
      write(*,'(a,f4.2,a,i2)')
     &        'binning = ', ibinning*0.01, ' PWO option = ', ipwo
      print*, 'coh. model = ', icmodel,
     &        ' ksi = ', sngl(ksi), ' alphae = ', sngl(alphae)
*
*     set up parameters
*
      if(ipwo.eq.1) then
        ntetdv = 3300
c       thmax  = 3.3d0*DEGRAD
        if(ibinning.eq.1) then
          nbinr = 330
          nsub  = 10
        elseif(ibinning.eq.2) then
          nbinr = 165
          nsub  = 20
        elseif(ibinning.eq.3) then
          nbinr = 110
          nsub  = 30
        endif
      else
        ntetdv = 5460
c       thmax  = 5.46d0*DEGRAD
        if(ibinning.eq.1) then
          nbinr = 546
          nsub  = 10
        elseif(ibinning.eq.2) then
          nbinr = 273
          nsub  = 20
        elseif(ibinning.eq.3) then
          nbinr = 182
          nsub  = 30
        endif
      endif
*
      dth = 1.d-3*DEGRAD
*
*     read input tables
*
c      DO IA = 28, 30
c        write(cha(1:2),'(i2)') IA
*
        open(20, file='tables/pi0_coul_c12.dat',
     &           form='formatted', status='old')
        do i = 1, nech
          do j = 1, maxtht
            read(20,'(i3,1x,i4,1x,2(1x,E11.5))') ic, jc, ffd1, ffd2
            if(i.ne.ic.or.j.ne.jc) then
              print*, 'Bad data given for coul ff ', i, ic, j, jc
              stop
            endif
*
            fc_a1(i,j) = ffd1
            fc_a2(i,j) = ffd2
*
          enddo
        enddo
        close(20)
*
        open(20, file='tables/pi0_stro_c12.dat',
     &           form='formatted', status='old')
        do i = 1, nech
          do j = 1, maxtht
            read(20,'(i3,1x,i4,1x,4(1x,E11.5))') ic, jc, ffd1, ffd2,
     &                                                   ffd3, ffd4
            if(i.ne.ic.or.j.ne.jc) then
              print*, 'Bad data given for stro ff ', i, ic, j, jc
              stop
            endif
*
            fs_b1(i,j) = ffd1
            fs_b2(i,j) = ffd2
            fs_c1(i,j) = ffd3
            fs_c2(i,j) = ffd4
*
          enddo
        enddo
        close(20)
c      ENDDO
*
      OPEN(20, file='tables/eche.dat', form='formatted', status='old')
      do i = 1, nech
        read(20,'(f7.5)') etc(i)
      enddo
      close(20)
*
      open(20, file='tables/tmt_c12.dat',form='formatted', status='old')
*
      do j = 1, nech
        do i = 1, ntetdv
          read(20,'(i3,1x,i4,2(1x,i3))') jc, ic, imin(i,j), imax(i,j)
*
          if(i.ne.ic.or.j.ne.jc) then
            print*, 'Bad data given for tmt'
            stop
          endif
*
          do i1 = 1, nbinr
            rmat(i,i1,j) = 0.
          enddo
          do i1 = imin(i,j), imax(i,j)
            read(20,'(1x,f9.7)') rmat(i,i1,j)
          enddo
*
        enddo
      enddo
      close(20)
*
      open(20, file='tables/inc_c12.dat',form='formatted', status='old')
      do j = 1, nech
        do i = 1, maxtht
          read(20,'(i3,1x,i4,1x,E12.6)') jc, ic, dsdoic(j,i)
          if(i.ne.ic.or.j.ne.jc) then
            print*, 'Bad data given for inc_cs ', i, ic, j, jc
            stop
          endif
          dsdoic(j,i) = dsdoic(j,i)*1.d-3         ! from ubarn/sr to mbarn/sr
        enddo
      enddo
      close(20)
*
*     get folded with res. unitary c.s. array
*
      do j = 1, nech
        do i1 = 1, nbinr
          do k = 1, nterms
            dfprob(k,i1,j) = 0.
          enddo
        enddo
      enddo
*
      do j = 1, nech
        ebeam = etc(j)
        do i = 1, ntetdv
*
          th  = (DBLE(i)-0.5d0)*dth
          dom = TWOPI * (DCOS(th-dth*0.5d0) - DCOS(th+dth*0.5d0))
          CALL CRSEC6(ebeam,th,dsdo,j,i)
          dsdo(4) = dsdoic(j,i)
*
          do i1 = imin(i,j), imax(i,j)
            do k = 1, nterms
              dfprob(k,i1,j) = dfprob(k,i1,j)+dsdo(k)*dom*rmat(i,i1,j)
            enddo
          enddo
        enddo
      enddo
*
*     dump folded with res. unitary c.s. array into output file
*
      open(20, file='tables/dfp_c12.dat',form='formatted', status='new')
*
      do j = 1, nech
        do i1 = 1, nbinr
          do k = 1, nterms
            write(20,'(i3,1x,i3,1x,i1,1x,E12.5)')
     &        j, i1, k, dfprob(k,i1,j)
          enddo
        enddo
      enddo
      close(20)
      end
*
*     Calculate parametrized cross-section:
*
      SUBROUTINE CRSEC6(ebeam,thpi0,sumall,ie,it)
      implicit none
*
*     parameters
*
      double precision PI, TWOPI, DEGRAD
      PARAMETER (PI     = 3.141592653589793238d0)
      PARAMETER (TWOPI  = 2.0d0*PI)
      PARAMETER (DEGRAD = PI/180.d0)
*
      integer nech, maxtht, maxthr, nterms
      parameter(nech = 180, maxtht = 5500, maxthr = 550, nterms = 5)
*
*     arguments
*
      integer ie, it
      double precision ebeam, sumall(nterms), thpi0
*
*
*     common variables and blocks
*
      integer icmodel
      double precision    ksi, alphae
      common/cparam_com/  ksi, alphae, icmodel
*
      real    fc_a1, fc_a2, fs_b1, fs_b2, fs_c1, fs_c2
      common/ff_com/fc_a1(nech,maxtht), fc_a2(nech,maxtht),
     &              fs_b1(nech,maxtht), fs_b2(nech,maxtht),
     &              fs_c1(nech,maxtht), fs_c2(nech,maxtht)
*
*     specific parameters
*
      double precision PI0MAS, U, ALPHA, gev2mb, Gpi02g, Coul_coef, mN
      PARAMETER (PI0MAS = 0.1349766D0, U = 0.931494D0, mN = 0.93892d0)
      PARAMETER (ALPHA  = 1.D0/137.0359895D0)
      PARAMETER (gev2mb = 0.38937966D0)           ! GeV^-2 to millibarn
      PARAMETER (Gpi02g = 1.0D-9)                 ! eV to GeV
      PARAMETER (Coul_coef = 8.D0*gev2mb*ALPHA*Gpi02g)
*
*     local variables
*
      double precision dsdopr, dsdonc, dsdint, dsdin2, dcosth, dsinth,
     &  C1, A1, A2, A3, RECMAS, TKINRC, Delta2,
     &  PPI0DP, EPI0DP, PI0BET, PI0BEM, DISCRI, ARECOI, TCOULM, TNCOHE,
     &  qMeV, qfm, Z, FC1, FC2, FN1, FN2
*
      double precision a_c12
      data a_c12 /12.d0/
      save a_c12
*
      sumall(1) = 0.d0
      sumall(2) = 0.d0
      sumall(3) = 0.d0
      sumall(5) = 0.d0
*
c      DO IA = 28, 30
        Z      = 6.0D0                 ! C12
        ARECOI = a_c12
        RECMAS = ARECOI*U
        dcosth = DCOS(thpi0)
        dsinth = DSIN(thpi0)
*
        C1     = (2.D0*ebeam*RECMAS+PI0MAS*PI0MAS)*0.5d0
        A1     = ebeam*ebeam*dcosth**2-(ebeam+RECMAS)**2
        A2     = 2.D0*ebeam*C1*dcosth
        A3     = C1*C1 - PI0MAS*PI0MAS*(ebeam+RECMAS)**2
C
        DISCRI = A2*A2-4.D0*A1*A3
        if(DISCRI.Lt.0.D0) then
          print*, 'CRSEC6 ERR: Discriminant is negative'
          stop
        endif
        PPI0DP = (-A2-DSQRT(DISCRI))/(2.D0*A1)
C
        EPI0DP = DSQRT(PPI0DP*PPI0DP+PI0MAS*PI0MAS)
        TKINRC = ebeam-EPI0DP
C
cc      t      = PI0MAS*PI0MAS-2.D0*(ebeam*EPI0DP-ebeam*PPI0DP*dcosth)
cc      t1     = (ebeam*dsinth)**2 ! -t ! 4.d0 * ebeam * PPI0DP * dsin(thpi0*0.5d0)**2
C
        Delta2 = -2.D0*RECMAS*TKINRC      ! the same as t in Mandelstamm
C
        qMeV   = 1.D3*DSQRT(2.D0*RECMAS*TKINRC)
        qfm    = qMeV/197.327053D0        ! in 1/fm
C
        FC1    = fc_a1(ie,it) / (1.d0+1.729d0*abs(Delta2))
        FC2    = fc_a2(ie,it) / (1.d0+1.729d0*abs(Delta2))
        FN1    = fs_b1(ie,it) - ksi * fs_c1(ie,it)
        FN2    = fs_b2(ie,it) - ksi * fs_c2(ie,it)
C
        PI0BET = PPI0DP/EPI0DP
        PI0BEM = PI0BET/PI0MAS
C
        TCOULM = Z * DSQRT(Coul_coef*PI0BEM) * PI0BEM
     &             * ebeam**2 * dsinth / ABS(Delta2)
C
        if(icmodel.eq.1) then           ! Cornell (pure Cornell if alphae = 1):
*
          TNCOHE = 10.0d0*ebeam**(alphae-1.d0)
*
        elseif(icmodel.eq.2) then       ! Sibirtsev
*
          TNCOHE = 2.28d0*(mN*(2.d0*ebeam+mN))**0.641d0
*
        elseif(icmodel.eq.3) then       ! Laget
*
          TNCOHE = 5.30d0*(mN*(2.d0*ebeam+mN))**0.440d0
*
        else                            ! unknown value
*
          TNCOHE = 0.d0
*
        endif
*
        TNCOHE = TNCOHE * ARECOI * ebeam * DSQRT(1.D-3) * dsinth
*
        dsdopr = TCOULM*TCOULM*(FC1*FC1+FC2*FC2)    ! in milibarn/sr
        dsdonc = TNCOHE*TNCOHE*(FN1*FN1+FN2*FN2)    ! in milibarn/sr
        dsdint = 2.D0*TCOULM*TNCOHE*(FC1*FN1+FC2*FN2)
        dsdin2 = 2.D0*TCOULM*TNCOHE*(FC2*FN1-FC1*FN2)
C
        sumall(1) = sumall(1) + dsdopr
        sumall(2) = sumall(2) + dsdint
        sumall(3) = sumall(3) + dsdonc
        sumall(5) = sumall(5) + dsdin2
*
c      ENDDO
*
      return
      end
*
