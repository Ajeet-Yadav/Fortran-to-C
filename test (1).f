c alfa.1, beta.1 KAW and .05 of whistler,a0 1 and .08


      program dzhkVT

      implicit none
      
      integer*4 nt, m, n
      integer*4 ndg, ndg1, ndg2
      integer*4 ntd
      parameter (nt=500, m=11, n=2**m)
      parameter (ndg=2, ndg1=1, ndg2=1)
      parameter (ntd=nt*ndg2/ndg1)
       
      complex*16 e(n), b(n),b2(n), a(n),a2(n),efck(n)
      complex*16 d(n), dnck(n)
      complex*16 v(n)
      
      complex*16 ae(n), be(n)
      complex*16 ab(n), bb(n)
      complex*16 ad(n), bd(n)
      complex*16 ab2(n), bb2(n)
      complex*16 ew(n)
      complex*16 av(n), bv(n)
      complex*16 fe(n),fd(n)
      complex*16 ak(n), akC(n),ak2(n),akC2(n)
      
      complex*16 feC(n),fdC(n)
      complex*16 dc(n),vc(n)
      complex*16 xte(n),xtd(n)
      
      complex*16 te(n),td(n)
      complex*16 tb(n),tb1(n),tb2(n),tw(n)
      complex*16 xe(n),xd(n),xb(n),xb1(n),xb2(n)       
      complex*16 ykw(n),ykw2(n)
                      
      real*8 k(n),k2(n),kk(n)
      real*8 ef2(n), ef2k(n),den(n),denk(n),ef(n),denck(n)
      real*8 den2k(n)
      real*8 c1,c2,c3,c4,c5,c6,g

      real*8     a0
      real*8     beta, gam, alf
      real*8     vel
      real*8     dt
      integer*4  niter, it
      real*8     ttime, TMAX,DTMX,DTMN
      real*8     nrm, LONG, dx 
      real*8     pi,pi2,th
      complex*16 ui 
      real*8     plsn,wn
      logical precor
      integer*4  i
      integer*4  i20
      integer*4  idt
      integer*4  iali
      
      character*4 f0,f1,f2,f3,f4,f5,f6
    
      character*5  f7,f8

      common /kw/ k /k2/ k2 /ykw/ ykw /ykw2/ ykw2
      common /te/ te /td/ td /tb1/ tb1 /tb2/ tb2 
      common /ew/ ew /tw/ tw
      common /xe/ xe /xd/ xd /xb1/ xb1 /xb2/ xb2
      common /norm/ nrm
      common /step/ dx
      common /dato/ a0, beta, gam, alf,vel
      data idt    /0/
      pi=4.d0*datan(1.d0)
      pi2     = 2.d0*pi
    
      ui      = (0.d0,1.d0)
      dt    = .00005d0
      TMAX   = dfloat(nt) 
      data precor  /.TRUE./
 
      ttime = 0.d0
       
      niter = int(float(ndg2)/float(ndg1)+.5d0)  
      it    = 0
        
      DTMX  = 1.d-5
      DTMN   =DTMX/50.d0
       alf   = 0.1d0
      LONG = pi2/alf
      beta  = 0.1d0
       gam = 0.d0  
      a0 =1.0d0 
      vel = 0.d0 
      nrm     = 1.d0/dfloat(n)
      dx      = LONG/dfloat(n)

        c1=0.0011d0
        c2=9.9d0
        c3=160.26d0
        
        c5=0.0d0
        
  
       f0 = 'ef00'
       call filen(f0, 4)
       open(20, file = f0)

       f1 = 'dn00'
       call filen(f1, 4)
        open(21, file = f1)

      
       f2 = 'ek00'
       call filen(f2, 4)
      open(22, file = f2)

       f3 = 'dk00'
       call filen(f3, 4)
       open(23, file = f3)         
    
        f4 = 'dz00'
        call filen (f4 ,4)
        open (24,file = f4)

        f5 = 'eu00'
        call filen(f5, 4)
        open(25, file = f5)  
      
c      f6 = 'ev00'
c      call filen(f6, 4)
c      open(26, file = f6)

c      f8 = 'fpu00'
c      call filen(f8, 4)
c      open(30, file = f8)

       f7 = 'prm00'
      call filen(f7, 5)
      open(15, file = f7)


      write(15,9000) nt
9000  format(' nt   = ', i6)

      write(15,9001) n
9001  format(' n    = ', i6)

      write(15,9002) ndg, ndg1, ndg2
9002  format(' ndg  = ',i6,
     .       ' ndg1 = ',i6,
     .       ' ndg2 = ',i6)
      
      write(15,9003) LONG
9003  format(/,' LONG = ', 1pe10.3)

      write(15,9004) dt, dx
9004  format(' dt   = ',1pe10.3,
     .       ' dx   = ',1pe10.3)

      write(15,9005) a0
9005  format(/,' A0   = ',1pe10.3)

      write(15,9006) beta, gam, alf
9006  format(' beta = ',1pe10.3,
     .       ' gam  = ',1pe10.3,
     .       ' alf  = ',1pe10.3)
       
          

       do i = 1, n/2
         k(i) = -dfloat(i-1)*pi2/LONG
         k(n-i+1) = dfloat(i)*pi2/LONG
         kk(n/2-i+1)=k(i)         
      enddo
       do i=n/2+1,n
         kk(n-i+1+n/2)=k(i)
        end do 
       do i=1,n
          k2(i)=k(i)**2
          ykw(i)=ui*k2(i)
          ykw2(i)=(-ui*c1*k2(i)**2)+ui*c2*k2(i)

       enddo
      
      call coef (ak,akC,n,dt)
      call coef2(ak2,akC2,n,dt)
      
  
      call intdat(e,d,n,m)
      
      write(15,*) ' Number of Plasmons'
      call intcst(e,n,plsn,15)
      write (*,*) 'plasmon1'
      write (*,*) plsn
       
      
      do i=1,n
        xe(i)  = e(i)
        xd(i) = d(i)
      end do         
            
       call dcfftn(-1,xe,xte,m)
       call dcfftn(-1,xd,xtd,m)
      do i = 1, n
         ef2(i) = (nrm*abs(xte(i)))**2
         den(i) = (nrm*abs(xtd(i)))**2
  	   ef(i)=nrm*xte(i)
c?	 den(i)=nrm*xtd(i)
      enddo
        
      call nlterm(e,d,fe,fd,n,m,dt,c3)   
      do i = 1, n         
               b(i) = e(i)
	       e(i)=b(i)-ui*k2(i)*dt*e(i)+.5d0*fe(i)
	       b2(i)=d(i)
               d(i)=b2(i)+ui*c1*dt*(k2(i)**2)*d(i)
     .              -ui*c2*k2(i)*dt*d(i)+.5d0*fd(i)
       enddo
       call intcst (e,n,plsn,15)
       write (*,*) plsn
            
      call transl4(e,ae,b,ab,d,ad,b2,ab2,n,-1)          
      
      
100   continue

      if (ttime .le. TMAX) then
         

         do i20 = 1, 20
         do iali = 1, 2
         
            if (iali .eq. 2) then
               do i = 1, n
                  e(i)  = ae(i)
                  b(i)  = ab(i)
                  d(i)  = ad(i)
                  b2(i) = ab2(i)
               enddo
            endif
       
         call nlterm(e,d,fe,fd,n,m,dt,c3)   
            do i = 1, n
               a(i)=e(i)
               e(i)=b(i)+ ak(i)*e(i)+ fe(i)
               b(i)=a(i) 
               
               a2(i)=d(i)
               d(i)=b2(i)+ ak2(i)*d(i)+ fd(i)
               b2(i)=a2(i)
           enddo
	 
          if (precor) then
        call nlterm(e,d,feC,fdC,n,m,dt,c3)   
	  do i=1,n
           e(i)=b(i)+.50d0*akC(i)*(b(i)+e(i))+.25d0*(fe(i)+feC(i))
           d(i)=b2(i)+.50d0*akC2(i)*(b2(i)+d(i))+.25d0*(fd(i)+fdC(i))
	  enddo
	  endif


            if (iali .eq. 1) then
               do i = 1, n
                  be(i) = e(i)                 
                  bb(i) = b(i)
                  bd(i) = d(i)
                  bb2(i)=b2(i)
               enddo
            else
              do i = 1, n
                  ae(i) = e(i)
                  ab(i) = b(i)
                  ad(i) = d(i)
                  ab2(i)=b2(i)
               enddo
            endif

c  end iali: 1, 2
         enddo
         
            call antial4(e,be,ae,b,bb,ab,d,bd,ad,b2,bb2,ab2,n)  
            call transl4(e,ae,b,ab,d,ad,b2,ab2,n,-1)          

c  end i20: 1, 20
                   
         enddo
         
            call intcst(e,n,plsn,15)
         
         ttime = ttime + 20.d0*dt

           

         if (int(niter*ttime) .eq. it) then
            it = it + 1
             write(*,*)  ttime
c            write(26,*)  ttime
         
                     
            do i = 1, n/2
              ef2k(n/2-i+1)  = nrm*abs(e(i))
              den2k(n/2-i+1)  = nrm*abs(d(i))
               denck(n/2-i+1) = nrm*(d(i))
           	efck(n/2-i+1) = nrm*(e(i))
c                write(*,*) i,k(i)
            enddo
            do i = n/2+1, n
               ef2k(n-i+1+n/2)  = nrm*abs(e(i))
               den2k(n-i+1+n/2)  = nrm*abs(d(i))
               denck(n-i+1+n/2) = nrm*(d(i))
      	       efck(n-i+1+n/2) = nrm*(e(i))
c               write(*,*) i,k(i)
            enddo


            do i =1,n
                xe(i)  = e(i)
                xd(i)  = d(i)
            enddo
             call dcfftn(-1,xe,xte,m)
             call dcfftn(-1,xd,xtd,m)
            
            do i = 1, n
               ef2(i) = (nrm*abs(xte(i)))**2
               den(i) = (nrm*abs(xtd(i)))**2
       	       ef(i) = nrm*abs(xte(i))	      
c?	       den(i)=nrm*abs(xtd(i)
            enddo

            write(*,*)  it, plsn
            write(15,*) it, plsn
c            write(*,*) it, den

            call wrtb1(ef2,n,1,20)
            call wrtb1(den,n,1,21)
            call wrtb1(ef2k,n,1,22)
            call wrtb1(den2k,n,1,23)   
            write(24,*) ef2k(n/2)    !used for stem plot
            write(25,*) den2k(n/2)   !used for stem plot
           
         endif
         goto 100
      endif

      close(20)
      close(21)
      close(22)
       close(23)
      close(24)
      close(25)
        close(26)
      close(30)   
      close(15)
      
  
    
      
      write(*,*) '*** end ***'
      stop
      end

      subroutine wrtb1(tb,n,ndiag,uni)
      implicit none
      
      integer uni, n, ndiag, i
      real*8 tb(n)

      do i = 1, n, ndiag
         write(uni, 9000) tb(i)
      enddo

      return
9000  format(1pe15.5)
      end

      subroutine intcst(tb1,n,plas,uni)
      implicit none
      
      integer uni, n
      complex*16 tb1(n)
      
      real*8 plas
      real*8 nrm, nrm2
      real*8 t1
      integer*4 i
      
      common /norm/ nrm

      nrm2 = nrm**2
      plas = 0.d0
      do i = 1, n
         t1 = abs(tb1(i))**2
         plas = plas + t1
      enddo

      plas = nrm2*plas
      
      return
      end

      subroutine coef(ak,akC,n,dt)
      implicit none
      
      integer*4 n
      complex*16 ykw(1),ykw2(1)
      complex*16 ak(n),akC(n)
      
      real*8 dt
     
      complex*16 ui
      integer*4 i
      
      real*8 k(1), k2(1)
      
       common /kw/ k/k2/k2/ykw/ykw    
 
      ui = (0.d0,1.d0)
      do i = 1, n
         ak(i)  = (cdexp(-ykw(i)*dt)-cdexp(ykw(i)*dt))
         akC(i) = (cdexp(-.5d0*ykw(i)*dt)-cdexp(.5d0*ykw(i)*dt))
      enddo
      
      return
      end

      subroutine coef2(ak2,akC2,n,dt)
      implicit none
      
      integer*4 n
      complex*16 ykw2(1),ykw(1)
      complex*16 ak2(n),akC2(n)
      
      real*8 dt
     
      complex*16 ui
      integer*4 i
      
      real*8 k(1), k2(1)
      
      common /kw/ k /k2/ k2 /ykw/ ykw /ykw2/ ykw2     
 
      ui = (0.d0,1.d0)
      do i = 1, n
         ak2(i)  = (cdexp(-ykw2(i)*dt)-cdexp(ykw2(i)*dt))
         akC2(i) = (cdexp(-.5d0*ykw2(i)*dt)-cdexp(.5d0*ykw2(i)*dt))
      enddo
      
      return
      end

      subroutine nlterm(e,d,fe,fd,n,m,dt,c3)
      implicit none

      integer*4 n, m
      complex*16 e(n),b(n),b2(n)
      complex*16 d(n)
      complex*16 fe(n)
      complex*16 fd(n)
      complex*16 ew(n),tw(n)

      complex*16 te(1),td(1),tb1(1),tb2(1)
      complex*16 xe(1),xd(1),xb1(1),xb2(1)
      real*8 k(1),k2(1),g
      real*8 alf,beta, dt,x
      real*8 nrm, nrm2
      real*8   c3,dx
      complex*16 ui
      integer*4 i
      
      common /te/ te /td/ td /tb1/ tb1 /tb2/ tb2
      common /xe/ xe /xd/ xd /xb1/ xb1 /xb2/ xb2
      common /kw/ k  /k2/ k2
      common /norm/ nrm
      common /step/ dx

 
      ui   = (0.d0,1.d0)
      nrm2 = nrm**2
  
       do i = 1, n
         xe(i) = e(i)
         xb1(i)= b(i)
         xd(i) = d(i)
         xb2(i)= b2(i)
      enddo

       call dcfftn(-1,xe,te,m)
       
       call dcfftn(-1,xd,td,m)
      
      do i = 1,n
	     xe(i) = nrm*te(i)*(dcmplx((abs(nrm*
     .             te(i)))**2,0.d0))
 	     xd(i) = nrm*td(i)*(dcmplx((abs(nrm*
     .             te(i)))**2,0.d0))

c                xd(i)=(dexp((c5*abs(nrm*te(i)))**2+
c     .       (abs(nrm*td(i)))**2)-1.d0)*nrm*td(i)
c                xe(i)=(dexp((abs(nrm*te(i)))**2+
c     .       (c5*abs(nrm*td(i)))**2)-1.d0)*nrm*te(i)

      enddo
    
        call dcfftn(1,xe,te,m)
      
        call dcfftn(1,xd,td,m)
       

      do i = 1, n
         fd(i)=2.d0*ui*c3*dt*td(i)
         fe(i)=2.d0*ui*dt*te(i)
      enddo
      
      return
      end

      
      subroutine antial4(h1,f1,g1,h2,f2,g2,h3,f3,g3,h4,f4,g4,n)
      implicit none
      
      integer*4 n
      complex*16 h1(n), f1(n), g1(n)
      complex*16 h2(n), f2(n), g2(n) 
      complex*16 h3(n), f3(n), g3(n)
      complex*16 h4(n), f4(n), g4(n)
      real*8 k(1)
      real*8 dx
      real*8 dx2
      real*8 a
      complex*16 c
      complex*16 ui
      integer*4 i
      
      common /kw/ k
      common /step/ dx

      ui = (0.d0,1.d0)
      
      dx2 = .5d0*dx
          
      do i = 1, n
        a = dx2*k(i)
         c = dcmplx( dcos(a), -dsin(a) )
         h1(i) = .5d0*( f1(i) + c*g1(i) )
         h2(i) = .5d0*( f2(i) + c*g2(i) )
         h3(i) = .5d0*( f3(i) + c*g3(i) )
         h4(i) = .5d0*( f4(i) + c*g4(i) )
      enddo

      return
      end

      subroutine transl4(h1,f1,h2,f2,h3,f3,h4,f4,n,isgn)
      implicit none
      
      integer*4 n, isgn
      complex*16 h1(n), f1(n)
      complex*16 h2(n), f2(n) 
      complex*16 h3(n), f3(n)
      complex*16 h4(n), f4(n)
      real*8 k(1),k2(1)
      real*8 dx
      real*8 dx2
      real*8 sg, a
      complex*16 ui
      complex*16 c
      integer*4 i
 
      common /kw/ k /k2/ k2
      common /step/ dx

      ui = (0.d0,1.d0)
      dx2 = .5d0*dx
      sg  = dfloat(isgn)
      
      do i = 1, n
          a = sg*dx2*k(i)
        
         c = dcmplx(dcos(a), -dsin(a))
         f1(i) = c*h1(i)
         f2(i)  =c*h2(i)
         f3(i) = c*h3(i)
         f4(i) = c*h4(i)
      enddo

      return
      end         

c  File names
      subroutine filen(fname, name)
      implicit none
      integer name
      integer i
      logical fexist
      character fname*(*)
      save

      do i = 1, 100
        inquire(file = fname, exist = fexist)
        if (.not. fexist) goto 16
        call bumpchar(fname, name)
      enddo
      stop
16    continue

      return
      end

      subroutine bumpchar(string, n)
      implicit none
      integer m, n
      character string*(*)
      m = n
      do m = n,1,-1
         if(string(m:m).ne.'9') then
            string(m:m) = char( ichar(string(m:m)) + 1 )
            return
         endif
         string(m:m) = '0'
      enddo
      return
      end

      SUBROUTINE dcfftn(ISIGN,DATA,FFT,NN)
      real*8 WR,WI,WPR,WPI,WTEMP,THETA
      real*8 TEMPR, TEMPI
      real*8 DATA(*)
      complex*16 FFT(*)
      N=2*(2**NN)
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/DFLOAT(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=WR*DATA(J)-WI*DATA(J+1)
            TEMPI=WR*DATA(J+1)+WI*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      DO I=1,N/2
        FFT(I) = DCMPLX(DATA(2*I-1),DATA(2*I))
      ENDDO
      RETURN
      END

      subroutine intdat(e,d,n,m)
      implicit none
      real *8 RAND, th
      integer*4 n,m
      complex*16 e(n)



      complex*16 d(n)
      complex*16 v(n)
      
           
      complex*16 te(1),td(1)
      complex*16 tv(1)
      real*8 pi2
      complex*16 ui
      
      real*8 nrm,ran
      real*8 dx,a1
      real*8 a0,b
      real*8 beta, gam, alf
      real*8 vel
      real*8 x,r
      real*8 r0
      integer*4 i,sd, r10
      common /te/ te /td/ td
      common /tv/ tv
      
      common /norm/ nrm
      common /step/ dx
      common /dato/ a0, beta, gam, alf, vel

      pi2  = 8.d0*datan(1.d0)
      ui   = (0.d0,1.d0)


    

         do i = 1, n, j
            x = dfloat(i-n/2)*dx
             te(i) = dcmplx(
     .             a0*(1.d0+beta*dcos(alf*x)),
     .             0.d0)
             td(i) = dcmplx(
     .             0.08d0*(1.d0+0.05d0*dcos(alf*x)),
     .             0.d0)



c             td(i)=0.01d0*dexp(-x**2/r0**2)                  
                                         
           enddo        

      call dcfftn(1,te,e,m)
      call dcfftn(1,td,d,m)

      
      return
      end



