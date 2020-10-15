

      program dzhkVT
 
   
      implicit none
      
      integer*4 nt, mx,my,nx,ny,jj
      integer*4 ndg, ndg1, ndg2
      integer*4 ntd
      parameter (nt =20, mx =6, nx = 2**mx)
      parameter ( my =6, ny = 2**my)

      parameter (ndg =2 , ndg1 = 1, ndg2 =2 )
      parameter (ntd = nt*ndg2/ndg1)
       
      complex*16 e(nx,ny), b(nx,ny), a(nx,ny),efck(nx,ny)
      complex*16 d(nx,ny),dnck(nx,ny)
      complex*16 v(nx,ny)
      
      complex*16 ae(nx,ny), be(nx,ny)
      complex*16 ab(nx,ny), bb(nx,ny)
      complex*16 ad(nx,ny),bd(nx,ny)
      complex*16 av(nx,ny), bv(nx,ny)
      complex*16 fe(nx,ny), fv(nx,ny)
      complex*16 ak(nx,ny), akC(nx,ny)
      
      complex*16 feC(nx,ny),fvC(nx,ny)
      complex*16 dc(nx,ny),vc(nx,ny)
      complex*16 xte(nx,ny), xtd(nx,ny), denz(nx,ny)
      
      complex*16 te(nx,ny), tb(nx,ny),td(nx,ny),tv(nx,ny)

      complex*16 xe(nx,ny), xb(nx,ny),xd(nx,ny)        
       complex*16 ykxw(nx),ykyw(ny), ef(nx,ny), den(nx,ny)
      real*8 kx(nx),kx2(nx),ky(ny),ky2(ny)
      real*8 ef2(nx,ny), ef2k(nx,ny),efk(nx,ny)
      real*8 den2k(nx,ny),denk(nx,ny), den2(nx,ny)

       complex*16 ffty1(nx,ny),fftx(nx),ffty(ny)
       complex*16 fftxy(nx), datay(ny),ezf(nx,ny)       
 
      real*8     a0,x, y,ran,th
      real*8     beta,alfx,alfy
      real*8     dt
      integer*4  niter, it
      real*8     ttime, TMAX,DTMX,DTMN
      real*8     nrmx, LONGx, dx
      real*8     nrmy, LONGY, dy
      real*8     nrmxy
      real*8     pi,pi2,c1,c2,c3,c4,c5,c6,c7,c8
      complex*16 ui 
      real*8     plsn
      logical    precor
      integer*4  i,j
      integer*4  i20
      integer*4  idt
      integer*4  iali
      character*4 f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13

      common /kw/ kx,ky/k2/kx2,ky2/ykw/ykxw,ykyw
      common /te/ te 
      common /xe/ xe 
      common /norm/ nrmx,nrmy
      common /step/ dx,dy
      common /dato/ a0, beta,alfx,alfy
      data idt    /0/
      data precor /.TRUE./
      pi=4.d0*datan(1.d0)
      pi2     = 2.d0*pi
      ui      = (0.d0,1.d0)
      dt    = 0.0005d0
      TMAX   = dfloat(nt) 
      
      c1=2.0d0
      c2=0.0132d0 
      c3=1.04d0

      c4=183.62d0
      c5=0.47d0
      c6=17.9d0
      c7=1.0d0
      c8=0.99d0

      ttime = 0.d0
       
      niter = int(float(ndg2)/float(ndg1)+.5d0)  
      it    = 0
                   
      DTMX  = 1.d-5
      DTMN   =DTMX/50.d0       
       alfx   = 0.1d0
       alfy   = 0.1d0
      LONGx = pi2/alfx
      LONGy = pi2/alfy

      beta = 0.1d0 
      a0    = 1.0d0
       
      nrmx     = 1.d0/ dfloat(nx)
      nrmy     = 1.d0/ dfloat(ny)
      nrmxy    =nrmx*nrmy 
      dx      = LONGx/dfloat(nx)
      dy      = LONGy/dfloat(ny)

      f0 = 'ef00'
      call filen(f0, 4)
      open(20, file = f0)


      f1 = 'dn00'
      call filen(f1, 4)
      open(21, file = f1)

      f2 = 'ekz0'
      call filen(f2, 4)
      open(22, file = f2)

      f3 = 'dkz0'
      call filen(f3, 4)
      open(23, file = f3)
      

      f4 = 'tt00'
      call filen(f4, 4)
      open(24, file = f4)

      f5 = 'pm00'
      call filen(f5, 5)
      open(15, file = f5)

      f6 = 'ez00'
      call filen(f6, 4)
      open(25, file = f6)

      f7 = 'dz00'
      call filen(f7, 4)
      open(26, file = f7)

      f8 = 'ex00'
      call filen(f8, 4)
      open(27, file = f8)
    
      f9 = 'dx00'
      call filen(f9, 4)
      open(28, file = f9)

      f10 = 'ekx0'
      call filen(f10, 4)
      open(29, file = f10)

      f11 = 'dkx0'
      call filen(f11, 4)
      open(30, file = f11)
      
      f12 = 'fx00'
      call filen(f12, 4)
      open(31, file = f12)

      f13 = 'fz00'
      call filen(f13, 4)
      open(32, file = f13)                                                                          

      write(15,9000) nt
9000  format(' nt   = ', i6)

      write(15,9001) nx,ny
9001  format(' nx    = ', i6,
     .       ' ny    = ', i6)  

      write(15,9002) ndg, ndg1, ndg2
9002  format(' ndg  = ',i6,
     .       ' ndg1 = ',i6,
     .       ' ndg2 = ',i6)
      
      write(15,9003) LONGx, longy
9003  format(/,' LONGx = ', 1pe10.3,
     .         /,' LONGy = ', 1pe10.3)

      write(15,9004) dt, dx,dy
9004  format(' dt   = ',1pe10.3,
     .       ' dx   = ',1pe10.3,
     .       ' dy   = ',1pe10.3)

      write(15,9005) a0
9005  format(/,' A0   = ',1pe10.3)

      write(15,9006) beta,  alfx,alfy
9006  format(' beta = ',1pe10.3,
     .       ' alfx  = ',1pe10.3,
     .       ' alfy  = ',1pe10.3)
             


      do i = 1, nx/2
         kx(i) = -dfloat(i-1)*pi2/LONGx
         kx(nx-i+1) = dfloat(i)*pi2/LONGx
      enddo
      do j = 1, ny/2
         ky(j) = -dfloat(j-1)*pi2/LONGy
         ky(ny-j+1) = dfloat(j)*pi2/LONGy
       end do

      do i=1,nx
         kx2(i)=kx(i)**2
         ykxw(i) = -ui*(kx2(i)**2)-ui*c1*(kx2(i))
       enddo

       do j=1,ny
         ky2(j)=ky(j)**2
         ykyw(j) = ui*ky2(j)-ui*c2*ky(j)
       
       enddo
 
      call coef(ak,akC,nx,ny,dt)

         do i = 1, nx
         x = dfloat(i-nx/2)*dx
         do j= 1,ny
            y= dfloat(j-ny/2)*dy
      
            te(i,j) = dcmplx(a0*(1.d0+
     .                .1d0*dcos(alfx*x))*(1.d0
     .                +.1d0*dcos(alfy*y)),0.d0)
      
	     td(i,j)=dcmplx(abs(te(i,j))**2, 0.d0)
	     tv(i,j)=(0.d0, 0.d0)
      
         enddo
         enddo
       
      call dcfft2n(1,te,nx,ny,mx,my,e)
      call dcfft2n(1,td,nx,ny,mx,my,d)
      call dcfft2n(1,tv,nx,ny,mx,my,v)
           
      write(15,*) ' Number of Plasmons'

      call intcst(e,nx,ny,plsn,15)
c      write (*,*) plsn
      do j = 1, ny
      do 21 i = 1, nx/2
         efk(nx/2-i+1,j)  = nrmxy*abs(e(i,j))
         denk(nx/2-i+1,j)  = nrmxy*abs(d(i,j))
  21  continue
      do 22 i=nx/2+1, nx
         efk(nx-i+1+nx/2,j)=nrmxy*abs(e(i,j))
         denk(nx-i+1+nx/2,j)=nrmxy*abs(d(i,j))
  22  continue
      end do
      do i = 1, nx
      do 23 j =1, ny/2
         ef2k(i,ny/2-j+1)  = efk(i,j)
         den2k(i,ny/2-j+1)  = denk(i,j) 
  23  continue
      do 24 j=ny/2+1, ny
          ef2k(i,ny-j+1+ny/2)  = efk(i,j)
          den2k(i,ny-j+1+ny/2)  = denk(i,j) 
  24   continue
       enddo

      do i=1,nx
      do j=1,ny
        xe(i,j) = e(i,j)
        xd(i,j) = d(i,j)                                                        
  
      enddo
      enddo


      call dcfft2n (-1,xe,nx,ny,mx,my,xte)
       call dcfft2n(-1,xd,nx,ny,mx,my,xtd)

    
      do i = 1, nx
      do j=1,ny
         ef2(i,j) = (nrmxy*abs(xte(i,j)))**2
         den(i,j) = nrmxy*dble(xtd(i,j)) 
      enddo
      enddo
         
c      call wrtb1(ef2,nx,ny,ndg,20)
c      call wrtb1(den,nx,ny,ndg,21)
                                                     
c         call wrtb1(ef2k,nx,ny,ndg,22)
c          call wrtb1(den2k,nx,ny,ndg,23)

c  nlterm

      do i = 1, nx 
      do j=1,ny
          xe(i,j) = nrmxy*xte(i,j)*nrmxy*
     .             dble(xtd(i,j))
	xd(i,j)=dcmplx((abs(nrmxy*xte(i,j)))**2, 0.d0)
	enddo
        enddo 
      
      call dcfft2n(1,xe,nx,ny,mx,my,te)

       call dcfft2n (1,xd,nx,ny,mx,my,td)
      
      do i = 1, nx
      do j=1,ny
         fe(i,j) = 2.d0*ui*dt*te(i,j)
         fv(i,j)=(dt/(1+(kx2(i))*c3))*((-c4*kx2(i)*ky2(j)
     .          -c5*ky2(j))*d(i,j)+(-c6*kx2(i)*ky2(j)
     .          -c7*kx2(i)*ky2(j)-c8*kx2(i)*(ky(j)**3))*td(i,j))
          
c     .   ky(j))*td(i,j))

      enddo
      enddo
     
      do i = 1, nx
      do j=1,ny   
         b(i,j) = e(i,j)
         e(i,j) = b(i,j)+ui*dt*e(i,j)*((kx2(i)**2)
     .             +c1*kx2(i)-ky2(j)+c2*ky(j))
c	ui*dt*e(i,j)*(kx2(i)+ky(j))
     .             +.5d0*fe(i,j)
         v(i,j)=v(i,j)+fv(i,j)
         d(i,j)=d(i,j)+dt*v(i,j)
      enddo
      end do

c       call intcst(e,nx,ny,plsn,15)
c       write (*,*) plsn

      call transl4(e,ae,b,ab,d,ad,v,av,nx,ny,-1)          
      

100   continue

      if (ttime .le. TMAX) then
         
         

         do i20 = 1, 20
         do iali = 1, 2
         
            if (iali .eq. 2) then
          
               do i = 1, nx
               do j = 1, ny
                  e(i,j)  = ae(i,j)
                  b(i,j)  = ab(i,j)
	          d(i,j)  = ad(i,j)
	          v(i,j)  = av(i,j)
	       end do
               enddo
            endif

 
c  nlterm


       do i = 1, nx
       do j=1,ny
      
         xe(i,j) = e(i,j)
         xd(i,j)=d(i,j)
       enddo
       enddo
      
      call dcfft2n(-1,xe,nx,ny,mx,my,te)
      call dcfft2n (-1,xd,nx,ny,mx,my,td)
      
      
      do i = 1, nx
      do j=1,ny
                                                 
          xe(i,j) = nrmxy*te(i,j)*nrmxy*
     .     dble(td(i,j))
          xd(i,j)=dcmplx((abs(nrmxy*te(i,j)))**2,0.d0)
       end do
       end do

      call dcfft2n(1,xe,nx,ny,mx,my,te)
      call dcfft2n (1,xd,nx,ny,mx,my,td)


      
      do i = 1, nx
      do j=1,ny
         fe(i,j) = 2.d0*ui*dt*te(i,j)
c         fv(i,j)=dt*(-ky2(j)*d(i,j)+(c1*kx2(i)+c2*kx(i)*
c     .   ky(j))*td(i,j))
	fv(i,j)=(dt/(1+(kx2(i))*c3))*((-c4*kx2(i)*ky2(j)
     .          -c5*ky2(j))*d(i,j)+(-c6*kx2(i)*ky2(j)
     .          -c7*kx2(i)*ky2(j)-c8*kx2(i)*(ky(j)**3))*td(i,j))

      enddo
      enddo

           
            do i = 1, nx
            do j = 1, ny
               a(i,j) = e(i,j)
               e(i,j) = b(i,j) + ak(i,j)*e(i,j) + fe(i,j)
               b(i,j) = a(i,j)
	       vc(i,j)=v(i,j)
	       dc(i,j)=d(i,j)
	       v(i,j)=v(i,j)+fv(i,j)
	       d(i,j)=d(i,j)+dt*v(i,j)
	     end do
             enddo
       if (precor) then  
                           
c  nlterm

       do i = 1, nx
       do j=1,ny
         xe(i,j) = e(i,j)
         xd(i,j)=d(i,j)
       enddo
       enddo

      call dcfft2n(-1,xe,nx,ny,mx,my,te)
      call dcfft2n (-1,xd,nx,ny,mx,my,td)

      do i = 1, nx
      do j=1,ny
          xe(i,j) = nrmxy*te(i,j)*nrmxy*
     .     dble(td(i,j))
          xd(i,j)=dcmplx((abs(nrmxy*te(i,j)))**2,0.d0)
       end do
       end do

      call dcfft2n(1,xe,nx,ny,mx,my,te)
      call dcfft2n (1,xd,nx,ny,mx,my,td)
            
      do i = 1, nx
      do j=  1, ny
         feC(i,j) = 2.d0*ui*dt*te(i,j)
c         fvC(i,j)=dt*(-ky2(j)*d(i,j)+(c1*kx2(i)+c2*kx(i)*
c     .   ky(j))*td(i,j))
	fvC(i,j)=(dt/(1+(kx2(i)*c3)))*((-c4*kx2(i)*ky2(j)
     .          -c5*ky2(j))*d(i,j)+(-c6*kx2(i)*ky2(j)
     .          -c7*kx2(i)*ky2(j)-c8*kx2(i)*(ky(j)**3))*td(i,j))



      end do
      end do



            do i = 1, nx
            do j= 1,ny
               e(i,j) = b(i,j)  + 
     .                    .5d0*akC(i,j)*(b(i,j) + e(i,j)) +
     .                    .25d0*(fe(i,j) + feC(i,j))
               
	       v(i,j)=vc(i,j)+.5d0*(fv(i,j)+fvC(i,j))
	       d(i,j)=dc(i,j)+dt*v(i,j)
	     end do
             enddo

            endif           


            if (iali .eq. 1) then
               do i = 1, nx
               do j=1,ny
                  be(i,j) = e(i,j)
                  bb(i,j) = b(i,j)
		  bd(i,j)=d(i,j)
		  bv(i,j)=v(i,j)
		 end do
                 end do
            else
               do i = 1, nx
               do j = 1, ny
                  ae(i,j) = e(i,j)
                  ab(i,j) = b(i,j)
	          ad(i,j) = d(i,j)
	          av(i,j) = v(i,j)
	       end do
               end do 

            endif
c  end iali: 1, 2
         enddo
         
            call antial4(e,be,ae,b,bb,ab,d,bd,ad,v,bv,av,nx,ny)  
            call transl4(e,ae,b,ab,d,ad,v,av,nx,ny,-1)
             
c  end i20: 1, 20
                   
         enddo

            call intcst(e,nx,ny,plsn,15)
         
         ttime = ttime + 20.d0*dt
         
           

           if (int(niter*ttime) .eq. it) then
            it = it + 1
c            write (*,*) ttime       
            write (24,*)  ttime

            do j = 1,ny
            do 31 i = 1, nx/2
               efk(nx/2-i+1,j)  = nrmxy*abs(e(i,j))
                denk(nx/2-i+1,j)  = nrmxy*abs(d(i,j))
               
  31       continue
          do 32 i=nx/2+1, nx
               efk(nx-i+1+nx/2,j)  = nrmxy*abs(e(i,j))
               denk(nx-i+1+nx/2,j)  = nrmxy*abs(d(i,j))
   32    continue
            enddo
                                                                                
            do i = 1, nx
            do 33 j =1, ny/2
               ef2k(i,ny/2-j+1)  = efk(i,j)
              den2k(i,ny/2-j+1)  =denk(i,j)
   33     continue
           do 34 j=ny/2+1, ny
                ef2k(i,ny-j+1+ny/2)  = efk(i,j)
                 den2k(i,ny-j+1+ny/2)  =denk(i,j)
  34        continue
            enddo

            do i =1,nx
            do j =1,ny
               xe(i,j)  = e(i,j)
               xb(i,j)  = b(i,j)
	       xd(i,j)  = d(i,j)
            end do
            end do
            call dcfft2n(-1,xe,nx,ny,mx,my,xte)
            call dcfft2n(-1,xd,nx,ny,mx,my,xtd)             
             do i = 1, nx
             do j= 1,ny
               ef2(i,j) = (nrmxy*abs(xte(i,j)))**2
              den(i,j)=  nrmxy*(xtd(i,j))
	      end do
              end do
c              do i = 1, nx
c             do j= 1,ny
c               write(20,*) real(nrmxy*xte(i,j)), imag(nrmxy*xte(i,j))
c               write(21,*) real(nrmxy*xtd(i,j)), imag(nrmxy*xtd(i,j))
c              end do
c              end do

             do j=1,ny
             write(25,*) nrmxy*abs(xte(nx/2,j))
             write(26,*)   nrmxy*(xtd(nx/2,j))
             end do
             
             do i=1,nx
             write(27,*) nrmxy*abs(xte(i,ny/2))
             write(28,*)   nrmxy*(xtd(i,ny/2))
             end do


            do j=1,ny
              write(22,*) ef2k(nx/2,j)
              write(23,*) den2k(nx/2,j)
            end do

           do i=1,nx
              write(29,*) ef2k(i,ny/2)
              write(30,*) den2k(i,ny/2)
            end do      

            write(31,*) ef2k(nx/2,ny/2),ef2k((nx/2)+1,ny/2),
     .                  ef2k((nx/2)+2,ny/2)

            write(32,*) ef2k(nx/2,ny/2),ef2k(nx/2,(ny/2)+1),
     .                  ef2k(nx/2,(ny/2)+2)

            write (*,*) it, plsn
            write(15,*) it, plsn
            
       call wrtb1(ef2,nx,ny,1,20)
       call wrtb1(den,nx,ny,1,21)

c       call wrtb1(ef2k,nx,ny,1,22)
c       call wrtb1(den2k,nx,ny,1,23)      
            
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
      close(27)
      close(28)
      close(15)
     
  
    
      
      write(*,*) '*** end ***'
      stop
      end

      subroutine wrtb1(tb,nx,ny,ndiag,uni)
      implicit none
      
      integer uni, nx,ny, ndiag, i,j
      real*8 tb(nx,ny)

      do i = 1, nx, ndiag
      do j = 1, ny, ndiag

         write(uni, 9000) tb(i,j)
      enddo
      enddo

      return
9000  format(1pe15.5)
      end


      subroutine intcst(tb1,nx,ny,plas,uni)
      implicit none
      
      integer uni, nx,ny
      complex*16 tb1(nx,ny)
      
      real*8 plas
      real*8 nrmx,nrmy,nrm2x,nrm2y
      real*8 t1
      integer*4 i,j
      
      common /norm/ nrmx,nrmy

      nrm2x = nrmx**2
      nrm2y = nrmy**2

      plas = 0.d0
      do i = 1, nx
      do j=1,ny
         t1 = abs(tb1(i,j))**2
         plas = plas + t1
      enddo
      enddo

      plas = nrm2x*nrm2y*plas
    
      return
      end

      subroutine coef(ak,akC,nx,ny,dt)
      implicit none
      
      integer*4 nx,ny
      complex*16 ykxw(1),ykyw(1)
      complex*16 ak(nx,ny),akC(nx,ny)
      
      real*8 dt

      complex*16 ui
      integer*4 i,j
      
      real*8 kx(1), kx2(1),ky(1)
      
      common /kw/ kx,ky/k2/ kx2
      common /ykw/ykxw,ykyw  
        
      ui = (0.d0,1.d0)
      do i = 1, nx
      do j=   1,ny
         ak(i,j)  =cdexp(-(ykyw(j)*dt+ykxw(i)*dt))
     .           -cdexp(ykyw(j)*dt+ykxw(i)*dt)
        akc(i,j) =cdexp(-.5d0*(ykyw(j)*dt+ykxw(i)*dt))
     .             -cdexp(.5d0*(ykyw(j)*dt+ykxw(i)*dt))
         
      enddo
      enddo
         

      return
      end

      

     

      subroutine antial4(h1,f1,g1,h2,f2,g2,h3,f3,g3,h4,f4,g4,nx,ny)
      implicit none
      
      integer*4 nx,ny
      complex*16 h1(nx,ny), f1(nx,ny), g1(nx,ny)
      complex*16 h2(nx,ny),f2(nx,ny),g2(nx,ny)
      complex*16 h3(nx,ny),f3(nx,ny),g3(nx,ny)
      complex*16 h4(nx,ny),f4(nx,ny),g4(nx,ny)
      real*8 kx(1),ky(1)
      real*8 dx,dy
      real*8 dx2,dy2

      real*8 a,ay
      complex*16 c,cy
      complex*16 ui
      integer*4 i,j
      
      common /kw/ kx,ky
      common /step/ dx,dy

      ui = (0.d0,1.d0)
      
      dx2 = .5d0*dx
      dy2 = .5d0*dy

      do i = 1, nx   
      do j= 1,ny
     
         a = (dx2*kx(i)+dy2*ky(j))
         c = dcmplx( dcos(a), -dsin(a) )
         h1(i,j) = .5d0*( f1(i,j) + c*g1(i,j) )
         h2(i,j) = .5d0*( f2(i,j) + c*g2(i,j) )
         h3(i,j)=.5d0*(f3(i,j)+c*g3(i,j))
         h4(i,j)=.5d0*(f4(i,j)+c*g4(i,j))
      enddo
      enddo
      return
      end


      subroutine transl4(h1,f1,h2,f2,h3,f3,h4,f4,nx,ny,isgn)

      implicit none
      
      integer*4 nx,ny, isgn
      complex*16 h1(nx,ny), f1(nx,ny) 
      complex*16 h2(nx,ny),f2(nx,ny)  
      complex*16 h3(nx,ny), f3(nx,ny)
      complex*16 h4(nx,ny), f4(nx,ny)
      real*8 kx(1),ky(1),kx2(1)
      real*8 dx,dy
      real*8 dx2,dy2
      real*8 sg, a,a1
      complex*16 ui
      complex*16 c,c1
      integer*4 i,j
      
      common /kw/ kx,ky /k2/ kx2
      common /step/ dx,dy
      
      ui = (0.d0,1.d0)
      dx2 = .5d0*dx
      dy2= .5d0*dy
      sg  = dfloat(isgn)
      
      do i = 1, nx
      do j = 1, ny
          
          a =sg*( dx2*kx(i)+dy2*ky(j))
         c = dcmplx(dcos(a),- dsin(a))
         f1(i,j) = c*h1(i,j)
         f2(i,j)  =c*h2(i,j)

       	f3(i,j)=c*h3(i,j)
	f4(i,j)=c*h4(i,j)
        enddo
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


      subroutine dcfft2n(isign,DATA2,ngx,ngy,mx,my,FFT2)
      integer*4  ngy, ngx,mx,my
      complex*16 DATA2
      complex*16 DATAY
      complex*16 FFT2
      dimension  DATA2(ngx,ngy)
      dimension DATAY(ngy)
      dimension FFTY1(ngx,ngy)
      dimension FFTX(ngx),FFTY(ngy),FFTXY(ngx)
      dimension FFT2(ngx,ngy)
      complex*16 FFTY
      complex*16 FFTY1
      complex*16 FFTX
      complex*16 FFTXY
      integer*4  isign
      do ix=1,ngx

      do 11 iy =1,ngy
         DATAY(iy)=DATA2(ix,iy)
 11   continue

      call dcfftn(isign,DATAY,FFTY,my)
       do 12 iy1=1,ngy
       FFTY1(ix,iy1)=FFTY(iy1)
 12   continue

      end do
      do iy=1,ngy
        do 13 ix=1,ngx
        FFTX(ix)=FFTY1(ix,iy)
  13  continue

      call dcfftn(isign,FFTX,FFTXY,mx)

      do 14 ix1=1,ngx
      FFT2(ix1,iy)=FFTXY(ix1)

 14   continue

      end do
      return
      end

 











