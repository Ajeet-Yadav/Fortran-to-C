
      implicit none
      
      integer*4 nt, mx,my,nx,ny,nz,mz
      integer*4 ndg, ndg1, ndg2
      integer*4 ntd,sd
      integer*4  niter, it
      parameter (nt = 15, mx =5, nx = 2**mx)
      parameter ( my =5, ny = 2**my)
      parameter (mz=5, nz=2**mz)

      parameter (ndg =2 , ndg1 = 1, ndg2 =1)
      parameter (ntd = nt*ndg2/ndg1)
       
      complex*16 a(nx,ny,nz),a2(nx,ny,nz)
      complex*16 e(nx,ny,nz),be(nx,ny,nz),ae(nx,ny,nz)
      complex*16 b(nx,ny,nz),bb(nx,ny,nz),ab(nx,ny,nz)
      complex*16 d(nx,ny,nz),ad(nx,ny,nz), bd(nx,ny,nz)
      complex*16 b2(nx,ny,nz),ab2(nx,ny,nz), bb2(nx,ny,nz)
      complex*16 fe(nx,ny,nz),fd(nx,ny,nz)
      complex*16 fe(nx,ny,nz),fd1(nx,ny,nz)
      complex*16 fe(nx,ny,nz),fd2(nx,ny,nz)
      complex*16 ak(nx,ny,nz), ak2(nx,ny,nz)
      complex*16 akC(nx,ny,nz),akC2(nx,ny,nz)
      complex*16 feC(nx,ny,nz),fdC(nx,ny,nz)
      complex*16 xte(nx,ny,nz),xe(nx,ny,nz)
      complex*16 xtd(nx,ny,nz),xd(nx,ny,nz)
      complex*16 te(nx,ny,nz),td(nx,ny,nz)
      complex*16 tb(nx,ny,nz),tb1(nx,ny,nz)
      complex*16 xb(nx,ny,nz),xb1(nx,ny,nz),xb2(nx,ny,nz)
      complex*16 ykxw(nx),ykyw(ny),ykzw(nz)
      complex*16 ykxw2(nx),ykyw2(ny),ykzw2(nz)
      complex*16 ffty1(nx,ny,nz),fftx(nx),ffty(ny),fftz(nz)
      complex*16 fftxy(nx), datay(ny),ezf(nx,ny,nz)
      complex*16 den(nx,ny,nz),ef(nx,ny,nz)
      real*8     kx(nx),kx2(nx),ky(ny),ky2(ny),kz(nz),kz2(nz)
      real*8     den2k(nx,ny,nz),denck(nx,ny,nz)
      real*8     den2(nx,ny,nz),ef2(nx,ny,nz)
      real*8     efk(nx,ny,nz),denk(nx,ny,nz),den3k(nx,ny,nz)
      real*8     ef2k(nx,ny,nz),efck(nx,ny,nz),ef3k(nx,ny,nz)
      real*8     c1,c2,c3,c4,c5,c6,c7,c8,g,alf
      real*8     a0,x,y,z
      real*8     beta,gam,alfx,alfy,alfz
      real*8     plsn,wn
      real*8     vel
      real*8     dt
      real*8     ttime, TMAX,DTMX,DTMN
      real*8     nrmx, LONGx, dx
      real*8     nrmy, LONGY, dy
      real*8     nrmz, LONGz, dz
      real*8     nrmxyz
      real*8     pi,pi2,th
      complex*16 ui
      
      logical    precor
      integer*4  i,j,k
      integer*4  i20
      integer*4  idt
      integer*4  iali
                     
    

      character*4 f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10
      character*4 f11,f12,f13,f14,f15,f16,f17,f18,f19,f20
      common /kw/ kx,ky,kz/k2/kx2,ky2,kz2/ykw/ykxw,ykyw,ykzw
      common /ykw2/ ykxw2,ykyw2,ykzw2
      common /te/ te /td/ td / fd1
      common /xe/ xe /xd/ xd
      common /norm/ nrmx,nrmy,nrmz
      common /step/ dx,dy,dz
      common /dato/ a0, beta,gam,alf,alfx,alfy,alfz
      data precor  /.TRUE./

      data precor /.TRUE./
      pi=4.d0*datan(1.d0)
      pi2     = 2.d0*pi
      ui      = (0.d0,1.d0)
      dt    = 0.0005d0
      TMAX   = dfloat(nt)
      ttime = 0.d0
      niter = int(float(ndg2)/float(ndg1)+.5d0)  
      it    = 0
      DTMX  = 1.d-5
      DTMN   =DTMX/50.d0
        c1=1.04d0
c        c2=1.0d0
c        c3=1.0d0
        c4=0.454d0
        c5=0.00175d0
        c6=14.181d0
        c7=51.48d0
        c8=0.d0

       alfx   = 0.2d0
       alfy   = 0.2d0
       alfz = 0.2d0
      LONGx = pi2/alfx
      LONGy = pi2/alfy
      LONGZ=  pi2/alfz
      beta  = 0.1d0
      gam = 0.d0
      a0 =1.2d0
      vel = 0.d0
      nrmx     = 1.d0/ dfloat(nx)
      nrmy     = 1.d0/ dfloat(ny)
      nrmz     = 1.d0/ dfloat(nz)
      nrmxyz    =nrmx*nrmy*nrmz 
      dx      = LONGx/dfloat(nx)
      dy      = LONGy/dfloat(ny)
      dz      = LONGZ/dfloat(nz)
      
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

      f14 = 'ne00'
      call filen(f14,4)
      open(35, file = f14)

      f15 = 'nh00'
      call filen(f15,4)
      open(36, file = f15)

      f16 = 'np00'
      call filen(f16,4)
      open(11, file = f16)
      
      f17 = 'ns00'
      call filen(f17,4)
      open(12,file = f17)

      f18 = 'sg00'
      call filen(f18,4)
      open(33,file = f18)

      f19 = 'nd00'
      call filen(f19,4)
      open(34, file = f19)
      
      fd1 = 'nfd00'
      call filen(f20,4)
      open(35, file = f20)
      
      fd2 = 'nfd100'
      call filen(f21,4)
      open(35, file = f21)
      
                       


      write(15,9000) nt
9000  format(' nt   = ', i6)

      write(15,9001) nx,ny,nz
9001  format(' nx    = ', i6,
     .       ' ny    = ', i6,
     .       ' nz    = ', i6)  

      write(15,9002) ndg, ndg1, ndg2
9002  format(' ndg  = ',i6,
     .       ' ndg1 = ',i6,
     .       ' ndg2 = ',i6)
      
      write(15,9003) LONGx, longy,LONGZ
9003  format(/,' LONGx = ', 1pe10.3,
     .         /,' LONGy = ', 1pe10.3,
     .         /,'LONGZ  = ', 1pe10.3)

      write(15,9004) dt, dx,dy,dz
9004  format(' dt   = ',1pe10.3,
     .       ' dx   = ',1pe10.3,
     .       ' dy   = ',1pe10.3,
     .       ' dz   = ',1pe10.3)
      write(15,9005) a0
9005  format(/,' A0   = ',1pe10.3)

      write(15,9006) beta,  alfx,alfy,alfz
9006  format(' beta = ',1pe10.3,
     .       ' alfx  = ',1pe10.3,
     .       ' alfy  = ',1pe10.3,
     .       ' alfz  = ',1pe10.3)
             


      do i = 1, nx/2
         kx(i) = -dfloat(i-1)*pi2/LONGx
         kx(nx-i+1) = dfloat(i)*pi2/LONGx
      enddo
      do j = 1, ny/2
         ky(j) = -dfloat(j-1)*pi2/LONGy
         ky(ny-j+1) = dfloat(j)*pi2/LONGy
       end do
      do k=1, nz/2
         kz(k)= -dfloat(k-1)*pi2/LONGZ
         Kz(nz-k+1) = dfloat(k)*pi2/LONGZ
      end do 
      do i=1,nx
         kx2(i)=kx(i)**2
         ykxw(i) = ui*(kx2(i))
         ykxw2(i)= ui*(c6*kx2(i))

       enddo

       do j=1,ny
         ky2(j)=ky(j)**2
         ykyw(j) = ui*(ky2(j))
         ykyw2(j)=ui*(c8*ky(j))

       enddo
       do k=1,nz
         kz2(k)=kz(k)**2
         ykzw(k)= ui*kz(k)
         ykzw2(k)=ui*c4*kz(k)+ui*c5*kz2(k)

       enddo
          
 
      call coef(ak,akC,nx,ny,nz,dt)
      call coef2 (ak2,akC2,nx,ny,nz,dt)

         do i = 1, nx
         x = dfloat(i-nx/2)*dx
         do j= 1,ny
         y= dfloat(j-ny/2)*dy
         do k= 1,nz
         z= dfloat(k-nz/2)*dz
         te(i,j,k) = dcmplx(a0*(4/(pi**0.5))*((x**2+y**2)**1.5      
     .                 /((1+z**2)**2)*dexp(-((x**2+y**2)
     .                /(1+z**2)))),0.d0)
      
         td(i,j,k)=dcmplx(0.08*(1.d0+
     .                .05d0*dcos(alfx*x))*(1.d0
     .                +.05d0*dcos(alfy*y))*(1.d0
     .                +.05d0*dcos(alfz*z)),0.d0)


         enddo
         enddo
         enddo

      call dcfft3n(1,te,nx,ny,nz,mx,my,mz,e)
      call dcfft3n(1,td,nx,ny,nz,mx,my,mz,d)
      write(15,*) ' Number of Plasmons'

      call intcst(e,nx,ny,nz,plsn,15)
      write (*,*) plsn

      do k=1,nz
      do j = 1, ny
      do 21 i = 1, nx/2
         efk(nx/2-i+1,j,k)  = nrmxyz*abs(e(i,j,k))
         denk(nx/2-i+1,j,k)  = nrmxyz*abs(d(i,j,k))

  21  continue
      do 22 i=nx/2+1, nx
         efk(nx-i+1+nx/2,j,k)=nrmxyz*abs(e(i,j,k))
         denk(nx-i+1+nx/2,j,k)=nrmxyz*abs(d(i,j,k))
  22  continue
      end do
      end do
      write(*,*) ' Number of Plasmons'

      do  k=1,nz
      do  i = 1, nx
      do 23 j =1, ny/2
         ef2k(i,ny/2-j+1,k)  = efk(i,j,k)
         den2k(i,ny/2-j+1,k)  = denk(i,j,k)

  23  continue
      do 24 j=ny/2+1, ny
          ef2k(i,ny-j+1+ny/2,k)  = efk(i,j,k)
          den2k(i,ny-j+1+ny/2,k)  = denk(i,j,k)

  24   continue
       enddo
       enddo

      do i=1,nx
      do j=1,ny
      do 25 k=1,nz/2
        ef3k(i,j,nz/2-k+1)  = ef2k(i,j,k)
        den3k(i,j,nz/2-k+1) = den2k(i,j,k)

 25   continue
      do 26 k=nz/2+1,nz
           ef3k(i,j,nz-k+1+nz/2)= ef2k(i,j,k)
           den3k(i,j,nz-k+1+nz/2) = den2k(i,j,k)

 26   continue
      enddo
      enddo



   
      do i=1,nx
      do j=1,ny
      do k = 1,nz
        xe(i,j,k)  = e(i,j,k)
        xd(i,j,k) = d(i,j,k)    
      enddo
      enddo
      enddo
      call dcfft3n (-1,xe,nx,ny,nz,mx,my,mz,xte)
      call dcfft3n(-1,xd,nx,ny,nz,mx,my,mz,xtd)


    
      do i = 1, nx
      do j=1,ny
      do k=1,nz
         ef2(i,j,k) = (nrmxyz*abs(xte(i,j,k)))**2
          ef(i,j,k)= nrmxyz*(xte(i,j,k))
          den(i,j,k) = nrmxyz*(xtd(i,j,k))

      enddo
      enddo
      enddo


     

       write(*,*) ' Number of Plasmons'
c  nlterm

      do i = 1, nx 
      do j=1,ny
      do k=1,nz
          xe(i,j,k) = nrmxyz*xte(i,j,k)*(dcmplx((abs(nrmxyz*
     .             xte(i,j,k)))**2,0.d0))

          xd(i,j,k) = nrmxyz*xtd(i,j,k)*(dcmplx((abs(nrmxyz*
     .             xte(i,j,k)))**2,0.d0))

      enddo
      enddo
      enddo
      
      call dcfft3n(1,xe,nx,ny,nz,mx,my,mz,te)
       call dcfft3n (1,xd,nx,ny,nz,mx,my,mz,td)

      do i = 1, nx
      do j=1,ny
      do k=1,nz
       fe(i,j,k) = ui*c1*dt*te(i,j,k)
       fd(i,j,k)=ui*c7*dt*td(i,j,k)

      enddo
      enddo
      enddo
  
     
      do i = 1, nx
      do j=1,ny 
      do k=1,nz  
         b(i,j,k) = e(i,j,k)
         e(i,j,k) = b(i,j,k)
     .             -ui*dt*e(i,j,k)*(kx2(i)+ky2(j)+kz(k))
     .             +fe(i,j,k)

         b2(i,j,k)=d(i,j,k)
         d(i,j,k)=b2(i,j,k)-ui*dt*d(i,j,k)*(c6*kx2(i)
     .             +c5*kz2(j)+c4*kz(k))+fd(i,j,k)

      enddo
      enddo
      enddo

       call intcst(e,nx,ny,nz,plsn,15)
       write (*,*)plsn

       call transl4(e,ae,b,ab,d,ad,b2,ab2,nx,ny,nz,-1)


 100   continue

      if (ttime .le. TMAX) then
         
         

         do i20 = 1, 20
         do iali = 1, 2

         
            if (iali .eq. 2) then
          
               do i = 1, nx
               do j = 1, ny
               do k = 1, nz
                  e(i,j,k)  = ae(i,j,k)
                  b(i,j,k)  = ab(i,j,k)
                  d(i,j,k)  = ad(i,j,k)
                  b2(i,j,k)  = ab2(i,j,k)

               enddo
               enddo
               enddo
            endif

 
c  nlterm



       do i = 1, nx
       do j=1,ny
       do k=1,nz
         xe(i,j,k) = e(i,j,k)
         xd(i,j,k)=d(i,j,k)

       enddo   
       enddo
       enddo
      
       call dcfft3n(-1,xe,nx,ny,nz,mx,my,mz,te)
       call dcfft3n(-1,xd,nx,ny,nz,mx,my,mz,td)
      
       do i = 1, nx
       do j=1,ny
       do k=1,nz
          xe(i,j,k) = nrmxyz*te(i,j,k)*(dcmplx((abs(nrmxyz*
     .             te(i,j,k)))**2,0.d0))

          xd(i,j,k) = nrmxyz*td(i,j,k)*(dcmplx((abs(nrmxyz*
     .             te(i,j,k)))**2,0.d0))

       enddo
       enddo
       enddo

      call dcfft3n(1,xe,nx,ny,nz,mx,my,mz,te)
      call dcfft3n (1,xd,nx,ny,nz,mx,my,mz,td)


      
      do i = 1, nx
      do j=1,ny
      do k=1,nz
         fe(i,j,k) = ui*c1*dt*te(i,j,k)
         fd(i,j,k)=ui*c7*dt*td(i,j,k)

         
      enddo
      enddo
      enddo

           
            do i = 1, nx
            do j = 1, ny
            do k = 1, nz
               a(i,j,k) = e(i,j,k)
               e(i,j,k) = b(i,j,k) + ak(i,j,k)*e(i,j,k) +fe(i,j,k)
               b(i,j,k) = a(i,j,k)
               
               a2(i,j,k)=d(i,j,k)
               d(i,j,k)=b2(i,j,k)+ ak2(i,j,k)*d(i,j,k)+ fd(i,j,k)
               b2(i,j,k)=a2(i,j,k)

            enddo
            enddo
            enddo


       if (precor) then       

       do i = 1, nx
       do j=1,ny
       do k=1,nz
         xe(i,j,k) = e(i,j,k)
         xd(i,j,k)=d(i,j,k)

       enddo 
       enddo
       enddo
                  
      call dcfft3n(-1,xe,nx,ny,nz,mx,my,mz,te)
      call dcfft3n(-1,xd,nx,ny,nz,mx,my,mz,td)
   
            
      do i = 1, nx
      do j=1,ny
      do k=1,nz
        xe(i,j,k) = nrmxyz*te(i,j,k)*(dcmplx((abs(nrmxyz*
     .             te(i,j,k)))**2,0.d0))
       
        xd(i,j,k) = nrmxyz*td(i,j,k)*(dcmplx((abs(nrmxyz*
     .             te(i,j,k)))**2,0.d0))
       enddo
       enddo
       enddo

      call dcfft3n(1,xe,nx,ny,nz,mx,my,mz,te)
      call dcfft3n(1,xd,nx,ny,nz,mx,my,mz,td)


       
            
      do i = 1, nx
      do j=1,ny
      do k=1,nz
            feC(i,j,k)=ui*c1*dt*te(i,j,k)
            fdC(i,j,k)=ui*c7*dt*td(i,j,k)

      enddo
      enddo
      enddo



            do i = 1, nx
            do j= 1,ny
            do k= 1,nz
               e(i,j,k) = b(i,j,k)  + 
     .                    0.50d0*akC(i,j,k)*(b(i,j,k) + e(i,j,k)) 
     .                    +0.5d0*(fe(i,j,k) + feC(i,j,k))
               d(i,j,k)=b2(i,j,k) +
     .                    0.50d0*akC2(i,j,k)*(b2(i,j,k)+d(i,j,k))
     .                   +0.5d0*(fd(i,j,k)+fdC(i,j,k))

            enddo
            enddo 
            enddo
            endif           


            if (iali .eq. 1) then
               do i = 1, nx
               do j=1,ny
               do k=1,nz
                  be(i,j,k) = e(i,j,k)
                  bb(i,j,k) = b(i,j,k)
                  bd(i,j,k) = d(i,j,k)
                  bb2(i,j,k)=b2(i,j,k)

               enddo
               enddo
               enddo
            else
               do i = 1, nx
               do j = 1, ny
               do k = 1, nz
                  ae(i,j,k) = e(i,j,k)
                  ab(i,j,k) = b(i,j,k)
                  ad(i,j,k) = d(i,j,k)
                  ab2(i,j,k)=b2(i,j,k)
               enddo
               enddo
               enddo
            endif
c  end iali: 1, 2
         enddo
         
            call antial4(e,be,ae,b,bb,ab,d,bd,ad,b2,bb2,ab2,nx,ny,nz)
            call transl4(e,ae,b,ab,d,ad,b2,ab2,nx,ny,nz,-1)

c  end i20: 1, 20
                   
         enddo

            call intcst(e,nx,ny,nz,plsn,15)
         
         ttime = ttime + 20.d0*dt
         
           

           if (int(niter*ttime) .eq. it) then
            it = it + 1
            write (*,*) ttime       
          
            do k = 1,nz
            do j = 1,ny
            do 31 i = 1, nx/2
               efk(nx/2-i+1,j,k)  = nrmxyz*abs(e(i,j,k))
               denk(nx/2-i+1,j,k)  = nrmxyz*abs(d(i,j,k))

  31       continue
          do 32 i=nx/2+1, nx
               efk(nx-i+1+nx/2,j,k)  = nrmxyz*abs(e(i,j,k))
               denk(nx-i+1+nx/2,j,k)  = nrmxyz*abs(d(i,j,k))

   32    continue
            enddo
            enddo
 
                                                                                             
            do k = 1, nz
            do i = 1, nx
            do 33 j =1, ny/2
               ef2k(i,ny/2-j+1,k)  = efk(i,j,k)
               den2k(i,ny/2-j+1,k)  =denk(i,j,k)

   33     continue
           do 34 j=ny/2+1, ny
                ef2k(i,ny-j+1+ny/2,k)  = efk(i,j,k)
                den2k(i,ny-j+1+ny/2,k)  =denk(i,j,k)

  34        continue
            enddo
            enddo
            


            do i =1,nx
            do j =1,ny
            do 35 k =1,nz/2
            ef3k(i,j,nz/2-k+1)=ef2k(i,j,k)
            den3k(i,j,nz/2-k+1)=den2k(i,j,k)

  35        continue
             Do 36 k=nz/2+1,nz
             ef3k(i,j,nz-k+nz/2+1)=ef2k(i,j,k)
             den3k(i,j,nz-k+1+nz/2)=den2k(i,j,k)

  36        continue
            Enddo
            Enddo

               
             do i = 1, nx
             do j= 1,ny
             do k=1,nz
	       xe(i,j,k)  = e(i,j,k)
               xb(i,j,k)  = b(i,j,k)
              xd(i,j,k)  = d(i,j,k)

            enddo
            enddo
            enddo
            call dcfft3n(-1,xe,nx,ny,nz,mx,my,mz,xte)
            call dcfft3n(-1,xd,nx,ny,nz,mx,my,mz,xtd)
             do i = 1, nx
             do j= 1,ny
             do k=1,nz
               ef2(i,j,k) = (nrmxyz*abs(xte(i,j,k)))**2
               ef(i,j,k)=nrmxyz*xte(i,j,k)
               den(i,j,k)=  nrmxyz*(xtd(i,j,k))
               den2(i,j,k)=(nrmxyz*abs(xtd(i,j,k)))**2

            enddo
            enddo
            enddo
            
           
             do j=1,ny
             do k=1,nz
             write(25,*) nrmxyz*abs(xte(nx/2,j,k))
             write(26,*)   nrmxyz*(xtd(nx/2,j,k))
             end do
             end do

             do i=1,nx
             do k=1,nz
             write(27,*) nrmxyz*abs(xte(i,ny/2,k))
             write(28,*) nrmxyz*(xtd(i,ny/2,k))
             end do
             end do
              
               do i=1,nx
               do j=1,ny
             write(35,*) nrmxyz*abs(xte(i,j,nz/2))
             write(36,*)   nrmxyz*(xtd(i,j,nz/2))
             end do
             end do

             do k=1,nz
              write(22,*) ef3k(nx/2,ny/2,k)
              write(23,*) den3k(nx/2,ny/2,k)
            end do


           do i=1,nx
              write(29,*) ef3k(i,ny/2,nz/2)
              write(30,*) den3k(i,ny/2,nz/2)
           end do



            do j=1,ny
              write(11,*) ef3k(nx/2,j,nz/2)
              write(12,*) den3k(nx/2,j,nz/2)
            end do

             do i = 1, nx
             do j= 1,ny
             do k = 1,nz
               write(24,2002) ef(i,j,k)
 2002         format(1pe15.5,1pe15.5)
              end do
              end do
              end do

            write(31,*) ef3k(nx/2,ny/2,nz/2),ef3k((nx/2)+1,ny/2,nz/2),
     .                  ef3k((nx/2)+2,ny/2,nz/2)
            write(32,*) ef3k(nx/2,ny/2,nz/2),ef3k(nx/2,(ny/2)+1,nz/2),
     .                  ef3k(nx/2,(ny/2)+2,nz/2)
            write(33,*) ef3k(nx/2,ny/2,nz/2),ef3k(nx/2,ny/2,(nz/2)+1),
     .                  ef3k(nx/2,ny/2,(nz/2)+2)
            write (*,*) it, plsn
            write(15,*) it, plsn

          call wrtb1(ef2,nx,ny,nz,1,20)
          call wrtb1(den2,nx,ny,nz,1,21)
             do i = 1, nx
             do j= 1,ny
             do k = 1,nz
               write(34,2003) den(i,j,k)
 2003         format(1pe15.5,1pe15.5)
              end do
              end do
              end do


           

          
         

      call wrtb1(ef2,nx,ny,nz,1,20)
      call wrtb1(ef3k,nx,ny,nz,1,21)
            
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
      close(11)
      close(12)
      close(35)
      close(36)
      close(33)
      close(31)
      close(29)
      close(30)
      close(32)


    
      
      write(*,*) '*** end ***'
      stop
      end

      subroutine wrtb1(tb,nx,ny,nz,ndiag,uni)
      implicit none
      
      integer uni, nx,ny,nz, ndiag, i,j,K
      real*8 tb(nx,ny,nz)

      do i = 1, nx, ndiag
      do j = 1, ny, ndiag
      do k = 1, nz, ndiag

         write(uni, 9000) tb(i,j,k)
      enddo
      enddo
      enddo

      return
9000  format(1pe15.5)
      end

      subroutine intcst(tb1,nx,ny,nz,plas,uni)
      implicit none
      
      integer uni, nx,ny,nz
      complex*16 tb1(nx,ny,nz)
      
      real*8 plas
      real*8 nrmx,nrmy,nrmz,nrm2x,nrm2y,nrm2z
      real*8 t1
      integer*4 i,j,k
      
      common /norm/ nrmx,nrmy,nrmz

      nrm2x = nrmx**2
      nrm2y = nrmy**2
      nrm2z = nrmz**2
      plas = 0.d0
      do i = 1, nx
      do j=1,ny
      do k=1,nz
         t1 = abs(tb1(i,j,k))**2
         plas = plas + t1
      enddo
      enddo
      enddo
      plas = nrm2x*nrm2y*nrm2z*plas
    
      return
      end

      subroutine coef(ak,akC,nx,ny,nz,dt)
      implicit none
      
      integer*4 nx,ny,nz
      complex*16 ykxw(1),ykyw(1),ykzw(1)
      complex*16 ak(nx,ny,nz),akC(nx,ny,nz)
c      complex*16 ykw(n)
c      complex*16 ak(n),akC(n)
      real*8 dt
      real*8 yi(1,1)

      complex*16 ui
      integer*4 i,j,k
      
      real*8 kx(1), kx2(1),ky(1),ky2(1),kz2(1),kz(1)
      
      common /kw/ kx,ky,kz/k2/ kx2,ky2,kz2
      common /ykw/ykxw,ykyw,ykzw 
c       common/ykw/ykw 
        
      ui = (0.d0,1.d0)
      do i = 1, nx
      do j=   1,ny
      do k= 1,nz
         ak(i,j,k)  =cdexp(-(ykyw(j)*dt+ykxw(i)*dt+ykzw(k)*dt))
     .           -cdexp(ykyw(j)*dt+ykxw(i)*dt+ykzw(k)*dt)
        akC(i,j,k) =cdexp(-.5d0*(ykyw(j)*dt+ykxw(i)*dt+ykzw(k)*dt))
     .             -cdexp(.5d0*(ykyw(j)*dt+ykxw(i)*dt+ykzw(k)*dt))

         
      enddo
      enddo
      enddo
      
      return
      end

      
      subroutine coef2(ak2,akC2,nx,ny,nz,dt)
      implicit none

      integer*4 nx,ny,nz
      complex*16 ykxw2(1),ykyw2(1),ykzw2(1)
      complex*16 ak2(nx,ny,nz),akC2(nx,ny,nz)

      real*8 dt

      complex*16 ui
      integer*4 i,j,k


      real*8 kx(1), kx2(1),ky(1),ky2(1),kz2(1),kz(1)


      common /kw/ kx,ky,kz/k2/ kx2,ky2,kz2
      common /ykw2/ ykxw2,ykyw2,ykzw2

      ui = (0.d0,1.d0)
      do i = 1, nx
      do j=   1,ny
      do k= 1,nz
         ak2(i,j,k)  =cdexp(-(ykyw2(j)*dt+ykxw2(i)*dt+ykzw2(k)*dt))
     .           -cdexp(ykyw2(j)*dt+ykxw2(i)*dt+ykzw2(k)*dt)
        akC2(i,j,k) =cdexp(-.5d0*(ykyw2(j)*dt+ykxw2(i)*dt+ykzw2(k)*dt))
     .             -cdexp(.5d0*(ykyw2(j)*dt+ykxw2(i)*dt+ykzw2(k)*dt))

      enddo
      enddo
      enddo


      return
      end          
     

      subroutine antial4(h1,f1,g1,h2,f2,g2,h3,f3,g3,h4,f4,g4,nx,ny,nz)

      implicit none
      
      integer*4 nx,ny,nz
      complex*16 h1(nx,ny,nz), f1(nx,ny,nz), g1(nx,ny,nz)
      complex*16 h2(nx,ny,nz),f2(nx,ny,nz),g2(nx,ny,nz)
      complex*16 h3(nx,ny,nz),f3(nx,ny,nz),g3(nx,ny,nz)
      complex*16 h4(nx,ny,nz),f4(nx,ny,nz),g4(nx,ny,nz)
      real*8 kx(1),ky(1),kz(1)
      real*8 dx,dy,dz
      real*8 dx2,dy2,dz2
      real*8 a,ay
      complex*16 c,cy
      complex*16 ui
      integer*4 i,j,k

      common /kw/ kx,ky,kz
      common /step/ dx,dy,dz


      ui = (0.d0,1.d0)
      
      dx2 = .5d0*dx
      dy2 = .5d0*dy
      dz2 = .5d0*dz

      do i = 1, nx   
      do j= 1,ny
      do k = 1,nz
         a = (dx2*kx(i)+dy2*ky(j)+dz2*kz(k))
         c = dcmplx( dcos(a), -dsin(a) )
         h1(i,j,k) = .5d0*( f1(i,j,k) + c*g1(i,j,k) )
         h2(i,j,k) = .5d0*( f2(i,j,k) + c*g2(i,j,k) )
         h3(i,j,k)=.5d0*(f3(i,j,k)+c*g3(i,j,k))
         h4(i,j,k)=.5d0*(f4(i,j,k)+c*g4(i,j,k))

      enddo   
      enddo
      enddo
      return
      end


      subroutine transl4(h1,f1,h2,f2,h3,f3,h4,f4,nx,ny,nz,isgn)


      implicit none
      
      integer*4 nx,ny,nz, isgn
      complex*16 h1(nx,ny,nz), f1(nx,ny,nz)
      complex*16 h2(nx,ny,nz),f2(nx,ny,nz)
      complex*16 h3(nx,ny,nz), f3(nx,ny,nz)
      complex*16 h4(nx,ny,nz), f4(nx,ny,nz)
      real*8 kx(1),ky(1),kx2(1),kz(1),ky2(1),kz2(1)
      real*8 dx,dy,dz
      real*8 dx2,dy2,dz2
      real*8 sg, a,a1
      complex*16 ui
      complex*16 c,c1
      integer*4 i,j,k

      common /kw/ kx,ky,kz /k2/ kx2,ky2,kz2
      common /step/ dx,dy,dz

      
     
     
      
      ui = (0.d0,1.d0)
      dx2 = .5d0*dx
      dy2= .5d0*dy
      dz2 = .5d0*dz
      sg  = dfloat(isgn)
      
      do i = 1, nx
      do j = 1, ny
      do k=1,nz
          
          a =sg*( dx2*kx(i)+dy2*ky(j)+dz2*kz(k))
         c = dcmplx(dcos(a),- dsin(a))
         f1(i,j,k) = c*h1(i,j,k)
         f2(i,j,k)  =c*h2(i,j,k)
         f3(i,j,k)=c*h3(i,j,k)
	 f4(i,j,k)=c*h4(i,j,k)
      enddo
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

      subroutine dcfft3n(ISIGN,DATA3,ngx,ngy,ngz,mx,my,mz,FFT3)
C       real*8 DATA3(ngx,ngy,ngz)
c       real*8 DATAZ(ngz)
      complex*16 DATA3(ngx,ngy,ngz)
      complex*16 DATAZ(ngz)
      complex*16 FFTZ1(ngx,ngy,ngz)
      complex*16 FFT3(ngx,ngy,ngz)
      complex*16 FFTXYZ(ngx)
      complex*16 FFTX(ngx)
      complex*16 FFTY(ngy)
      complex*16 FFTZ(ngz)
      complex*16 FFTZY(ngy)
      complex*16 FFT2(ngx,ngy,ngz)
       do ix=1,ngx
         do iy=1,ngy
            do 11 iz=1,ngz
            DATAZ(iz)=DATA3(ix,iy,iz)
 11    continue
            call dcfftn(ISIGN,DATAZ,FFTZ,mz)
            do 12 iz1=1,ngz
              FFTZ1(ix,iy,iz1)=FFTZ(iz1)
 12    continue
      enddo
      enddo
      do ix=1,ngx
      do iz=1,ngz
      do 13 iy=1,ngy
      FFTY(iy)=FFTZ1(ix,iy,iz)
 13    continue
      call dcfftn(ISIGN,FFTY,FFTZY,my)
      do 14 iy1=1,ngy
      FFT2(ix,iy1,iz)=FFTZY(iy1)
 14    continue
      enddo
      enddo
      do iz=1,ngz
      do iy=1,ngy
      do 15 ix=1,ngx
      FFTX(ix)=FFT2(ix,iy,iz)
 15    continue
      call dcfftn(ISIGN,FFTX,FFTXYZ,mx)
      do 16 ix1=1,ngx
      FFT3(ix1,iy,iz)=FFTXYZ(ix1)
 16    continue
      enddo
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


      
