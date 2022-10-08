c$Id:$
      subroutine elmt07(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2009: Regents of the University of California
c                               All Rights Reserved

c-----[--+---------+---------+---------+---------+---------+---------+-]
c      Purpose: Elastic, DKT plate bending element
c
c     D.O.F. order
c         1 = w
c         2 = theta-x
c         3 = theta-y
c
c     Input parameters
c         ndm = 2
c         ndf = 3
c         nen = 3 (or more)
c
c     Material parameters
c         E   - modulus of elasticity
c         nu  - Poisson ratio
c         h   - thickness
c         q   - uniform loading

c     Input for material data:
c          MATErial #
c            USER   6
c              ELAStic  E  nu h q

c---------------------------------------------------------------
c     Ref: Batoz, Bathe, Ho, "A Study of 3-Node Triangular Plate
c          Bending Elements," IJNME, v 15, no. 12, pp 1771-1812,
c          1980
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'prstrs.h'
      include   'comblk.h'

      real*8        b   ,c   ,bd   ,cd   ,aa   ,bb   ,cc   ,dd   ,ee
      common /elm7/ b(3),c(3),bd(3),cd(3),aa(3),bb(3),cc(3),dd(3),ee(3)

      real*8        xsj,bm
      common /shp7/ xsj,bm(3,9)

      logical    errck, tinput
      character  txt*15
      integer    ndm,ndf,nst,isw, ix(*), is(9)
      integer    i,j, l
      real*8     d1,d2,d3, dn1,dn2,dn3, xx,yy, thick
      real*8     d(*),xl(ndm,*),ul(*),tl(*),s(nst,*),p(nst)
      real*8     el(3,3),eps(3),sig(3)

      save 

c     Area coordinates stored in 'el'

      data el/0.5d0,0.0d0,0.5d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0,0.5d0/

c     Input material properties

      if(isw.eq.1) then

        errck = tinput(txt,1,d(4),4)
        thick = d(6)
        d(6)  = d(7)
        d(1)  = d(4)/(1.-d(5)*d(5))*thick**3/12.0
        d(2)  = d(5)*d(1)
        d(3)  = 0.5*(d(1)-d(2))
        write(iow,2000) d(4),d(5),thick,d(6)
        if(ior.lt.0) then
          write(*,2000) d(4),d(5),thick,d(6)
        endif

        do i = 1,3
          is(3*i-2) = ndf*(i-1) + 1
          is(3*i-1) = ndf*(i-1) + 2
          is(3*i  ) = ndf*(i-1) + 3
        end do ! i

c     Compute element tangent array

      elseif(isw.eq.3 .or. isw.eq.6) then

        call jactri(xl,ndm)
        xsj = xsj/6.0d0

c       Compute weighted jacobian and d-matrix constants

        d1 = d(1)*xsj
        d2 = d(2)*xsj
        d3 = d(3)*xsj

c       Compute element load vector

        p(1) = d(6)*xsj
        p(4) = p(1)
        p(7) = p(1)
        do l = 1,3
          call dktbm(el(1,l))

c         Compute strains

          do i = 1,3
            eps(i) = 0.0
            do j = 1,9
              eps(i) = eps(i) + bm(i,j)*ul(is(j))
            end do ! j
          end do ! i
          dn1 = d1*eps(1) + d2*eps(2)
          dn2 = d2*eps(1) + d1*eps(2)
          dn3 = d3*eps(3)
          do j = 1,9
            p(is(j)) = p(is(j)) - bm(1,j)*dn1
     &                          - bm(2,j)*dn2
     &                          - bm(3,j)*dn3
          end do ! j

c         Compute contribution to element stiffness for this point

          do i = 1,9
            dn1 = d1*bm(1,i) + d2*bm(2,i)
            dn2 = d2*bm(1,i) + d1*bm(2,i)
            dn3 = d3*bm(3,i)
            do j = 1,9
              s(is(i),is(j)) = s(is(i),is(j)) + dn1*bm(1,j)
     &                                        + dn2*bm(2,j)
     &                                        + dn3*bm(3,j)
            end do ! j
          end do ! i
        end do ! l

c     Compute and output element variables

      elseif(isw.eq.4) then

        call jactri(xl,ndm)
        eps(1) = 1.d0/3.d0
        eps(2) = eps(1)
        eps(3) = eps(1)
        call dktbm(eps)
        do i = 1,3
          eps(i) = 0.0
          do j = 1,9
            eps(i) = eps(i) + bm(i,j)*ul(j)
          end do ! j
        end do ! i
        sig(1) = d(1)*eps(1) + d(2)*eps(2)
        sig(2) = d(2)*eps(1) + d(1)*eps(2)
        sig(3) = d(3)*eps(3)*0.5d0
        xx = (xl(1,1)+xl(1,2)+xl(1,3))/3.d0
        yy = (xl(2,1)+xl(2,2)+xl(2,3))/3.d0
        mct = mct - 1
        if(mct.le.0) then
          mct = 50
          write(iow,2002) o,head
          if(ior.lt.0) then
            write(*,2002) o,head
          endif
        endif
        write(iow,2003) n,ma,xx,yy,sig
        if(ior.lt.0) then
          write(*,2003) n,ma,xx,yy,sig
        endif

c     Stress contours

      elseif(isw.eq.8) then

        call strsc7(ix,is,d,xl,ul,hr(nph),hr(nph+numnp),ndm,numnp)

      endif ! isw

c     Formats

2000  format(5x,'DKT Triangular Plate Element'/9x,'Material Constants'
     & //10x,'Modulus    =',e15.5/10x,'Poisson-nu =',e15.5/
     &   10x,'Thickness  =',e15.5/10x,'q-loading  =',e15.5/)

2002  format(a1,20a4//5x,'DKT Triangular Plate Moments'//
     & 1x,'elmt  mat     1-coord     2-coord',
     & 5x,'11-moment     22-moment     12-moment'/)

2003  format(2i5,2f12.3,3e14.5)

      end

      subroutine jactri(xl,ndm)

      implicit   none

      real*8        b   ,c   ,bd   ,cd   ,aa   ,bb   ,cc   ,dd   ,ee
      common /elm7/ b(3),c(3),bd(3),cd(3),aa(3),bb(3),cc(3),dd(3),ee(3)

      real*8        xsj,bm
      common /shp7/ xsj,bm(3,9)

      integer    ndm, i,j,k
      real*8     xl(ndm,*), cs,bs,sql

c     Form terms for jacobian of a triangle

      do i = 1,3
        j = mod(i,3) + 1
        k = mod(j,3) + 1
        b(i) = xl(2,k) - xl(2,j)
        c(i) = xl(1,j) - xl(1,k)
        sql = (b(i)+b(i))*(b(i)+b(i)) + (c(i)+c(i))*(c(i)+c(i))
        cs = c(i)/sql
        bs = b(i)/sql
        aa(i) = 6.*cs
        bb(i) = 3.*bs*c(i)
        cc(i) = c(i)*cs - b(i)*(bs+bs)
        dd(i) =-6.*bs
        ee(i) = b(i)*bs - c(i)*(cs+cs)
      end do ! i
      xsj = xl(1,1)*b(1) + xl(1,2)*b(2) + xl(1,3)*b(3)
      xsj = - xsj
      do i = 1,3
        bd(i) = b(i)/xsj
        cd(i) = c(i)/xsj
      end do ! i

      end

      subroutine dktbm(el)

      implicit   none

      real*8        b   ,c   ,bd   ,cd   ,aa   ,bb   ,cc   ,dd   ,ee
      common /elm7/ b(3),c(3),bd(3),cd(3),aa(3),bb(3),cc(3),dd(3),ee(3)

      real*8        xsj,bm
      common /shp7/ xsj,bm(3,9)

      integer    i,j,k, i1
      real*8     shm1,shm2
      real*8     a1(3),a2(3),b1(3),b2(3),c1(3),c2(3),d1(3),d2(3)
      real*8     e1(3),e2(3),el(3),shn(2,3)

c     Form shape function terms for strains

      do i = 1,3
        j = mod(i,3) + 1
        k = mod(j,3) + 1
        shn(1,i) = bd(i)*(el(i)*4.d0 - 1.0)
        shn(2,i) = cd(i)*(el(i)*4.d0 - 1.0)
        shm1     =(bd(j)*el(k) + bd(k)*el(j))*4.d0
        shm2     =(cd(j)*el(k) + cd(k)*el(j))*4.d0
        a1(i)    = aa(i)*shm1
        a2(i)    = aa(i)*shm2
        b1(i)    = bb(i)*shm1
        b2(i)    = bb(i)*shm2
        c1(i)    = cc(i)*shm1
        c2(i)    = cc(i)*shm2
        d1(i)    = dd(i)*shm1
        d2(i)    = dd(i)*shm2
        e1(i)    = ee(i)*shm1
        e2(i)    = ee(i)*shm2
      end do ! i

c     Form strain-displacement matrix

      i1 = 1
      do i = 1,3
        j = mod(i,3) + 1
        k = mod(j,3) + 1
        bm(1,i1  ) = a1(k) - a1(j)
        bm(1,i1+1) = b1(k) + b1(j)
        bm(1,i1+2) = c1(k) + c1(j) - shn(1,i)
        bm(2,i1  ) = d2(k) - d2(j)
        bm(2,i1+1) =-e2(k) - e2(j) + shn(2,i)
        bm(2,i1+2) =-b2(k) - b2(j)
        bm(3,i1  ) = a2(k) - a2(j) + d1(k) - d1(j)
        bm(3,i1+1) =-e1(k) - e1(j) + shn(1,i) - bm(2,i1+2)
        bm(3,i1+2) = c2(k) + c2(j) - shn(2,i) - bm(1,i1+1)
        i1 = i1 + 3
      end do ! i

      end

      subroutine strsc7(ix,is,d,xl,ul,dt,st,ndm,numnp)

      implicit   none

      real*8        xsj,bm
      common /shp7/ xsj,bm(3,9)

      integer    ndm,numnp,ix(*),is(*)
      integer    i,ii,j,l
      real*8     d(*),xl(ndm,*),ul(*),st(numnp,*),dt(numnp)
      real*8     eps(3),sig(3),el(3,3)

c     Area coordinates stored in 'el'

      data el/0.5d0,0.0d0,0.5d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0,0.5d0/

c     Compute center stress

      call jactri(xl,ndm)

        do l = 1,3
          call dktbm(el(1,l))

c         Compute strains

          do i = 1,3
            eps(i) = 0.0
            do j = 1,9
              eps(i) = eps(i) + bm(i,j)*ul(is(j))
            end do ! j
          end do ! i
          sig(1) = d(1)*eps(1) + d(2)*eps(2)
          sig(2) = d(2)*eps(1) + d(1)*eps(2)
          sig(3) = d(3)*eps(3)
          do j = 1,3
            ii       = ix(j)
            dt(ii)   = dt(ii)   + xsj
            st(ii,1) = st(ii,1) + sig(1)*xsj
            st(ii,2) = st(ii,2) + sig(2)*xsj
            st(ii,4) = st(ii,4) + sig(3)*xsj
          end do ! j
        end do ! l

      end
