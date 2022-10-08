c$Id:$
      subroutine elmt06(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2009: Regents of the University of California
c                               All Rights Reserved

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose: Elastic, DKQ plate bending element.

c     Input for material data:
c          MATErial #
c            USER   6
c              ELAStic  E  nu thick q

c-----[--+---------+---------+---------+---------+---------+---------+-]

      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'mdata.h'
      include   'prstrs.h'

      include   'comblk.h'

      real*8        b   ,c   ,aa   ,bb   ,cc   ,dd   ,ee
      common/elmdkq/b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)

      real*8          bm
      common /shpdkq /bm(3,12)

      logical    errck,tinput
      character  txt*15
      integer    i,j,k,l,lint, ndm,ndf,nst,isw, ix(*),is(12)
      real*8     xsj, xx,yy, d1,d2,d3,d6, dn1,dn2,dn3, thick
      real*8     d(*),xl(ndm,*),ul(*),tl(*),s(nst,*),p(ndf,*)
      real*8     sg(3,9),shp(3,8),eps(3),sig(3),td(4),s0(2)

      save

      data       s0 / 2*0.0d0 /

c     Input material properties

      if(isw.eq.1) then
        errck     = tinput(txt,1,td,4)
        d(4)      = td(1)
        d(5)      = td(2)
        thick     = td(3)
        d(6)      = td(4)
        d(1)      = d(4)/(1.0d0-d(5)*d(5))*thick**3/12.0d0
        d(2)      = d(5)*d(1)
        d(3)      = 0.5d0*(d(1)-d(2))
        write(iow,2000) d(4),d(5),thick,d(6)
        d(6)      = d(6)*0.25d0
        ia(1,iel) = 2
        ia(2,iel) = 3

c       Set assembly vector

        do i = 1,4
          is(3*i-2) = ndf*(i-1) + 1
          is(3*i-1) = ndf*(i-1) + 2
          is(3*i  ) = ndf*(i-1) + 3
        end do ! i

c     Compute element tangent array

      elseif(isw.eq.3 .or. isw.eq.6) then

        call jacqud(xl,ndm)
        l = 3
        call int2d(l,lint,sg)
        do l = 1,lint
          call qushp8(sg(1,l),shp,xl,ndm,xsj)
          xsj = xsj*sg(3,l)
          call dktqbm(shp(1,5),shp)

c         Compute weighted jacobian and d-matrix constants

          d1 = d(1)*xsj
          d2 = d(2)*xsj
          d3 = d(3)*xsj
          d6 = d(6)*xsj

c         Compute strains and bending moments

          do i = 1,3
            eps(i) = 0.0d0
            do j = 1,12
              eps(i) = eps(i) + bm(i,j)*ul(is(j))
            end do ! j
          end do ! i
          sig(1)   = d1*eps(1) + d2*eps(2)
          sig(2)   = d2*eps(1) + d1*eps(2)
          sig(3)   = d3*eps(3)
c         sig(3)   = d3*eps(3)*0.5d0

c         Compute element load vector

          p(1,1) = p(1,1) + d6*(1.d0 - sg(1,l))*(1.d0 - sg(2,l))
          p(1,2) = p(1,2) + d6*(1.d0 + sg(1,l))*(1.d0 - sg(2,l))
          p(1,3) = p(1,3) + d6*(1.d0 + sg(1,l))*(1.d0 + sg(2,l))
          p(1,4) = p(1,4) + d6*(1.d0 - sg(1,l))*(1.d0 + sg(2,l))

c         Complete element residual

          k = 0
          do i = 1,4
            do j = 1,3
              k      = k + 1
              p(j,i) = p(j,i) - bm(1,k)*sig(1)
     &                        - bm(2,k)*sig(2)
     &                        - bm(3,k)*sig(3)
            end do ! j
          end do ! i

c         Compute contribution to element stiffness for this point

          do i = 1,12
            dn1 = d1*bm(1,i) + d2*bm(2,i)
            dn2 = d2*bm(1,i) + d1*bm(2,i)
            dn3 = d3*bm(3,i)
            do j = 1,12
              s(is(i),is(j)) = s(is(i),is(j)) + dn1*bm(1,j)
     &                                        + dn2*bm(2,j)
     &                                        + dn3*bm(3,j)
            end do ! j
          end do ! i
        end do ! l

c     Compute and output element variables

      elseif(isw.eq.4) then

        call jacqud(xl,ndm)
        call qushp8(s0,shp,xl,ndm,xsj)
        call dktqbm(shp(1,5),shp)
        do i = 1,3
          eps(i) = 0.0
          do j = 1,12
            eps(i) = eps(i) + bm(i,j)*ul(j)
          end do ! j
        end do ! i
        sig(1) = d(1)*eps(1) + d(2)*eps(2)
        sig(2) = d(2)*eps(1) + d(1)*eps(2)
        sig(3) = d(3)*eps(3)/2.0
        xx = 0.25*(xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))
        yy = 0.25*(xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))
        mct = mct - 1
        if(mct.le.0) then
          mct = 50
          write(iow,2002) o,head
        endif
        write(iow,2003) n,ma,xx,yy,sig

c     Stress contours

      elseif(isw.eq.8) then

        call strsc6(ix,d,xl,ul,shp,hr(nph),hr(nph+numnp),ndm,numnp)

      endif

c     Formats

2000  format(5x,'DKQ Plate Bending Element'/10x,18hmaterial constants//
     & 10x,12hmodulus    = ,e15.5/10x,12hpoisson-nu = ,e15.5/
     & 10x,12hthickness  = ,e15.5/10x,12hq-loading  = ,e15.5/)

2002  format(a1,20a4//5x,'DKQ Plate Moments'//
     & 1x,'elmt',2x,'mat',5x,'1-coord',5x,'2-coord',
     & 5x,'11-moment',5x,'22-moment',5x,'12-moment'/)

2003  format(2i5,2f12.3,3e14.5)

      end

      subroutine jacqud(xl,ndm)

      implicit   none

      real*8        b   ,c   ,aa   ,bb   ,cc   ,dd   ,ee
      common/elmdkq/b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)

      real*8          bm
      common /shpdkq /bm(3,12)

      integer    ndm, i,k
      real*8     xl(ndm,*), sql

      do i = 1,4
        k     =  mod(i,4) + 1
        b(i)  =  xl(2,k) - xl(2,i)
        c(i)  =  xl(1,i) - xl(1,k)
        sql   =  1.0d0/(b(i)*b(i) + c(i)*c(i))
        aa(i) =  1.50d0*c(i)*sql
        bb(i) =  0.75d0*b(i)*c(i)*sql
        cc(i) = (0.25d0*c(i)*c(i) - 0.5d0*b(i)*b(i))*sql
        dd(i) = -1.50d0*b(i)*sql
        ee(i) = (0.25d0*b(i)*b(i) - 0.5d0*c(i)*c(i))*sql
      end do ! i

      end
      subroutine dktqbm(shm,shn)

      implicit   none

      real*8        b   ,c   ,aa   ,bb   ,cc   ,dd   ,ee
      common/elmdkq/b(4),c(4),aa(4),bb(4),cc(4),dd(4),ee(4)

      real*8          bm
      common /shpdkq /bm(3,12)

      integer    i,i1,i2,i3, j
      real*8     shm(3,4),shn(3,4)

      i1 = 1
      do i = 1,4
        j = mod(i+2,4) + 1
        i2 = i1 + 1
        i3 = i2 + 1
        bm(1,i1) = aa(i)*shm(1,i) - aa(j)*shm(1,j)
        bm(1,i2) = bb(i)*shm(1,i) + bb(j)*shm(1,j)
        bm(1,i3) = cc(i)*shm(1,i) + cc(j)*shm(1,j) - shn(1,i)
        bm(2,i1) = dd(i)*shm(2,i) - dd(j)*shm(2,j)
        bm(2,i2) =-ee(i)*shm(2,i) - ee(j)*shm(2,j) + shn(2,i)
        bm(2,i3) =-bb(i)*shm(2,i) - bb(j)*shm(2,j)
        bm(3,i1) = aa(i)*shm(2,i) - aa(j)*shm(2,j)
     &           + dd(i)*shm(1,i) - dd(j)*shm(1,j)
        bm(3,i2) =-ee(i)*shm(1,i) - ee(j)*shm(1,j) + shn(1,i) - bm(2,i3)
        bm(3,i3) = cc(i)*shm(2,i) + cc(j)*shm(2,j) - shn(2,i) - bm(1,i2)
        i1 = i1 + 3
      end do ! i

      end

      subroutine qushp8(s,shp,xl,ndm,xsj)

c     Shape function routine for 8-node serendipity elements

      implicit   none

      integer    ndm,i
      real*8     xsj, xs,xt,ys,yt, ss,tt, sn,tn
      real*8     s(2),shp(3,8),xl(ndm,1),si(4),ti(4)

      data si/-1.d0,1.d0,1.d0,-1.d0/,ti/-1.d0,-1.d0,1.d0,1.d0/

      xs = 0.0d0
      xt = 0.0d0
      ys = 0.0d0
      yt = 0.0d0
      do i = 1,4
        ss = si(i)*s(1)
        tt = ti(i)*s(2)
        sn = si(i)*(1.d0 + tt)
        tn = ti(i)*(1.d0 + ss)
        xs = xs + sn*xl(1,i)
        xt = xt + tn*xl(1,i)
        ys = ys + sn*xl(2,i)
        yt = yt + tn*xl(2,i)
        shp(1,i) = 0.25d0*sn*(ss + ss + tt)
        shp(2,i) = 0.25d0*tn*(ss + tt + tt)
        shp(3,i) = 0.25d0*(1.d0 + ss)*(1.d0 + tt)*(-1.d0 + ss + tt)
      end do ! i
      xsj = (xs*yt - ys*xt)*0.25d0
      xs  = xs/xsj
      xt  = xt/xsj
      ys  = ys/xsj
      yt  = yt/xsj
      xsj = xsj*0.25d0
      do i = 5,7,2
        ss         =  si(i-4)*s(1)
        tt         =  ti(i-4)*s(2)
        shp(1,i)   = -s(1) *(1.d0 + tt)
        shp(2,i)   =  0.5d0*ti(i-4)*(1.d0 - s(1)*s(1))
        shp(3,i)   =  0.5d0*(1.d0 - s(1)*s(1))*(1.d0 + tt)
        shp(1,i+1) = -0.5d0*si(i-4)*(1.d0 - s(2)*s(2))
        shp(2,i+1) = -s(2) *(1.d0 - ss)
        shp(3,i+1) =  0.5d0*(1.d0 - ss)*(1.d0 - s(2)*s(2))
      end do ! i
      do i = 1,8
        sn       = yt*shp(1,i) - ys*shp(2,i)
        shp(2,i) = xs*shp(2,i) - xt*shp(1,i)
        shp(1,i) = sn
      end do ! i

      end

      subroutine strsc6(ix,d,xl,ul,shp,dt,st,ndm,numnp)

      implicit   none

      real*8          bm
      common /shpdkq/ bm(3,12)

      integer    ndm,numnp,ix(*)
      integer    i,ii,j,l
      real*8     d(*),xl(ndm,*),ul(*),shp(3,*),st(numnp,*),dt(numnp)
      real*8     eps(3),sig(3),ss(2,4), xsj

      data ss/-1.d0,-1.d0, 1.d0,-1.d0,  1.d0, 1.d0, -1.d0, 1.d0/

      call jacqud(xl,ndm)

      do l = 1,4
        call qushp8(ss(1,l),shp,xl,ndm,xsj)
        call dktqbm(shp(1,5),shp)
        do i = 1,3
          eps(i) = 0.0d0
          do j = 1,12
            eps(i) = eps(i) + bm(i,j)*ul(j)
          end do ! j
        end do ! i
        sig(1)   = d(1)*eps(1) + d(2)*eps(2)
        sig(2)   = d(2)*eps(1) + d(1)*eps(2)
        sig(3)   = d(3)*eps(3)*0.5d0
        ii       = ix(l)
        dt(ii)   = dt(ii)   + xsj
        st(ii,1) = st(ii,1) + sig(1)*xsj
        st(ii,2) = st(ii,2) + sig(2)*xsj
        st(ii,4) = st(ii,4) + sig(3)*xsj
      end do ! l

      end
