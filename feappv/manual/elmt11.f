c$Id:$
      subroutine elmt11(d,ul,xl,ix,tl,s,r,ndf,ndm,nst,isw)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2009: Regents of the University of California
c                               All Rights Reserved

c-----[--+---------+---------+---------+---------+---------+---------+-]
c     Purpose:  St. Venant torsion problem solution by warping and
c               Prandtl stress function.

c     Input data forms:
c               'elas' = Elastic shear modulus (G)
c               'quad' = Quadrature order (default = 2)
c               'warp' = warping function solution
c                        Assign to: dof = 1 (default if not input)
c               'stre' = stress function solution
c                        Assign to: dof = 1 if no warping
c                                   dof = 2 if warping also done
c     Typical Material Data Input:
c         MATErial MA
c           USER 11
c             ELAStic    G
c             QUADrature NQ
c             WARPing
c             STREss
c
c      N.B. Parts in upper case are required inputs.
c-----[--+---------+---------+---------+---------+---------+---------+-]
      implicit   none

      include   'bdata.h'
      include   'cdata.h'
      include   'eldata.h'
      include   'iofile.h'
      include   'print.h'
      include   'prstrs.h'
      include   'strnum.h'

      include   'comblk.h'

      integer    ndf,ndm,nst,isw

      logical    errck,tinput, pcomp, warp,stre
      character  text*15
      integer    i,j,l, ii,jj, ll,lint, inc
      real*8     jg1,jg2, con, a1,a2, dd
      real*8     gw(3,16), sw(3,16), gs(3,16), ss(3,16)
      real*8     sg(3,16),shp(3,9,16), xsj(16), td(2)
      real*8     xx(16),yy(16)

      integer    ix(*)
      real*8     d(*),xl(ndm,*),ul(ndf,*),tl(*),s(nst,*),r(ndf,*)

c     Input the material properties

      if(isw.eq.1) then
        i     = 2
        text  = 'start'
        do while(.not.pcomp(text,'    ',4))
          errck = tinput(text,1,td,2)

          if(pcomp(text,'elas',4)) then
            d(1) = td(1)
          elseif(pcomp(text,'quad',4)) then
            i    = td(1)
          elseif(pcomp(text,'warp',4)) then
            d(4) = 1.d0
          elseif(pcomp(text,'stre',4)) then
            d(5) = 1.d0
          endif
        end do ! while

        i    = max(2,min(4,i))
        d(2) = i
        write(iow,2000) d(1),i

c       Set warping function solution if nothing specified

        if(max(d(4),d(5)).eq.0.0d0) then
          d(4) = 1.d0
        endif

c       Warping function solution

        if(d(4).gt.0.0d0) then
          write(iow,2001) 
        endif

c       Stress function solution

        if(d(5).gt.0.0d0) then
          write(iow,2002) 
        endif

        istv = 12    ! Maximum number of projected quantities
        lint = 0
        con = 45.0d0/atan(1.0d0)

c     Check element for errors

      elseif(isw.eq.2) then

c     Compute the element tangent array

      elseif(isw.eq.3 .or. isw.eq.6) then

        warp = d(4).gt.0.0d0
        stre = d(5).gt.0.0d0

        if(warp.and.stre) then
          inc = 2
        else
          inc = 1
        endif
        l = nint(d(2))
        call int2d(l,lint,sg)
        do ll = 1,lint
          call shp2d(sg(1,ll),xl,shp(1,1,ll),xsj(ll),
     &               ndm,nel,ix,.false.)
          xsj(ll) = xsj(ll)*sg(3,ll)

c         Compute shearing strains

          xx(ll) = 0.0d0
          yy(ll) = 0.0d0
          do i = 1,nel
            xx(ll) = xx(ll) + shp(3,i,ll)*xl(1,i)
            yy(ll) = yy(ll) + shp(3,i,ll)*xl(2,i)
          end do ! i

c         Warping function solution

          if(warp) then
            gw(1,ll) = 0.0d0
            gw(2,ll) = 0.0d0
            do i = 1,nel
              gw(1,ll) = gw(1,ll) + shp(1,i,ll)*ul(1,i)
              gw(2,ll) = gw(2,ll) + shp(2,i,ll)*ul(1,i)
            end do ! i
            dd =  d(1)
            gw(1,ll) =  gw(1,ll) - yy(ll)
            gw(2,ll) =  gw(2,ll) + xx(ll)
            sw(1,ll) =  dd*gw(1,ll)
            sw(2,ll) =  dd*gw(2,ll)

            ii = 1
            do i = 1,nel
c             Residual

              r(1,i) = r(1,i) - (shp(1,i,ll)*sw(1,ll)
     &                        +  shp(2,i,ll)*sw(2,ll))*xsj(ll)

c             Stiffness

              a1 = shp(1,i,ll)*dd*xsj(ll)
              a2 = shp(2,i,ll)*dd*xsj(ll)

              jj = 1
              do j = 1,nel
                s(ii,jj) = s(ii,jj) + a1*shp(1,j,ll)
     &                              + a2*shp(2,j,ll)
                jj       = jj + ndf
              end do ! j
              ii = ii + ndf
            end do ! i
          endif

c         Stress function solution

          if(stre) then
            ss(1,ll) = 0.0d0
            ss(2,ll) = 0.0d0
            do i = 1,nel
              ss(1,ll) = ss(1,ll) + shp(1,i,ll)*ul(inc,i)
              ss(2,ll) = ss(2,ll) + shp(2,i,ll)*ul(inc,i)
            end do ! i
            dd =  1.d0/d(1)
            gs(1,ll) =  dd*ss(1,ll) + xx(ll)
            gs(2,ll) =  dd*ss(2,ll) + yy(ll)
            ii = inc
            do i = 1,nel
c             Residual

              r(inc,i) = r(inc,i) - (shp(1,i,ll)*gs(1,ll)
     &                            +  shp(2,i,ll)*gs(2,ll))*xsj(ll)

c             Stiffness

              a1 = shp(1,i,ll)*dd*xsj(ll)
              a2 = shp(2,i,ll)*dd*xsj(ll)

              jj = inc
              do j = 1,nel
                s(ii,jj) = s(ii,jj) + a1*shp(1,j,ll)
     &                              + a2*shp(2,j,ll)
                jj       = jj + ndf
              end do ! j
              ii = ii + ndf
            end do ! i
          endif
        end do ! ll

c     Compute and output the element variables

      elseif(isw.eq.4 .or. isw.eq.8) then

        warp = d(4).gt.0.0d0
        stre = d(5).gt.0.0d0

        if(warp.and.stre) then
          inc = 2
        else
          inc = 1
        endif

        if(n.eq.1) then
         jg1 = 0.0d0
         jg2 = 0.0d0
        endif
        l = nint(d(2))
        call int2d(l,lint,sg)
        do ll = 1,lint
          call shp2d(sg(1,ll),xl,shp(1,1,ll),xsj(ll),
     &               ndm,nel,ix,.false.)
          xx(ll) = 0.0d0
          yy(ll) = 0.0d0
          do i = 1,nel
            xx(ll) = xx(ll) + shp(3,i,ll)*xl(1,i)
            yy(ll) = yy(ll) + shp(3,i,ll)*xl(2,i)
          end do ! i

c         Print header

          if(isw.eq.4 .and. prt) then
            mct = mct - 1
            if(mct.le.0) then
              write(iow,4000) o,head
              mct = 17
            endif
          endif

c         Warping function

          if(warp) then
            gw(1,ll) = 0.0d0
            gw(2,ll) = 0.0d0
            do i = 1,nel
              gw(1,ll) = gw(1,ll) + shp(1,i,ll)*ul(1,i)
              gw(2,ll) = gw(2,ll) + shp(2,i,ll)*ul(1,i)
            end do ! i
            gw(1,ll) = gw(1,ll) - yy(ll)
            gw(2,ll) = gw(2,ll) + xx(ll)
            gw(3,ll) = sqrt(gw(1,ll)*gw(1,ll)+gw(2,ll)*gw(2,ll))
            sw(1,ll) = d(1)*gw(1,ll)
            sw(2,ll) = d(1)*gw(2,ll)
            sw(3,ll) = sqrt(sw(1,ll)*sw(1,ll)+sw(2,ll)*sw(2,ll))
            if(isw.eq.4 .and. prt) then
              a1     = con*atan2(gw(2,ll),gw(1,ll))
              a2     = con*atan2(sw(2,ll),sw(1,ll))
              write(iow,4001) n,ma,xx(ll),yy(ll),(sw(ii,ll),ii=1,3),a2,
     &                                 'Warping',(gw(ii,ll),ii=1,3),a1
              if(ior.lt.0) then
                write(*,4001) n,ma,xx(ll),yy(ll),(sw(ii,ll),ii=1,3),a2,
     &                                 'Warping',(gw(ii,ll),ii=1,3),a1
              endif
            endif
          endif

c         Stress function

          if(stre) then
            ss(1,ll) = 0.0d0
            ss(2,ll) = 0.0d0
            do i = 1,nel
              ss(1,ll) = ss(1,ll) + shp(2,i,ll)*ul(inc,i)
              ss(2,ll) = ss(2,ll) - shp(1,i,ll)*ul(inc,i)
            end do ! i
            ss(3,ll) = sqrt(ss(1,ll)*ss(1,ll)+ss(2,ll)*ss(2,ll))
            gs(1,ll) = ss(1,ll)/d(1)
            gs(2,ll) = ss(2,ll)/d(1)
            gs(3,ll) = sqrt(gs(1,ll)*gs(1,ll)+gs(2,ll)*gs(2,ll))
            if(isw.eq.4 .and. prt) then
              a1     = con*atan2(gs(2,ll),gs(1,ll))
              a2     = con*atan2(ss(2,ll),ss(1,ll))
              write(iow,4001) n,ma,xx(ll),yy(ll),(ss(ii,ll),ii=1,3),a2,
     &                                 'Stress ',(gs(ii,ll),ii=1,3),a1
              if(ior.lt.0) then
                write(*,4001) n,ma,xx(ll),yy(ll),(ss(ii,ll),ii=1,3),a2,
     &                                 'Stress ',(gs(ii,ll),ii=1,3),a1
              endif
            endif
          endif
        end do ! ll

        if(warp) then
          do ll = 1,lint
            jg1 = jg1 + (gw(1,ll)*gw(1,ll)+gw(2,ll)*gw(2,ll))
     &                *  xsj(ll)*d(1)*sg(3,ll)
          end do ! ll
        endif
        if(stre) then
          do ll = 1,lint
            jg2 = jg2 - (ss(1,ll)*(gs(1,ll) + 2.d0*yy(ll))
     &                +  ss(2,ll)*(gs(2,ll) - 2.d0*xx(ll)))
     &                *  xsj(ll) * sg(3,ll)
          end do ! ll
        endif

c       Output torsional stiffness

        if(n.eq.numel) then
          if(warp) then
            write(iow,4002) jg1
            if(ior.lt.0) then
              write(*,4002) jg1
            endif
          endif
          if(stre) then
            write(iow,4003) jg2
            if(ior.lt.0) then
              write(*,4003) jg2
            endif
          endif
        endif

c       Project to nodes

        if(isw.eq.8) then
          call sct2d(ix,sw,gw,ss,gs,shp,xsj,
     &               hr(nph),hr(nph+numnp),hr(ner),erav,
     &               lint,nel,numnp,warp,stre)
        endif

c     Compute element mass arrays

      elseif(isw.eq.5) then

      endif

c     Formats

2000  format(5x,'T o r s i o n   E l e m e n t'//
     &      10x,'Material constants'//
     &      10x,'Shear modulus  =',1p,1e15.5/
     &      10x,'Quadrature points  =',i5/)

2001  format(5x,'Warping function solution')

2002  format(5x,'Stress function solution')

4000  format(a1,20a4//5x,'Torsion element stresses and strains'//
     &  ' Elmt Matl   x-Coord   y-Coord    1-Stress    2-Stress ',
     &  ' Max-Stress  Angle'/34x,
     &  '1-Strain    2-Strain  Max-Strain  Angle')

4001  format(2i5,0p,2f10.3,1p,3e12.4,0p,1f8.2/
     &      10x,a,' function',4x,1p,3e12.4,0p,1f8.2/)

4002  format(5x,'Warping function torsion constant (GJ) = ',1p,1e15.7)

4003  format(5x,'Stress  function torsion constant (GJ) = ',1p,1e15.7)

      end

      subroutine sct2d(ix,sw,gw,ss,gs,shp,xsj,dt,st,ser,erav,
     &                 lint,nel,numnp,warp,stre)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Project torsion stresses to nodes

c      Inputs:

c      Outputs:
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      logical   warp, stre
      integer   lint,nel,numnp
      integer   i,l,ll, inc
      real*8    xg, erav

      integer   ix(*)
      real*8    dt(numnp),st(numnp,*),ser(*)
      real*8    xsj(*),shp(3,9,*),sw(3,*),gw(3,*),ss(3,*),gs(3,*)

      save

      if(warp.and.stre) then
        inc = 6
      else
        inc = 0
      endif

c     Lumped projection routine

      do l = 1,lint

c       Compute lumped projection and assemble stress integrals

        do i = 1,nel
          ll = ix(i)
          if(ll.gt.0) then

            xg     = shp(3,i,l)*xsj(l)
            dt(ll) = dt(ll) + xg

c           Stress projections

            if(warp) then
              st(ll,1) = st(ll,1) + sw(1,l)*xg
              st(ll,2) = st(ll,2) + sw(2,l)*xg
              st(ll,3) = st(ll,3) + sw(3,l)*xg
              st(ll,4) = st(ll,4) + gw(1,l)*xg
              st(ll,5) = st(ll,5) + gw(2,l)*xg
              st(ll,6) = st(ll,6) + gw(3,l)*xg
            endif

            if(stre) then
              st(ll,1+inc) = st(ll,1+inc) + ss(1,l)*xg
              st(ll,2+inc) = st(ll,2+inc) + ss(2,l)*xg
              st(ll,3+inc) = st(ll,3+inc) + ss(3,l)*xg
              st(ll,4+inc) = st(ll,4+inc) + gs(1,l)*xg
              st(ll,5+inc) = st(ll,5+inc) + gs(2,l)*xg
              st(ll,6+inc) = st(ll,6+inc) + gs(3,l)*xg
            endif

c           Error estimation projection

            ser(ll)  = ser(ll)  + erav*xg

          endif
        end do
      end do

      end
