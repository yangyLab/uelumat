       include "umat_mc.for" 
       SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
c
c
c     initialize rhs and amatrx to zero
c
      parameter(zero=0.d0,one=1.0d0)
      parameter (ndim=2, ndof=2, ndi=3, nshr=1, nnodemax=4,
     1     ntens=4, ninpt=4, nsvint=8)
c
      dimension  stiff(ndof*nnodemax,ndof*nnodemax),
     1 force(ndof*nnodemax), shape(nnodemax), dshape(ndim,nnodemax),
     2 xjaci(ndim,ndim), bmat(nnodemax*ndim), statevLocal(nsvint),
     3 stress(ntens), ddsdde(ntens, ntens),
     4 stran(ntens), dstran(ntens), wght(ninpt)
c  
      data wght /one, one, one, one/

      do k1=1, ndof*nnode
        rhs(k1, 1)= zero
        do k2=1, ndof*nnode 
          amatrx(k1, k2)= zero
        end do
      end do
c
c     loop over integration points
c
      do kintk = 1, ninpt
c
c       evaluate shape functions and derivatives
c
        call shapefcn(kintk,ninpt,nnode,ndim,shape,dshape)
c
c       compute coordinates at the integration point
c
c       form B-matrix
c
c        djac = one
        call jacobian(jelem,mcrd,ndim,nnode,coords,dshape,
     1   djac,xjaci,pnewdt)
c
        call bmatrix(xjaci,dshape,nnode,ndim,bmat)
c
c       calculate incremental strains
c 
        call straininc(ntens,ndof,ndim,nnode,mlvarx,bmat,du,dstran)
c        
c       call constitutive routine
c
        call statevar(kintk,nsvint,svars,statevLocal,1)
c
c       material routine for computing state variables
c
c
        call umat(stress,statevLocal,ddsdde,sse,spd,scd,rpl,ddsddt,
     1  drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred,
     2  cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     3  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)
c
        call statevar(kintk,nsvint,svars,statevLocal,0)
c        
c       form stiffness matrix and internal force vector
c
        call stiffmatrix(ntens,nnode,ndim,ndof,
     1       wght(kintk),djac,ddsdde,stress,bmat,stiff,force)
c
c       assemble rhs and lhs
c
        do k1=1, ndof*nnode
            rhs(k1, 1) = rhs(k1, 1) - force(k1) 
          do k2=1, ndof*nnode
            amatrx(k1, k2) = amatrx(k1, k2) + stiff(k1,k2)
          end do
        end do
      end do     
c
      return
      end
c
c shape functions
c
      subroutine shapefcn(kintk,ninpt,nnode,ndim,dN,dNdz)
c
      include 'aba_param.inc'
c
      parameter (dmone=-1.0d0,one=1.0d0,four=4.0d0,eight=8.0d0,
     *     gaussCoord=0.577350269d0)
      parameter (maxElemNode=8,maxDof=3,i2d4node=24,i3d8node=38)
      dimension dN(*),dNdz(ndim,*),coord24(2,4),coord38(3,8)
c     
      data  coord24 /dmone, dmone,
     2                 one, dmone,
     3                 one,  one,
     4               dmone,  one/
c     
C  3D 8-nodes
C
c         determine (g,h)
c
        g = coord24(1,kintk)*gaussCoord
        h = coord24(2,kintk)*gaussCoord
c
c       shape functions
        dN(1) = (one - g)*(one - h)/four;
        dN(2) = (one + g)*(one - h)/four;
        dN(3) = (one + g)*(one + h)/four;
        dN(4) = (one - g)*(one + h)/four;
c
c       derivative d(Ni)/d(g)
        dNdz(1,1) = -(one - h)/four;
        dNdz(1,2) =  (one - h)/four;
        dNdz(1,3) =  (one + h)/four;
        dNdz(1,4) = -(one + h)/four;
c
c       derivative d(Ni)/d(h)
        dNdz(2,1) = -(one - g)/four;
        dNdz(2,2) = -(one + g)/four;
        dNdz(2,3) =  (one + g)/four;
        dNdz(2,4) =  (one - g)/four;
c      
      return
      end
c*****************************************************************
      subroutine jacobian(jelem,mcrd,ndim,nnode,
     1 coords,dshape,djac,xjaci,pnewdt)
c
c     Notation:  ndim ....... element dimension
c                nnode  ..... number of nodes
c                coords ..... coordinates of nodes
c                dshape ..... derivs of shape fcn
c                djac ....... determinant of Jacobian
c                xjaci ...... inverse of Jacobian matrix
c
      include 'aba_param.inc'
      parameter(zero=0.d0, fourth=0.25d0, maxDof=3)
      dimension  xjac(maxDof,maxDof), xjaci(ndim,*), coords(mcrd,*)
      dimension  dshape(ndim,*)
c
      do i = 1, ndim
        do j = 1, ndim
           xjac(i,j)  = zero
           xjaci(i,j) = zero
        end do
      end do
c
      do inod= 1, nnode
          do idim = 1, ndim
            do jdim = 1, ndim
              xjac(jdim,idim) = xjac(jdim,idim) + 
     1        dshape(jdim,inod)*coords(idim,inod)    
            end do
          end do 
      end do
C

        djac = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
        if (djac .gt. zero) then
          ! jacobian is positive - o.k.
          xjaci(1,1) =  xjac(2,2)/djac
          xjaci(2,2) =  xjac(1,1)/djac
          xjaci(1,2) = -xjac(1,2)/djac
          xjaci(2,1) = -xjac(2,1)/djac
        else
          ! negative or zero jacobian
          write(7,*)'WARNING: element',jelem,'has neg. Jacobian'
          pnewdt = fourth
        endif
      return
      end
c*****************************************************************
      subroutine bmatrix(xjaci,dshape,nnode,ndim,bmat)
c
c     Notation:  
c                bmat(i) .....dN1/dx, dN1/dy, dN2/dx, dN2/dy..
c                xjaci ...... inverse Jabobian matrix
c                dshape ......derivative of shape functions
c
      include 'aba_param.inc'
      parameter (zero=0.d0)
      dimension bmat(*),  dshape(ndim,*)
      dimension xjaci(ndim,*)
c
      do i = 1, nnode*ndim
          bmat(i) = zero
      end do
c
      do inod = 1, nnode
        do ider = 1, ndim
          do idim = 1, ndim
             irow = idim + (inod - 1)*ndim
             bmat(irow) = bmat(irow) + 
     1       xjaci(idim,ider)*dshape(ider,inod)      
          end do
        end do
      end do 
c
      return
      end

c*****************************************************************
      subroutine straininc(ntens,ndof,ndim,nnode,
     1 mlvarx,bmat,du,dstran)
c
c     Notation:  
c       dstran(i)  incremental strain component 
c       note:      i = 1   xx direction 
c                    = 2   yy direction 
c                    = 3   zz direction
c                    = 4   xy direction
c    u() - displacement
c   du() - increment of displacement in the last inc.
c   
c
      include 'aba_param.inc'
      parameter(zero=0.d0, one=1.d0)
c
c      dimension dstran(*), bmat(ndim,*), 
c     1     du(mlvarx, *), xdu(3), xx1(3,*), 
c     2      u(mlvarx, *), utmp(3),
c     3      utmpOld(3),xx1Old(3,*),eps(3,3),dInvFold(3,3)
      dimension dstran(*), bmat(ndim,*),du(mlvarx,*),xdu(3)
      do i = 1, ntens
        dstran(i) = zero
      end do
c
c************************************
c    Compute incremental strains
c************************************
c
      do nodi = 1, nnode
c           
         incr_row = (nodi - 1)*ndof
c
         do i = 1, ndof
           xdu(i)= du(i + incr_row,1)
c           utmp(i) = u(i + incr_row,1)
c           utmpOld(i) = utmp(i)-xdu(i)
         end do
c
         dNidx = bmat(1,nodi)
         dNidy = bmat(2,nodi)
c
           dstran(1) = dstran(1) + dNidx*xdu(1)
           dstran(2) = dstran(2) + dNidy*xdu(2)
           dstran(4) = dstran(4) + 
     1          dNidy*xdu(1) + 
     2          dNidx*xdu(2)  
c
      end do
c
      return
      end

c*****************************************************************
      subroutine statevar(npt,nsvint,statev,statev_ip,icopy)
c
c     Transfer data to/from element-level state variable array from/to
c     material-point level state variable array.
c
      include 'aba_param.inc'

      dimension statev(*),statev_ip(*)

      isvinc= (npt-1)*nsvint     ! integration point increment

      if (icopy .eq. 1) then
c
c     Prepare arrays for entry into umat
c
        do i = 1, nsvint
          statev_ip(i)=statev(i+isvinc)
        end do
c
      else
c
c     Update element state variables upon return from umat
c
        do i = 1, nsvint
          statev(i+isvinc)=statev_ip(i)
        end do
      end if
c
      return
      end

c*****************************************************************
      subroutine stiffmatrix(ntens,nnode,ndim,ndof,
     1 weight,djac,ddsdde,stress,bmat,stiff,force)
c
c     Stiffness matrix and internal force contributions at 
c     material integration point
c
      include 'aba_param.inc'

      parameter(zero=0.d0,maxDof=3)

      dimension stiff(ndof*nnode,*), stiff_p(maxDof,maxDof)
      dimension force(*), force_p(maxDof)
      dimension stress(*),bmat(ndim,*),ddsdde(ntens,*)

      dNjdx = zero
      dNjdy = zero
      do i = 1, ndof*nnode
          force(i) = zero
        do j = 1, ndof*nnode
          stiff(j,i) = zero
        end do
      end do
c
      dvol= weight*djac
      do nodj = 1, nnode
c
          incr_col = (nodj - 1)*ndof
c
          dNjdx = bmat(1,nodj)
          dNjdy = bmat(2,nodj)
            force_p(1) = dNjdx*stress(1) + dNjdy*stress(4)
            force_p(2) = dNjdy*stress(2) + dNjdx*stress(4)
c        
          do jdof = 1, ndof
c
            jcol = jdof + incr_col
c
            force(jcol) = force(jcol) +
     &           force_p(jdof)*dvol
c            
          end do
c
          do nodi = 1, nnode
c
            incr_row = (nodi -1)*ndof
c
            dNidx = bmat(1,nodi)
            dNidy = bmat(2,nodi)
c            
              stiff_p(1,1) = dNidx*ddsdde(1,1)*dNjdx
     &             + dNidy*ddsdde(4,4)*dNjdy
     &             + dNidx*ddsdde(1,4)*dNjdy
     &             + dNidy*ddsdde(4,1)*dNjdx
c
              stiff_p(1,2) = dNidx*ddsdde(1,2)*dNjdy
     &             + dNidy*ddsdde(4,4)*dNjdx
     &             + dNidx*ddsdde(1,4)*dNjdx
     &             + dNidy*ddsdde(4,2)*dNjdy
c              
              stiff_p(2,1) = dNidy*ddsdde(2,1)*dNjdx
     &             + dNidx*ddsdde(4,4)*dNjdy
     &             + dNidy*ddsdde(2,4)*dNjdy
     &             + dNidx*ddsdde(4,1)*dNjdx
c
              stiff_p(2,2) = dNidy*ddsdde(2,2)*dNjdy
     &             + dNidx*ddsdde(4,4)*dNjdx
     &             + dNidy*ddsdde(2,4)*dNjdx
     &             + dNidx*ddsdde(4,2)*dNjdy
c              
            do jdof = 1, ndof
              icol = jdof + incr_col
              do idof = 1, ndof
                irow = idof + incr_row
                stiff(irow,icol) = stiff(irow,icol) +
     &          stiff_p(idof,jdof)*dvol
              end do
            end do
          end do
      end do
c
      return
      end
c
