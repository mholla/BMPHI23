************************************************************************
!
! User element for transient or steady state heat transfer
!  in 2D or 3D.  This is for planar, axisymetric, and 3D.
!
! Solution variables (or nodal variables) are the temperatures.
!
! This subroutine is for the following element types
!  > two-dimensional 4 node isoparametric element as shown below
!       with 1pt (reduced) or 4pt (full) gauss integration.
!  > three-dimensional 8 node isoparametric element as shown below
!       with 1pt (reduced) or 8pt (full) gauss integration.
!
! Surface flux boundary conditions are supported in the following
!  elements.  Based on our convention, the face on which the fliud
!  flux is applied is the "label", i.e.
!  - U1,U2,U3,U4,... refer to fluid fluxes applied to faces
!                     1,2,3,4,... respectively,
!
!
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
!
! Shuolun Wang, Feb, 2022
!
! advection term is added in suppressor
!
! use with
! 2dsquareAdv_sup.inp
! 2dsquareAdv_sup_unit.inp
! Larva2.inp
! 2dots.inp
! 3dots.inp
! 4dots.inp
! 5dots.inp
! 3dsquareAdv_sup.inp
!
! to verify the code:
! verification_2d.inp
! verification_3d.inp
***********************************************************************
!
! User element statement in the input file (set ? values as needed):
!
!  2D elements
!  *User Element,Nodes=4,Type=U?,Iproperties=2,Properties=9,Coordinates=2,Variables=?,Unsymm
!  1,2,11
!
!  3D elements
!  *User Element,Nodes=8,Type=U3,Iproperties=2,Properties=9,Coordinates=3,Variables=?,Unsymm
!  1,2,3,11
!
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!       1) polymer volume fraction (phi)
!
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nIntPt
!          svars(1+j) = phi ---- polymer volume fraction at integ pt k
!          j = j + nlSdv
!       end loop over k
!
!     In the input file, set 'User output variables'= number of global SDV's
!
!     In the input file, set 'ngSdv'= number of global SDV's
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'varibles'=(nlSdv*nIntPt)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Ktherm = props(1)  ! Thermal conductivity
!     Cheat  = props(2)  ! Specific heat
!     rho    = props(3)  ! Density
!     theta0 = props(4)  ! Initial temperature
!     Tconv  = props(5)  ! Ambient temperature for convection BC
!     Trad   = props(6)  ! Ambient temperature for radiation BC
!     nlSdv  = jprops(1) ! Number of local sdv's per integ pt
!     ngSdv  = jprops(2) ! Number of global sdv's per integ pt
!
!***********************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and
      !    that value must be set here the same.

      integer numElem,ElemOffset,err,nInt,nIntS

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=10000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=10000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of volume integration points
      parameter(nInt=4)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of surface integration points
      parameter(nIntS=1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c      real*8, allocatable :: globalSdv(:,:,:)
      real*8 globalSdv(numElem,nInt,4)

      end module global

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.

      use global

      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

      uvar(1) = globalSdv(noel-ElemOffset,npt,1) !
      uvar(2) = globalSdv(noel-ElemOffset,npt,2) !
      uvar(3) = globalSdv(noel-ElemOffset,npt,3) !
      uvar(4) = globalSdv(noel-ElemOffset,npt,4) !
c      do i=1,nUvarm
c         uvar(i) = globalSdv(noel-ElemOffset,npt,i)
c      enddo

      return
      end subroutine uvarm

****************************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      integer lenJobName,lenOutDir,nDim
      real*8 transient
      character*256 jobName,outDir,fileName


      !----------------------------------------------------------------
      !
      ! Perform initial checks
      !
      !
      ! Open the debug/error message file
      !
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      fileName = outDir(1:lenOutDir)//'\aaMSGS_'//
     +     jobName(1:lenJobName)//'.dat'
      open(unit=80,file=fileName,status='unknown')


      ! Check the procedure type, this should be a coupled
      !  temperature displacement or pore pressure displacement
      !  which are any of the following (62-65,71-73)
      !
      if((lflags(1).eq.64).or.(lflags(1).eq.65).or.
     +     (lflags(1).eq.62).or.(lflags(1).eq.63).or.
     +     (lflags(1).eq.71).or.
     +     (lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! all is good
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         write(80,*) 'Abaqus does not have the right procedure'
         write(80,*) 'go back and check the procedure type'
         write(80,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem, I think this only matters
      !  when it comes to output in viewer
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a small displacement analysis'
         write(80,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         write(80,*) 'Abaqus thinks you are doing'
         write(80,*) 'a linear purturbation step'
         call xit
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.0.0) return
      !
      ! Done with initial checks
      !
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !
      ! Check for steady state or transient analysis
      !
      if((lflags(1).eq.62).or.(lflags(1).eq.63).or.
     +     (lflags(1).eq.71)) then
         !
         ! This is steady state
         !
         transient = 0.d0
      else
         !
         ! This is transient
         !
         transient = 1.d0
      endif
      !
      ! Done with check for steady state or transient
      !
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !
      ! Call the paricular element to perform the analysis
      !
      if(jtype.eq.1) then
         !
         ! This is a plane strain analysis
         !
         nDim = 2
         call UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,transient)
         !
         !
      elseif(jtype.eq.2) then
         !
         ! This is an axisymmetric analysis
         !
         nDim = 2
         call UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,transient)
         !
         !
      elseif(jtype.eq.3) then
         !
         ! This is a 3D analysis
         !
         nDim = 3
         call U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +        PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +        TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +        PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +        NJPROP,PERIOD,
     +        nDim,transient)
         !
         !
      else
         !
         ! We have a problem...
         !
         write(*,*) 'Element type not supported, jtype=',jtype
         write(80,*) 'Element type not supported, jtype=',jtype
         call xit
         !
      endif
      !
      ! Done with this element, RHS and AMATRX already returned
      !  as output from the specific element routine called
      !
      !----------------------------------------------------------------


      return
      end subroutine uel

************************************************************************

      subroutine UPE4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,transient)

      use global
*
      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 thetaNew(nNode),thetaOld(nNode),dtheta(nNode)
      real*8 u(nNode,nDim),du(nNode,ndofel),uOld(nNode,nDim)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer ii,jj,a12,b12,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nDofN

      real*8 Iden(3,3),Le,theta0,Rt(nNode,1),Ktt(nNode,nNode),sh(nNode)
      real*8 w(nInt),ds,flux,detMapJ,dsh(nNode,nDim),theta_tau,theta_t
      real*8 dTdX(nDim,1),dTdt,xi(nInt,2),Ktherm,Cheat,dKdT,dCdT,rho
      real*8 wS(nIntS),xLocal(nIntS),yLocal(nIntS),Nvec(1,nNode),ResFac
      real*8 transient,heatGen,dHdT,dshxi(nNode,nDim),TanFac
      real*8 Ru(nDim*nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 Rv(nNode,1),Kvv(nNode,nNode),vNew(nNode),dv(nNode),vOld(nNode)
      real*8 krate,Diffv,ResFac2,dDvdv
      real*8 TanFac2,coordsC(mcrd,nNode)
      real*8 Rmu(nNode,1),Rphi(nNode,1),Kmumu(nNode,nNode)
      real*8 Kphiphi(nNode,nNode),Kmuphi(nNode,nNode)
      real*8 Kphimu(nNode,nNode)
      real*8 muNew(nNode),dmu(nNode),muOld(nNode)
      real*8 phiNew(nNode),dphi(nNode),phiOld(nNode)
      real*8 mu_tau,mu_t,dmudt,dmudX(nDim,1)
      real*8 phi_tau,phi_t,dphidt,dphidX(nDim,1)
      real*8 mobility,eps,dedphi,dmdmu
      real*8 Rtemp(nNode,1),Ktemptemp(nNode,nNode),Kphitemp(nNode,nNode)
      real*8 Ktempphi(nNode,nNode)
      real*8 tempNew(nNode),dtemp(nNode),tempOld(nNode)
      real*8 temp_tau,temp_t,dtempdt,dtempdX(nDim,1)
      real*8 Stefan,C2,C1
      real*8 R_u(nNode,1),R_v(nNode,1)
      real*8 K_uu(nNode,nNode),K_vv(nNode,nNode)
      real*8 K_uv(nNode,nNode),K_vu(nNode,nNode)
      real*8 actNew(nNode),dact(nNode),actOld(nNode)
      real*8 supNew(nNode),dsup(nNode),supOld(nNode)
      real*8 u_tau,u_t,dudt,dudX(nDim,1)
      real*8 v_tau,v_t,dvdt,dvdX(nDim,1)
      real*8 gamma,pa,pb,diffusivity,pe
      real*8 c_t,c_tau
      real*8 Diff_u,Diff_v,F_const,k_const
      real*8 alpha,beta
      real*8 R_phi(nNode,1),K_phiphi(nNode,nNode),K_uphi(nNode,nNode)
      real*8 K_vphi(nNode,nNode),K_phiu(nNode,nNode),K_phiv(nNode,nNode)
      real*8 permittivity,sigma


      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)






      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero   !  residual of displacement
      R_u = zero   !  residual of activator
      R_v = zero   !  residual of suppressor
      R_phi = zero !  residual of electric potenential



      ! diagonal
      Kuu  = zero     ! tangent of displacement
      K_uu = zero     ! tangent of activator
      K_vv = zero     ! tangent of suppressor
      K_phiphi = zero ! tangent of electric potential


      ! off-diagonal
      K_uv     = zero   ! tangent of activator - suppressor
      K_uphi   = zero   ! tangent of activator - potential
      K_vu     = zero   ! tangent of suppressor - activator
      K_vphi   = zero   ! tangent of suppressor - potential
      K_phiu   = zero   ! tangent of potential - activator
      K_phiv   = zero   ! tangent of potential - suppressor


      Energy = zero









      ! Obtain nodal degree of freedom
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         ! the activator
         k = k + 1
         actNew(i) = Uall(k)
         dact(i) = DUall(k,1)
         actOld(i) = actNew(i) - dact(i)
         ! the suppressor
         k = k + 1
         supNew(i) = Uall(k)
         dsup(i) = DUall(k,1)
         supOld(i) = supNew(i) - dsup(i)
         ! the electric potential
         k = k + 1
         phiNew(i) = Uall(k)
         dphi(i) = DUall(k,1)
         phiOld(i) = phiNew(i) - dphi(i)
      enddo

      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo




      ! Impose any time-stepping changes on the increments of
      !
      !
c      do i=1,nNode
c         if(dabs(dact(i)).gt.1.d6) then
c            pnewdt = 0.5
c            return
c         endif
c      enddo
      !
      !
      !
c      do i=1,nNode
c         if(dabs(dsup(i)).gt.1.d6) then
c            pnewdt = 0.5
c            return
c         endif
c      enddo






      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions
            !
            c_t  = 0.00001
            !
         else
            !
            ! this is not the first increment, read old values
            !
            c_t  = svars(1+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif



         ! Obtain activator and its derivative's at
         !  this intPt at the begining and end of the incrment
         !
         u_tau = zero
         u_t = zero
         dudt = zero
         dudX = zero
         do k=1,nNode
            u_tau = u_tau + actNew(k)*sh(k)
            u_t   = u_t + actOld(k)*sh(k)
            do i=1,nDim
               dudX(i,1) = dudX(i,1) + actNew(k)*dsh(k,i)
            enddo
         enddo
         dudt = (u_tau - u_t)/dtime



         ! Obtain the species v and its derivative's at
         !  this intPt at the begining and end of the incrment
         !
         v_tau = zero
         v_t = zero
         dvdt = zero
         dvdX = zero
         do k=1,nNode
            v_tau = v_tau + supNew(k)*sh(k)
            v_t   = v_t + supOld(k)*sh(k)
            do i=1,nDim
               dvdX(i,1) = dvdX(i,1) + supNew(k)*dsh(k,i)
            enddo
         enddo
         dvdt = (v_tau - v_t)/dtime



         ! Obtain the electric potential and its derivative's at
         !  this intPt at the begining and end of the incrment
         !
         phi_tau = zero
         phi_t = zero
         dphidt = zero
         dphidX = zero
         do k=1,nNode
            phi_tau = phi_tau + phiNew(k)*sh(k)
            phi_t   = phi_t + phiOld(k)*sh(k)
            do i=1,nDim
               dphidX(i,1) = dphidX(i,1) + phiNew(k)*dsh(k,i)
            enddo
         enddo
         dphidt = (phi_tau - phi_t)/dtime



         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! I Stopped here, I should update the integ later!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point
         !
         call integ(props,nprops,dtime,
     +     u_tau,u_t,v_tau,v_t,
     +     Diff_u,Diff_v,F_const,k_const,sigma,permittivity,
     +     c_t,c_tau)
         !
         ! I need sigma, permittivity
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



c$$$         ! Time stepping algorithim based on the constitutive response.
c$$$         !
c$$$         phiLmt = 0.005d0
c$$$         umeror = dabs((phi_tau - phi_t)/phiLmt)
c$$$         if(umeror.le.half) then
c$$$            pnewdt = 1.5d0
c$$$         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
c$$$            pnewdt = 1.25d0
c$$$         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
c$$$            pnewdt = 0.75d0
c$$$         else
c$$$            pnewdt = half
c$$$         endif

         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = c_tau
         jj = jj + nlSdv ! setup for the next intPt



         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = c_tau
         globalSdv(jelem,intPt,2) = - dphidX(1,1) ! x component of electric field
         globalSdv(jelem,intPt,3) = - dphidX(2,1) ! y component of electric field



         ! calculate some lengthy quantities
         alpha = - u_tau*v_tau**2.0 + F_const*(1.0 - u_tau)
         beta  = u_tau*v_tau**2.0 - (F_const+k_const)*v_tau

         ! Compute/update the activator residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         ! I need gamma, pa, pb
         ! I need sigma
c         R_u = R_u +  detmapJ*w(intpt)*
c     +         (
c     +         transpose(Nvec)*dudt
c     +       + matmul(dsh,dudX)*Diff_u
c     +       - transpose(Nvec)*alpha
c     +       - sigma*matmul(transpose(Nvec),matmul(transpose(dphidX),dudX))
c     +         )
         R_u = R_u +  detmapJ*w(intpt)*
     +         (
     +         transpose(Nvec)*dudt
     +       + matmul(dsh,dudX)*Diff_u
     +       - transpose(Nvec)*alpha
     +         )




         ! Compute/update the suppressor residual vector
         !
         ! I need diffusivity,pe
c         R_v = R_v + detmapJ*w(intpt)*
c     +         (
c     +         transpose(Nvec)*dvdt
c     +       + Diff_v*matmul(dsh,dvdX)
c     +       - transpose(Nvec)*beta
c     +         )
         R_v = R_v + detmapJ*w(intpt)*
     +         (
     +         transpose(Nvec)*dvdt
     +       + Diff_v*matmul(dsh,dvdX)
     +       - transpose(Nvec)*beta
     +       - sigma*matmul(transpose(Nvec),matmul(transpose(dphidX),dvdX))
     +         )



         ! Compute/update the eletric residual vector
         !
         ! I need permittivity
c         R_phi = R_phi + detmapJ*w(intpt)*
c     +         (
c     +         - permittivity*matmul(dsh,dphidX)
c     +         )
         R_phi = R_phi + detmapJ*w(intpt)*
     +         (
     +         matmul(dsh,dphidX)
     +         )



         ! Compute/update the K_uu
         !
         !
c         K_uu = K_uu + detmapJ*w(intPt)*
c     +          (
c     +          -matmul(transpose(Nvec),Nvec)*(1.0/dtime)
c     +          -matmul(dsh,transpose(dsh))*Diff_u
c     +          -matmul(transpose(Nvec),Nvec)*(v_tau**2.0 + F_const)
c     +          +sigma*matmul(matmul(transpose(Nvec),transpose(dphidX)),transpose(dsh))
c     +          )
         K_uu = K_uu + detmapJ*w(intPt)*
     +          (
     +          -matmul(transpose(Nvec),Nvec)*(1.0/dtime)
     +          -matmul(dsh,transpose(dsh))*Diff_u
     +          -matmul(transpose(Nvec),Nvec)*(v_tau**2.0 + F_const)
     +          )



         ! Compute/update the K_vv
         !
         !
c         K_vv = K_vv + detmapJ*w(intPt)*
c     +            (
c     +            -(1.0/dtime)*matmul(transpose(Nvec),Nvec)
c     +            -Diff_v*matmul(dsh,transpose(dsh))
c     +            +matmul(transpose(Nvec),Nvec)*(2.0*u_tau*v_tau-F_const-k_const)
c     +            )
         K_vv = K_vv + detmapJ*w(intPt)*
     +            (
     +            -(1.0/dtime)*matmul(transpose(Nvec),Nvec)
     +            -Diff_v*matmul(dsh,transpose(dsh))
     +            +matmul(transpose(Nvec),Nvec)*(2.0*u_tau*v_tau-F_const-k_const)
     +            +sigma*matmul(matmul(transpose(Nvec),transpose(dphidX)),transpose(dsh))
     +            )




         ! Compute/update the K_phiphi
         !
         !
c         K_phiphi = K_phiphi + detmapJ*w(intPt)*
c     +            (
c     +            +permittivity*matmul(dsh,transpose(dsh))
c     +            )
         K_phiphi = K_phiphi + detmapJ*w(intPt)*
     +            (
     +            -matmul(dsh,transpose(dsh))
     +            )




         ! Compute/update the K_uv
         !
         !
         K_uv  = K_uv  + detmapJ*w(intPt)*
     +           (
     +           -matmul(transpose(Nvec),Nvec)*2.0*u_tau*v_tau
     +           )


         ! Compute/update the K_uphi
         !
         !
c         K_uphi  = K_uphi  + detmapJ*w(intPt)*
c     +           (
c     +           +sigma*matmul(matmul(transpose(Nvec),transpose(dudX)),transpose(dsh))
c     +           )
         K_uphi  = zero



         ! Compute/update the K_vu
         !
         !
         K_vu = K_vu + detmapJ*w(intPt)*
     +           (
     +           +(v_tau**2.0)*matmul(transpose(Nvec),Nvec)
     +           )


         ! Compute/update the K_vphi
         !
         !
c         K_vphi = zero
         K_vphi  = K_vphi  + detmapJ*w(intPt)*
     +           (
     +           +sigma*matmul(matmul(transpose(Nvec),transpose(dvdX)),transpose(dsh))
     +           )



         K_phiu = zero
         K_phiv = zero



         ! The displacement residual and tangent
         !
c         do i=1,nNode
c            Ru(i,1) = zero
c            Kuu(i,i) = one
c         enddo
         Ru = zero
         Kuu = one


      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Start loop over surface fluid flux terms here
      !
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the flux
            !  acts on is the flux ``label''
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! magnitude

            if((face.ge.1).and.(face.le.4)) then
               !
               ! heat flux applied
               !
               if(face.eq.1) then
                  faceFlag = 1
               elseif(face.eq.2) then
                  faceFlag = 2
               elseif(face.eq.3) then
                  faceFlag = 3
               elseif(face.eq.4) then
                  faceFlag = 4
               endif
               !
               if(nIntS.eq.1) then
                  call xintSurf2D1pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.2) then
                  call xintSurf2D2pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.3) then
                  call xintSurf2D3pt(faceFlag,xLocal,yLocal,wS)
               else
                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
                  write(80,*) 'Invalid nIntS points, nIntS=',nIntS
                  call xit
               endif
               !
               ! loop over integ points
               !
               do ii=1,nIntS
                  !
                  ! Compute shape functions, derivatives, and the
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(nNode,xLocal(ii),yLocal(ii),faceFlag,
     +                 coords,sh,ds)
                  !
                  ! Modify the residual, loop over nodes, recall
                  !  sh(n)=0 when n is not on this face
                  !
                  do n=1,nNode
                     Rt(n,1) = Rt(n,1) - wS(ii)*ds*sh(n)*flux
                  enddo
                  !
                  ! No change to the tangent matrix since we do
                  !  not have temperature dependence here
                  !
               enddo ! loop over nIntS
               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over flux and traction terms
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the tangent matrix
      !

      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode

!      write(*,*) 'nDofN=',nDofN

      ! The residual vector (RHS)
      !
      do i=1,nNode
         A11 = nDofN*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! displacement
         !
         rhs(A11,1) = Ru(A12,1)
         rhs(A11+1,1) = Ru(A12+1,1)
         !
         ! activator
         !
         rhs(A11+2,1) = R_u(i,1)
         !
         ! suppressor
         !
         rhs(A11+3,1) = R_v(i,1)
         !
         ! potential
         !
         rhs(A11+4,1) = R_phi(i,1)
         !
         !
      enddo


      ! Assemble the element tangent matrix (AMATRX)
      !
      do i=1,nNode
         do j=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            B11 = nDofN*(j-1)+1
            B12 = nDim*(j-1)+1
            !
            ! displacement
            !
            amatrx(A11,B11) = Kuu(A12,B12)
            amatrx(A11,B11+1) = Kuu(A12,B12+1)
            amatrx(A11+1,B11) = Kuu(A12+1,B12)
            amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)

            !
            !  k_uu
            !
            amatrx(A11+2,B11+2) = K_uu(i,j)
            !
            !  k_vv
            !
            amatrx(A11+3,B11+3) = K_vv(i,j)
            !
            !  k_phiphi
            !
            amatrx(A11+4,B11+4) = K_phiphi(i,j)
            !
            !  K_uv
            !
            amatrx(A11+2,B11+3) = K_uv(i,j)
            !
            !  K_uphi
            !
            amatrx(A11+2,B11+4) = K_uphi(i,j)
            !
            !  K_vu
            !
            amatrx(A11+3,B11+2) = K_vu(i,j)
            !
            !  K_vphi
            !
            amatrx(A11+3,B11+4) = K_vphi(i,j)
            !
            !  K_phiu
            !
            amatrx(A11+4,B11+2) = K_phiu(i,j)
            !
            !  K_phiv
            !
            amatrx(A11+4,B11+3) = K_phiv(i,j)
            !
            enddo
         enddo
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine UPE4

************************************************************************

      subroutine UAX4(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,transient)

      use global

      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      real*8 u(nNode,2),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),muNew(nNode)
      real*8 muOld(nNode),dMU(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),u_t(nNode,ndofel),v(nNode,2)
      real*8 coordsC(mcrd,nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,faceFlag

      logical timeStart

      real*8 Iden(3,3),Le,theta0,phi0,Ru(2*nNode,1),Rc(nNode,1),body(3)
      real*8 Kuu(2*nNode,2*nNode),Kcc(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,2),dsh0(nNode,2),dshC0(nNode,2),detMapJ0C,AR0
      real*8 ARc,AR_t,Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),AR
      real*8 sh(nNode),detMapJ,phi_t,dsh(nNode,2),detMapJC,phiLmt,umeror
      real*8 dshC(nNode,2),mu_tau,mu_t,dMUdX(2,1),dMUdt,F_tau(3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,2),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3),phi_tau,dPdt,DphiDmu,DphidotDmu,Mfluid
      real*8 Smat(3,1),BodyForceRes(2*nNode,1),Vmol,SpCUMod(3,3,3),DmDJ
      real*8 SmatAx(4,1),BodyForceResAx(2*nNode,1),dTRdF(3,3,3,3),DmDmu
      real*8 BmatAx(4,2*nNode),Gmat(4,2*nNode),G0mat(4,2*nNode),flux,ds
      real*8 Amat(4,4),Qmat(4,4),AmatAx(5,5),QmatAx(5,5),yLocal(nIntS)
      real*8 G0matAx(5,2*nNode),GmatAx(5,2*nNode),xLocal(nIntS),detF_t
      real*8 wS(nIntS),Kuc(2*nNode,nNode),Kcu(nNode,2*nNode),ResFac
      real*8 TanFac,Nvec(1,nNode),AmatUC(4,1),SpCUModFac(3,3),transient
      real*8 AmatCU(2,5),SpUCMod(3,3),ARcrit,Kbc_t,Kbc,tRelax,mu0
      real*8 timeContact,timeContact_t


      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)


      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point


c$$$      ! Allocate memory for the globalSdv's used for viewing
c$$$      !  results on the dummy mesh
c$$$      !
c$$$      if(.not.allocated(globalSdv)) then
c$$$         !
c$$$         ! allocate memory for the globalSdv's
c$$$         !
c$$$         ! numElem needs to be set in the MODULE
c$$$         ! nInt needs to be set in the UEL
c$$$         !
c$$$         stat=0
c$$$c         allocate(globalSdv(numElem,nInt,ngSdv))
c$$$c         deallocate(globalSdv)
c$$$         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
c$$$         if(stat.ne.0) then
c$$$            write(*,*) '//////////////////////////////////////////////'
c$$$            write(*,*) 'error when allocating globalSdv'
c$$$            write(*,*) '//////////////////////////////////////////////'
c$$$            write(*,*) '   stat=',stat
c$$$            write(*,*) '  ngSdv=',ngSdv
c$$$            write(*,*) '   nInt=',nInt
c$$$            write(*,*) 'numElem=',numElem
c$$$            write(*,*) '  nNode=',nNode
c$$$            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
c$$$            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
c$$$            write(*,*) '//////////////////////////////////////////////'
c$$$            write(80,*) '//////////////////////////////////////////////'
c$$$            write(80,*) 'error when allocating globalSdv'
c$$$            write(80,*) '//////////////////////////////////////////////'
c$$$            write(80,*) '   stat=',stat
c$$$            write(80,*) '  ngSdv=',ngSdv
c$$$            write(80,*) '   nInt=',nInt
c$$$            write(80,*) 'numElem=',numElem
c$$$            write(80,*) '  nNode=',nNode
c$$$            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
c$$$            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
c$$$            write(80,*) '//////////////////////////////////////////////'
c$$$            call xit
c$$$         endif
c$$$         write(*,*) '-------------------------------------------------'
c$$$         write(*,*) '----------- globalSDV ALLOCATED -----------------'
c$$$      if((kstep.eq.1).and.(kinc.eq.1)) then
c$$$         write(*,*) '-------------------------------------------------'
c$$$         write(*,*) '----------YOU PUT NUMBER OF ELEMENTS -----------'
c$$$         write(*,*) '---------- numElem=',numElem
c$$$         write(*,*) '---------- UAX4 ELEMENTS ------------------------'
c$$$         write(*,*) '-------------------------------------------------'
c$$$         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
c$$$         write(*,*) '---------- nInt= ',nInt
c$$$         write(*,*) '---------- nIntS=',nIntS
c$$$         write(*,*) '-------------------------------------------------'
c$$$         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
c$$$         write(*,*) '---------- ngSdv=',ngSdv
c$$$         write(*,*) '-------------------------------------------------'
c$$$      endif
c$$$      if(kinc.eq.1) then
c$$$         if(transient.ne.zero) then
c$$$            write(*,*) '------- Step=',kstep
c$$$            write(*,*) '---------------------------------------------'
c$$$            write(*,*) '----------------- TRANSIENT -----------------'
c$$$         else
c$$$            write(*,*) '--------------- STEADY STATE ----------------'
c$$$         endif
c$$$      endif
c$$$         write(*,*) '-------------------------------------------------'
c$$$      endif


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain initial conditions
      !
      theta0 = props(8)
      phi0   = props(9)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rc = zero
      Kuu = zero
      Kcc = zero
      Kuc = zero
      Kcu = zero
      Energy = zero


      ! Body forces
      !
      body(1:3) = zero


      ! Obtain nodal displacements and chemical potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         k = k + 1
         muNew(i) = Uall(k)
         dMU(i) = DUall(k,1)
         muOld(i) = muNew(i) - dMU(i)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose any time-stepping changes on the increments of
      !  chemical potential or displacement if you wish
      !
      ! chemical potential increment
      !
      do i=1,nNode
         if(dabs(dMU(i)).gt.1.d6) then
            pnewdt = 0.5
            return
         endif
      enddo
      !
      ! displacement increment, based on element diagonal
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,3))**two) +
     +     ((coordsC(2,1)-coordsC(2,3))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo



      !----------------------------------------------------------------
      !
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the
      !  `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! Map shape functions from local to global current coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         write(80,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif


      ! For an axisymmetric problem, find the ``r'' that
      !  shows up in the integrals for axisymmetric
      !  and in the F(3,3), the factors of 2Pi are for integrals
      !  i.e., dV = 2 pi r dr dz
      !
      AR0  = zero
      ARc  = zero
      AR_t = zero
      do i=1,nNode
         ! radial coord in ref config at centroid
         AR0  = AR0 + sh0(i)*coords(1,i)
         ! radial coord in current config at centroid
         ARc  = ARc + sh0(i)*(coords(1,i) + u(i,1))
         ! radial coord in current config at centroid in previous step
         AR_t = AR_t + sh0(i)*(coords(1,i) + uOld(i,1))
      enddo



      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in
      !  the `F-bar' method. `Tau' represents the end of the increment
      !  and `t' the previous increment.
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      !
      ! modify for axisymmetric
      !
      Fc_tau(3,3) = ARc/AR0
      Fc_t(3,3) = AR_t/AR0
      !
      ! axisymmetric implementation detF
      !
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over body integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         elseif(nInt.eq.1) then
            call xint2D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            write(80,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         write(80,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif



      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions
            !
            phi_t  = phi0
            !
         else
            !
            ! this is not the first increment, read old values
            !
            phi_t  = svars(1+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            write(80,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! Map shape functions from local to global current coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            write(80,*) 'Unexpected error, mcrd=',mcrd
            call xit
         endif


         ! For an axisymmetric problem, find the ``r'' that
         !  shows up in the integrals for axisymmetric
         !  and in the F(3,3), the factors of 2Pi are for integrals
         !  i.e., dV = 2 pi r dr dz
         !
         !
         AR0  = zero
         AR   = zero
         AR_t = zero
         do i=1,nNode
            AR0 = AR0 + sh(i)*coords(1,i)
            AR  = AR  + sh(i)*(coords(1,i) + u(i,1))
            AR_t = AR_t + sh(i)*(coords(1,i) + uOld(i,1))
         enddo
         AR0  = two*Pi*AR0
         AR   = two*Pi*AR
         AR_t = two*Pi*AR_t


         ! Obtain the chemical potential and its derivative's at
         !  this intPt at the begining and end of the incrment
         !
         mu_tau = zero
         mu_t = zero
         dMUdt = zero
         dMUdX = zero
         do k=1,nNode
            mu_tau = mu_tau + muNew(k)*sh(k)
            mu_t   = mu_t + muOld(k)*sh(k)
            do i=1,nDim
               dMUdX(i,1) = dMUdX(i,1) + muNew(k)*dshC(k,i)
            enddo
         enddo
         dMUdt = (mu_tau - mu_t)/dtime



         ! Obtain, and modify the deformation gradient at this integration
         !  point.  Modify the deformation gradienet for use in the `F-bar'
         !  method.  Also, take care of plane-strain or axisymetric
         !
         F_tau = Iden
         F_t = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
                  F_t(i,j) = F_t(i,j) + dsh(k,j)*uOld(k,i)
               enddo
            enddo
         enddo
         !
         ! modify F(3,3) for axisymetric, give R/R0
         !
         F_tau(3,3) = AR/AR0
         F_t(3,3) = AR_t/AR0
         !
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 4 node fully integrated linear
         !  element, do not use the `F-bar' method for any other element
         !
         if((nNode.eq.4).and.(nInt.eq.4)) then
            call mdet(F_tau,detF_tau)
            call mdet(F_t,detF_t)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
            F_t = ((detFc_tau/detF_tau)**third)*F_t
         endif
         call mdet(F_tau,detF)


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point
         !
         call integ(props,nprops,dtime,
     +        F_tau,mu_tau,phi_t,theta0,
     +        T_tau,SpTanMod,
     +        phi_tau,dPdt,DphiDmu,DphidotDmu,
     +        Mfluid,DmDmu,DmDJ,Vmol,
     +        SpUCMod,SpCUModFac)
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = phi_tau
         jj = jj + nlSdv ! setup for the next intPt


         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = phi_tau   ! polymer volume fraction
         globalSdv(jelem,intPt,2) = T_tau(1,1) ! radial stress
         globalSdv(jelem,intPt,3) = T_tau(3,3) ! hoop stress
         globalSdv(jelem,intPt,4) = -(one/three)*
     +        (T_tau(1,1)+T_tau(2,2)+T_tau(3,3)) ! mean normal pressure



         ! Time stepping algorithim based on the constitutive response.
         !  Here based on the change in the polymer volume fraction change.
         !
         phiLmt = 0.005d0
         umeror = dabs((phi_tau - phi_t)/phiLmt)
         if(umeror.le.half) then
            pnewdt = 1.5d0
         elseif((umeror.gt.half).and.(umeror.le.0.8d0)) then
            pnewdt = 1.25d0
         elseif((umeror.gt.0.8d0).and.(umeror.le.1.25d0)) then
            pnewdt = 0.75d0
         else
            pnewdt = half
         endif


         ! Compute/update the displacement residual vector
         !
         SmatAx(1,1) = T_tau(1,1)
         SmatAx(2,1) = T_tau(2,2)
         SmatAx(3,1) = T_tau(1,2)
         SmatAx(4,1) = T_tau(3,3)
         !
         BmatAx = zero
         do kk=1,nNode
            BmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(2,2+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            BmatAx(3,2+nDim*(kk-1)) = dshC(kk,1)
            BmatAx(4,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo
         !
         BodyForceResAX = zero
         do kk=1,nNode
            BodyForceResAx(1+nDim*(kk-1),1) = sh(kk)*body(1)
            BodyForceResAx(2+nDim*(kk-1),1) = sh(kk)*body(2)
         enddo
         !
         Ru = Ru + detmapJC*w(intpt)*AR*
     +        (
     +        -matmul(transpose(BmatAx),SmatAx)
     +        + BodyForceResAx
     +        )


         ! Compute/update the chemical potential residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo
         !
         ResFac = (transient*dPdt)/(detF*Vmol*phi_tau*phi_tau)
         !
         Rc = Rc + detmapJC*w(intpt)*AR*
     +        (
     +        transpose(Nvec)*ResFac - Mfluid*matmul(dshC,dMUdX)
     +        )



         ! Compute/update the displacement tangent matrix
         !
         GmatAx = zero
         do kk=1,nNode
            GmatAx(1,1+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(2,2+nDim*(kk-1)) = dshC(kk,1)
            GmatAx(3,1+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(4,2+nDim*(kk-1)) = dshC(kk,2)
            GmatAx(5,1+nDim*(kk-1)) = sh(kk)/(AR/(two*Pi))
         enddo

         G0matAx = zero
         do kk=1,nNode
            G0matAx(1,1+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(2,2+nDim*(kk-1)) = dshC0(kk,1)
            G0matAx(3,1+nDim*(kk-1)) = dshC0(kk,2)
            G0matAx(4,2+nDim*(kk-1)) = dshC0(kk,2)
            G0matAX(5,1+nDim*(kk-1)) = sh0(kk)/ARc
         enddo

         AmatAx = zero
         AmatAx(1,1) = SpTanMod(1,1,1,1)
         AmatAx(1,2) = SpTanMod(1,1,2,1)
         AmatAx(1,3) = SpTanMod(1,1,1,2)
         AmatAx(1,4) = SpTanMod(1,1,2,2)
         AmatAx(1,5) = SpTanMod(1,1,3,3)
         AmatAx(2,1) = SpTanMod(2,1,1,1)
         AmatAx(2,2) = SpTanMod(2,1,2,1)
         AmatAx(2,3) = SpTanMod(2,1,1,2)
         AmatAx(2,4) = SpTanMod(2,1,2,2)
         AmatAx(2,5) = SpTanMod(2,1,3,3)
         AmatAx(3,1) = SpTanMod(1,2,1,1)
         AmatAx(3,2) = SpTanMod(1,2,2,1)
         AmatAx(3,3) = SpTanMod(1,2,1,2)
         AmatAx(3,4) = SpTanMod(1,2,2,2)
         AmatAx(3,5) = SpTanMod(1,2,3,3)
         AmatAx(4,1) = SpTanMod(2,2,1,1)
         AmatAx(4,2) = SpTanMod(2,2,2,1)
         AmatAx(4,3) = SpTanMod(2,2,1,2)
         AmatAx(4,4) = SpTanMod(2,2,2,2)
         AmatAx(4,5) = SpTanMod(2,2,3,3)
         AmatAx(5,1) = SpTanMod(3,3,1,1)
         AmatAx(5,2) = SpTanMod(3,3,2,1)
         AmatAx(5,3) = SpTanMod(3,3,1,2)
         AmatAx(5,4) = SpTanMod(3,3,2,2)
         AmatAx(5,5) = SpTanMod(3,3,3,3)

         QmatAx = zero
         QmatAx(1,1) = third*(AmatAx(1,1)+AmatAx(1,4)+AmatAx(1,5))
     +        - (two/three)*T_tau(1,1)
         QmatAx(2,1) = third*(AmatAx(2,1)+AmatAx(2,4)+AmatAx(2,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(3,1) = third*(AmatAx(3,1)+AmatAx(3,4)+AmatAx(3,5))
     +        - (two/three)*T_tau(1,2)
         QmatAx(4,1) = third*(AmatAx(4,1)+AmatAx(4,4)+AmatAx(4,5))
     +        - (two/three)*T_tau(2,2)
         QmatAx(5,1) = third*(AmatAx(5,1)+AmatAx(5,4)+AmatAx(5,5))
     +        - (two/three)*T_tau(3,3)
         QmatAx(1,4) = QmatAx(1,1)
         QmatAx(2,4) = QmatAx(2,1)
         QmatAx(3,4) = QmatAx(3,1)
         QmatAx(4,4) = QmatAx(4,1)
         QmatAx(5,4) = QmatAx(5,1)
         QmatAx(1,5) = QmatAx(1,1)
         QmatAx(2,5) = QmatAx(2,1)
         QmatAx(3,5) = QmatAx(3,1)
         QmatAx(4,5) = QmatAx(4,1)
         QmatAx(5,5) = QmatAx(5,1)


         if((nNode.eq.4).and.(nInt.eq.4)) then
            !
            ! This is the tangent using the F-bar method with the
            !  4 node fully integrated linear element
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           + matmul(transpose(GmatAx),matmul(QmatAx,
     +           (G0matAx-GmatAx)))
     +           )
         else
            !
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu + detMapJC*w(intpt)*AR*
     +           (
     +           matmul(matmul(transpose(GmatAx),AmatAx),GmatAx)
     +           )
         endif



         ! Compute/update the chemical potential tangent matrix
         !
         TanFac = transient*
     +        (
     +        (one/(detF*Vmol*phi_tau**two))*
     +        (two*(dPdt/phi_tau)*DphiDmu - DphidotDmu)
     +        )
         !
         Kcc = Kcc + detmapJC*AR*w(intPt)*
     +        (
     +        TanFac*matmul(transpose(Nvec),Nvec)
     +        + Mfluid*matmul(dshC,transpose(dshC))
     +        + DmDmu*matmul(matmul(dshC,dMUdX),Nvec)
     +        )



         ! Compute/update the chemical potential - displacement tangent matrix.
         !  The F-bar method will have some effect, however we neglect that here.
         !
         SpCUMod = zero
         do i=1,nDim
            do k=1,nDim
               do l=1,nDim
                  SpCUMod(i,k,l) = SpCUMod(i,k,l)
     +                 + dMUdX(k,1)*SpCUModFac(i,l)
               enddo
            enddo
         enddo
         !
         AmatCU = zero
         AmatCU(1,1) = SpCUMod(1,1,1)
         AmatCU(1,2) = SpCUMod(1,2,1)
         AmatCU(1,3) = SpCUMod(1,1,2)
         AmatCU(1,4) = SpCUMod(1,2,2)
         AmatCU(1,5) = SpCUMod(1,3,3)
         AmatCU(2,1) = SpCUMod(2,1,1)
         AmatCU(2,2) = SpCUMod(2,2,1)
         AmatCU(2,3) = SpCUMod(2,1,2)
         AmatCU(2,4) = SpCUMod(2,2,2)
         AmatCU(2,5) = SpCUMod(2,3,3)
         !
         Kcu = Kcu - detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(dshC,AmatCU),GmatAX)
     +        )


         ! Compute/update the displacement - chemical potential tangent matrix
         !  The F-bar method will have some effect, however we neglect that here.
         !
         AmatUC = zero
         AmatUC(1,1) = SpUCMod(1,1)
         AmatUC(2,1) = SpUCMod(2,2)
         AmatUC(3,1) = SpUCMod(1,2)
         AmatUC(4,1) = SpUCMod(3,3)
         !
         Kuc = Kuc + detMapJC*w(intpt)*AR*
     +        (
     +        matmul(matmul(transpose(BmatAX),AmatUC),Nvec)
     +        )


      enddo
      !
      ! End the loop over body integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface fluid flux terms here
      !
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the flux/traction
            !  acts on is the flux/traction ``label''
            !
            face = jdltyp(i,1) ! label
            flux = adlmag(i,1) ! flux magnitude

            if((face.ge.1).and.(face.le.4)) then
               !
               ! fluid flux applied
               !
               if(face.eq.1) then
                  faceFlag = 1
               elseif(face.eq.2) then
                  faceFlag = 2
               elseif(face.eq.3) then
                  faceFlag = 3
               else
                  faceFlag = 4
               endif
               !
               if(nIntS.eq.1) then
                  call xintSurf2D1pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.2) then
                  call xintSurf2D2pt(faceFlag,xLocal,yLocal,wS)
               elseif(nIntS.eq.3) then
                  call xintSurf2D3pt(faceFlag,xLocal,yLocal,wS)
               else
                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
                  write(80,*) 'Invalid nIntS points, nIntS=',nIntS
                  call xit
               endif
               !
               ! loop over integ points
               !
               do ii=1,nIntS
                  !
                  ! Compute shape functions, derivatives, and the
                  !  mapping jacobian (ds)
                  !
                  call computeSurf(nNode,xLocal(ii),yLocal(ii),faceFlag,
     +                 coordsC,sh,ds)
                  !
                  ! Axisymmetric ``radius''
                  !
                  AR = zero
                  do n=1,nNode
                     AR = AR + sh(n)*coordsC(1,n)
                  enddo
                  AR = two*Pi*AR
                  !
                  ! Modify the residual, loop over nodes, recall
                  !  sh(n)=0 when n is not on this face
                  !
                  do n=1,nNode
                     Rc(n,1) = Rc(n,1) - wS(ii)*ds*sh(n)*flux*AR
                  enddo
                  !
                  ! No change to the tangent matrix
                  !
               enddo ! loop over nIntS
               !
c$$$            elseif((face.ge.91).and.(face.le.94)) then
c$$$               !
c$$$               ! We have that special boundary condition, free swelling
c$$$               !  with mu=mu0, followed by contact and flux=zero
c$$$               !
c$$$               if(transient.lt.0.5d0) then
c$$$                  write(*,*) 'the BC switch is only for transient'
c$$$                  write(80,*) 'the BC switch is only for transient'
c$$$                  call xit
c$$$               endif
c$$$               !
c$$$               if(face.eq.91) then
c$$$                  faceFlag = 1
c$$$                  if(nNode.eq.4) then
c$$$                     k = 1 ! used to determine contact later
c$$$                     m = 2 ! used to determine contact later
c$$$                  elseif(nNode.eq.6) then
c$$$                     k = 1
c$$$                     m = 2
c$$$                     q = 4
c$$$                  endif
c$$$               elseif(face.eq.92) then
c$$$                  faceFlag = 2
c$$$                  if(nNode.eq.4) then
c$$$                     k = 2
c$$$                     m = 3
c$$$                  elseif(nNode.eq.6) then
c$$$                     k = 2
c$$$                     m = 3
c$$$                     q = 5
c$$$                  endif
c$$$               elseif(face.eq.93) then
c$$$                  faceFlag = 3
c$$$                  if(nNode.eq.4) then
c$$$                     k = 3
c$$$                     m = 4
c$$$                  elseif(nNode.eq.6) then
c$$$                     k = 1
c$$$                     m = 3
c$$$                     q = 6
c$$$                  endif
c$$$               else
c$$$                  faceFlag = 4
c$$$                  k = 2
c$$$                  m = 3
c$$$                  if(nNode.eq.6) then
c$$$                     write(*,*) 'error, no face 4'
c$$$                     write(80,*) 'error, no face 4'
c$$$                     call xit
c$$$                  endif
c$$$               endif
c$$$               !
c$$$               if(nIntS.eq.1) then
c$$$                  call xintSurf2D1pt(faceFlag,xLocal,yLocal,wS)
c$$$               elseif(nIntS.eq.2) then
c$$$                  call xintSurf2D2pt(faceFlag,xLocal,yLocal,wS)
c$$$               elseif(nIntS.eq.3) then
c$$$                  call xintSurf2D3pt(faceFlag,xLocal,yLocal,wS)
c$$$               else
c$$$                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
c$$$                  write(80,*) 'Invalid nIntS points, nIntS=',nIntS
c$$$                  call xit
c$$$               endif
c$$$               !
c$$$               ! loop over integ points
c$$$               !
c$$$               do ii=1,nIntS
c$$$
c$$$                  ! Obtain state variables from previous increment
c$$$                  !
c$$$                  if((kinc.le.1).and.(kstep.eq.1)) then
c$$$                  !
c$$$                  ! this is the first increment, of the first step
c$$$                  !  give initial conditions (or just anything)
c$$$                  !
c$$$                     Kbc_t = 1.d3
c$$$                  !
c$$$                  else
c$$$                  !
c$$$                  ! this is not the first increment, read old values
c$$$                  !
c$$$                     Kbc_t  = svars(nlSdv*nInt + 1)
c$$$                  !
c$$$                  endif
c$$$
c$$$
c$$$                  ! Compute shape functions, derivatives, and the
c$$$                  !  mapping jacobian (ds)
c$$$                  !
c$$$                  call computeSurf(nNode,xLocal(ii),yLocal(ii),faceFlag,
c$$$     +                 coordsC,sh,ds)
c$$$
c$$$
c$$$                  ! Axisymmetric ``radius''
c$$$                  !
c$$$                  AR = zero
c$$$                  do n=1,nNode
c$$$                     AR = AR + sh(n)*coordsC(1,n)
c$$$                  enddo
c$$$                  AR = two*Pi*AR
c$$$
c$$$
c$$$                  ! The actual boundary condition prior to contact
c$$$                  !
c$$$                  mu0 = flux
c$$$
c$$$
c$$$                  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$                  ! This only works for axisymmetric contact along the right side
c$$$                  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$                  !
c$$$                  !
c$$$                  ARcrit = props(10)  ! this is the wall location
c$$$                  tRelax = props(11)  ! relaxation time on the BC
c$$$
c$$$
c$$$
c$$$                  ! Check long before contact will be made so that the
c$$$                  !  bulk modulus made be reduced prior to contact
c$$$                  !
c$$$                  timeStart = .false.
c$$$                  !
c$$$                  if(nNode.eq.4) then
c$$$                     !
c$$$                     ! 2 nodes per face on the 4 node quadrillateral
c$$$                     !
c$$$                     if((coordsC(1,k).ge.ARcrit-0.5*Le).and.
c$$$     +                    (coordsC(1,m).ge.ARcrit-0.5*Le)) then
c$$$                        timeStart = .true.
c$$$                     endif
c$$$                     !
c$$$                  elseif(nNode.eq.6) then
c$$$                     !
c$$$                     ! 3 nodes per face on the 6 node triangle
c$$$                     !
c$$$                     if((coordsC(1,k).ge.ARcrit-0.5*Le).and.
c$$$     +                    (coordsC(1,m).ge.ARcrit-0.5*Le).and.
c$$$     +                    (coordsC(1,q).ge.ARcrit-0.5*Le)) then
c$$$                        timeStart = .true.
c$$$                     endif
c$$$                     !
c$$$                  endif
c$$$                  !
c$$$                  if((timeStart).and.(timeContact_t.lt.zero)) then
c$$$                     !
c$$$                     ! just made contact with the wall, save this time
c$$$                     !
c$$$                     timeContact = time(2)
c$$$                  else
c$$$                     timeContact = timeContact_t
c$$$                  endif
c$$$
c$$$
c$$$
c$$$                  ! Check if this surface is in contact with the wall
c$$$                  !  (this is really for the whole element, not for
c$$$                  !  each integration point on the surface)
c$$$                  !
c$$$                  ! this is solving the differential equation implicitly for
c$$$                  !  exponential decay, such that Kbc=Kbc0*exp(-t/tRelax)
c$$$                  !
c$$$                  if(nNode.eq.4) then
c$$$                     !
c$$$                     ! 2 nodes per face on the 4 node quadrillateral
c$$$                     !
c$$$                     if((coordsC(1,k).ge.ARcrit-1.d-6).and.
c$$$     +                    (coordsC(1,m).ge.ARcrit-1.d-6)) then
c$$$                        Kbc = Kbc_t/(one + (dtime/tRelax))
c$$$                     else
c$$$                        Kbc = Kbc_t
c$$$                     endif
c$$$                     !
c$$$                  elseif(nNode.eq.6) then
c$$$                     !
c$$$                     ! 3 nodes per face on the 6 node triangle
c$$$                     !
c$$$                     if((coordsC(1,k).ge.ARcrit-1.d-6).and.
c$$$     +                    (coordsC(1,m).ge.ARcrit-1.d-6).and.
c$$$     +                    (coordsC(1,q).ge.ARcrit-1.d-6)) then
c$$$                        Kbc = Kbc_t/(one + (dtime/tRelax))
c$$$                     else
c$$$                        Kbc = Kbc_t
c$$$                     endif
c$$$                     !
c$$$                  endif
c$$$                  if(Kbc.lt.1.d0) Kbc = zero
c$$$
c$$$
c$$$
c$$$
c$$$                  ! Save state variables
c$$$                  !
c$$$                  svars(nlSdv*nInt + 1) = Kbc
c$$$                  svars(nlSdv*nInt + 2) = timeContact
c$$$                  !
c$$$                  !
c$$$                  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$
c$$$
c$$$                  ! Obtain the chemical potential at this integ point
c$$$                  !
c$$$                  mu_tau = zero
c$$$                  do n=1,nNode
c$$$                     mu_tau = mu_tau + muNew(n)*sh(n)
c$$$                  enddo
c$$$
c$$$
c$$$                  ! Modify the residual
c$$$                  !
c$$$                  do n=1,nNode
c$$$                     Rc(n,1) = Rc(n,1) +
c$$$     +                    wS(ii)*AR*ds*sh(n)*Kbc*(mu_tau - mu0)
c$$$                  enddo
c$$$
c$$$
c$$$                  ! Modify the chemical potential tangent
c$$$                  !
c$$$                  do n=1,nNode
c$$$                     do m=1,nNode
c$$$                        Kcc(n,m) = Kcc(n,m) -
c$$$     +                       wS(ii)*AR*ds*sh(n)*sh(m)*Kbc
c$$$                     enddo
c$$$                  enddo
c$$$
c$$$               enddo ! loop over nIntS
c$$$               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over flux and traction terms
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.
      !
      call AssembleElement(nDim,nNode,nDofEl,
     +     Ru,Rc,Kuu,Kuc,Kcu,Kcc,
     +     rhs,amatrx)
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine UAX4

************************************************************************

      subroutine U3D8(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     +     PROPS,NPROPS,coords,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     +     TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,
     +     PREDEF,NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,
     +     NJPROP,PERIOD,
     +     nDim,transient)

      use global

      IMPLICIT NONE
*
*     VARIABLES DEFINED IN UEL, PASSED BACK TO ABAQUS
*
      REAL(8) :: RHS,AMATRX,SVARS,ENERGY
*
*     VARIABLES PASSED INTO UEL
*
      REAL(8) :: PROPS,coords,Uall,DUall,Vel,Accn,TIME,
     1  DTIME,PARAMS,ADLMAG,PREDEF,DDLMAG,PNEWDT,PERIOD
      INTEGER :: NDOFEL,NRHS,NSVARS,NPROPS,MCRD,NNODE,JTYPE,KSTEP,KINC,
     1  JELEM,NDLOAD,JDLTYP,NPREDF,LFLAGS,MLVARX,MDLOAD,JPROPS,NJPROP
*
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),coords(MCRD,NNODE),Uall(NDOFEL),
     2 DUall(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)


      real*8 u(nNode,3),du(nNode,ndofel),thetaNew(nNode)
      real*8 thetaOld(nNode),dtheta(nNode),muNew(nNode)
      real*8 muOld(nNode),dMU(nNode),uNew(nNode,ndofel)
      real*8 uOld(nNode,ndofel),v(nNode,3)
      real*8 coordsC(mcrd,nNode)

      integer i,j,k,l,m,n,nIntPt,nDim,intpt,pOrder,a1,b1,a11,b11,face
      integer ii,jj,a12,b12,pe,nSdv,stat,q,nIntV,nIntPtV,p,ngSdv
      integer nlSdv,kk,lenJobName,lenOutDir,faceFlag,nDofN

      real*8 Iden(3,3),Le,theta0,phi0,Ru(3*nNode,1),Rc(nNode,1),body(3)
      real*8 Kuu(3*nNode,3*nNode),Kcc(nNode,nNode),sh0(nNode),detMapJ0
      real*8 dshxi(nNode,3),dsh0(nNode,3),dshC0(nNode,3),detMapJ0C,Vmol
      real*8 Fc_tau(3,3),Fc_t(3,3),detFc_tau,detFc_t,w(nInt),DmDmu,DmDJ
      real*8 sh(nNode),detMapJ,dsh(nNode,3),detMapJC,phiLmt,umeror
      real*8 dshC(nNode,3),mu_tau,mu_t,dMUdX(3,1),dMUdt,F_tau(3,3)
      real*8 F_t(3,3),detF_tau,xi(nInt,3),detF,TR_tau(3,3),T_tau(3,3)
      real*8 SpTanMod(3,3,3,3),dPdt,DphiDmu,DphidotDmu,Mfluid
      real*8 Smat(6,1),Bmat(6,3*nNode),BodyForceRes(3*nNode,1),flux
      real*8 Gmat(9,3*nNode),G0mat(9,3*nNode),Amat(9,9),Qmat(9,9),dA
      real*8 xLocal(nIntS),yLocal(nIntS),zLocal(nIntS),wS(nIntS),detF_t
      real*8 Kuc(3*nNode,nNode),Kcu(nNode,3*nNode),Nvec(1,nNode),ResFac
      real*8 AmatUC(6,1),TanFac,AmatCU(3,9),SpUCMod(3,3),SpCUMod(3,3,3)
      real*8 SpCUModFac(3,3),transient
      real*8 R_u(nNode,1),R_v(nNode,1),R_phi(nNode,1)
      real*8 K_uu(nNode,nNode),K_vv(nNode,nNode),K_phiphi(nNode,nNode)
      real*8 K_uv(nNode,nNode),K_uphi(nNode,nNode),K_vu(nNode,nNode)
      real*8 K_vphi(nNode,nNode),K_phiu(nNode,nNode),K_phiv(nNode,nNode)
      real*8 actNew(nNode),dact(nNode),actOld(nNode)
      real*8 supNew(nNode),dsup(nNode),supOld(nNode)
      real*8 phiNew(nNode),dphi(nNode),phiOld(nNode)
      real*8 c_t,c_tau
      real*8 u_tau,u_t,dudt,dudX(nDim,1)
      real*8 v_tau,v_t,dvdt,dvdX(nDim,1)
      real*8 phi_tau,phi_t,dphidt,dphidX(nDim,1)
      real*8 Diff_u,Diff_v,F_const,k_const,sigma,permittivity
      real*8 alpha,beta







      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)

      character*256 jobName,outDir,fileName

      ! Get element parameters
      !
      nlSdv  = jprops(1) ! number of local sdv's per integ point
      ngSdv  = jprops(2) ! number of global sdv's per integ point




      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero   !  residual of displacement
      R_u = zero   !  residual of activator
      R_v = zero   !  residual of suppressor
      R_phi = zero !  residual of electric potenential


      ! diagonal
      Kuu  = zero     ! tangent of displacement
      K_uu = zero     ! tangent of activator
      K_vv = zero     ! tangent of suppressor
      K_phiphi = zero ! tangent of electric potential


      ! off-diagonal
      K_uv     = zero   ! tangent of activator - suppressor
      K_uphi   = zero   ! tangent of activator - potential
      K_vu     = zero   ! tangent of suppressor - activator
      K_vphi   = zero   ! tangent of suppressor - potential
      K_phiu   = zero   ! tangent of potential - activator
      K_phiv   = zero   ! tangent of potential - suppressor


      Energy = zero





      ! Obtain nodal degree of freedom
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
            uOld(i,j) = u(i,j) - du(i,j)
         enddo
         ! the activator
         k = k + 1
         actNew(i) = Uall(k)
         dact(i) = DUall(k,1)
         actOld(i) = actNew(i) - dact(i)
         ! the suppressor
         k = k + 1
         supNew(i) = Uall(k)
         dsup(i) = DUall(k,1)
         supOld(i) = supNew(i) - dsup(i)
         ! the electric potential
         k = k + 1
         phiNew(i) = Uall(k)
         dphi(i) = DUall(k,1)
         phiOld(i) = phiNew(i) - dphi(i)
      enddo



      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


c      ! Impose any time-stepping changes on the increments of
c      !  chemical potential or displacement if you want
c      !
c      ! chemical potential increment
c      !
c      do i=1,nNode
c         if(dabs(dMU(i)).gt.1.d6) then
c            pnewdt = 0.5
c            return
c         endif
c      enddo
      !
      ! displacement increment, based on element diagonal
      !
c      Le = dsqrt(((coordsC(1,1)-coordsC(1,7))**two) +
c     +     ((coordsC(2,1)-coordsC(2,7))**two) +
c     +     ((coordsC(3,1)-coordsC(3,7))**two))
c      !
c      do i=1,nNode
c         do j=1,nDim
c            if(dabs(du(i,j)).gt.10.d0*Le) then
c               pnewdt = 0.5
c               return
c            endif
c         enddo
c      enddo



      !----------------------------------------------------------------
      !
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Here, check for hourglass stabilization and get
      !  the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures, 33, 3277-3296.
      !
      !
      ! Obtain shape functions and their local gradients at the element
      !  centriod, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         write(80,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif



      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif




      ! Calculate the deformation gradient at the element centriod
      !  at the the begining and end of the increment for use in
      !  the `F-bar' method
      !
      Fc_tau = Iden
      Fc_t = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
               Fc_t(i,j) = Fc_t(i,j) + dsh0(k,j)*uOld(k,i)
            enddo
         enddo
      enddo
      call mdet(Fc_tau,detFc_tau)
      call mdet(Fc_t,detFc_t)
      !
      ! With the deformation gradient known at the element centriod
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nInt.eq.1) then
         call xint3D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
      elseif(nInt.eq.8) then
         call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
      else
         write(*,*) 'Invalid number of int points, nInt=',nInt
         write(80,*) 'Invalid number of int points, nInt=',nInt
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions (or just anything)
            !
            c_t  = 0.00001
            !
         else
            !
            ! this is not the first increment, read old values
            !
            c_t  = svars(1+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            write(80,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Obtain activator and its derivative's at
         !  this intPt at the begining and end of the incrment
         !
         u_tau = zero
         u_t = zero
         dudt = zero
         dudX = zero
         do k=1,nNode
            u_tau = u_tau + actNew(k)*sh(k)
            u_t   = u_t + actOld(k)*sh(k)
            do i=1,nDim
               dudX(i,1) = dudX(i,1) + actNew(k)*dsh(k,i)
            enddo
         enddo
         dudt = (u_tau - u_t)/dtime



         ! Obtain the species v and its derivative's at
         !  this intPt at the begining and end of the incrment
         !
         v_tau = zero
         v_t = zero
         dvdt = zero
         dvdX = zero
         do k=1,nNode
            v_tau = v_tau + supNew(k)*sh(k)
            v_t   = v_t + supOld(k)*sh(k)
            do i=1,nDim
               dvdX(i,1) = dvdX(i,1) + supNew(k)*dsh(k,i)
            enddo
         enddo
         dvdt = (v_tau - v_t)/dtime



         ! Obtain the electric potential and its derivative's at
         !  this intPt at the begining and end of the incrment
         !
         phi_tau = zero
         phi_t = zero
         dphidt = zero
         dphidX = zero
         do k=1,nNode
            phi_tau = phi_tau + phiNew(k)*sh(k)
            phi_t   = phi_t + phiOld(k)*sh(k)
            do i=1,nDim
               dphidX(i,1) = dphidX(i,1) + phiNew(k)*dsh(k,i)
            enddo
         enddo
         dphidt = (phi_tau - phi_t)/dtime



         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! I Stopped here, I should update the integ later!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive time integration at this integ. point
         !
         call integ(props,nprops,dtime,
     +     u_tau,u_t,v_tau,v_t,
     +     Diff_u,Diff_v,F_const,k_const,sigma,permittivity,
     +     c_t,c_tau)
         !
         ! I need sigma, permittivity
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = c_tau
         jj = jj + nlSdv ! setup for the next intPt



         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = c_tau
         globalSdv(jelem,intPt,2) = - dphidX(1,1) ! x component of electric field
         globalSdv(jelem,intPt,3) = - dphidX(2,1) ! y component of electric field
         globalSdv(jelem,intPt,4) = - dphidX(3,1) ! y component of electric field



         ! calculate some lengthy quantities
         alpha = - u_tau*v_tau**2.0 + F_const*(1.0 - u_tau)
         beta  = u_tau*v_tau**2.0 - (F_const+k_const)*v_tau


         ! Compute/update the activator residual vector
         !
         do kk=1,nNode
            Nvec(1,kk) = sh(kk)
         enddo



         R_u = R_u +  detmapJ*w(intpt)*
     +         (
     +         transpose(Nvec)*dudt
     +       + matmul(dsh,dudX)*Diff_u
     +       - transpose(Nvec)*alpha
     +         )



         ! Compute/update the suppressor residual vector
         !


         R_v = R_v + detmapJ*w(intpt)*
     +         (
     +         transpose(Nvec)*dvdt
     +       + Diff_v*matmul(dsh,dvdX)
     +       - transpose(Nvec)*beta
     +       - sigma*matmul(transpose(Nvec),matmul(transpose(dphidX),dvdX))
     +         )


         ! Compute/update the eletric residual vector
         !
         ! I need permittivity
c         R_phi = R_phi + detmapJ*w(intpt)*
c     +         (
c     +         - permittivity*matmul(dsh,dphidX)
c     +         )

         R_phi = R_phi + detmapJ*w(intpt)*
     +         (
     +         matmul(dsh,dphidX)
     +         )




         ! Compute/update the K_uu
         !


         K_uu = K_uu + detmapJ*w(intPt)*
     +          (
     +          -matmul(transpose(Nvec),Nvec)*(1.0/dtime)
     +          -matmul(dsh,transpose(dsh))*Diff_u
     +          -matmul(transpose(Nvec),Nvec)*(v_tau**2.0 + F_const)
     +          )



         ! Compute/update the K_vv
         !
         !
c         K_vv = K_vv + detmapJ*w(intPt)*
c     +            (
c     +            -(1.0/dtime)*matmul(transpose(Nvec),Nvec)
c     +            -Diff_v*matmul(dsh,transpose(dsh))
c     +            +matmul(transpose(Nvec),Nvec)*(2.0*u_tau*v_tau-F_const-k_const)
c     +            +sigma*matmul(matmul(transpose(Nvec),transpose(dphidX)),transpose(dsh))
c     +            )

         K_vv = K_vv + detmapJ*w(intPt)*
     +            (
     +            -(1.0/dtime)*matmul(transpose(Nvec),Nvec)
     +            -Diff_v*matmul(dsh,transpose(dsh))
     +            +matmul(transpose(Nvec),Nvec)*(2.0*u_tau*v_tau-F_const-k_const)
     +            +sigma*matmul(matmul(transpose(Nvec),transpose(dphidX)),transpose(dsh))
     +            )



         ! Compute/update the K_phiphi
         !
         !
c         K_phiphi = K_phiphi + detmapJ*w(intPt)*
c     +            (
c     +            +permittivity*matmul(dsh,transpose(dsh))
c     +            )

         K_phiphi = K_phiphi + detmapJ*w(intPt)*
     +            (
     +            -matmul(dsh,transpose(dsh))
     +            )


         ! Compute/update the K_uv
         !
         !
c         K_uv  = K_uv  + detmapJ*w(intPt)*
c     +           (
c     +           -matmul(transpose(Nvec),Nvec)*2.0*u_tau*v_tau
c     +           )

         K_uv  = K_uv  + detmapJ*w(intPt)*
     +           (
     +           -matmul(transpose(Nvec),Nvec)*2.0*u_tau*v_tau
     +           )


         ! Compute/update the K_uphi
         !
         !
         K_uphi  = zero





         ! Compute/update the K_vu
         !
         !
c         K_vu = K_vu + detmapJ*w(intPt)*
c     +           (
c     +           +(v_tau**2.0)*matmul(transpose(Nvec),Nvec)
c     +           )

         K_vu = K_vu + detmapJ*w(intPt)*
     +           (
     +           +(v_tau**2.0)*matmul(transpose(Nvec),Nvec)
     +           )

         ! Compute/update the K_vphi
         !
         !
c         K_vphi  = K_vphi  + detmapJ*w(intPt)*
c     +           (
c     +           +sigma*matmul(matmul(transpose(Nvec),transpose(dvdX)),transpose(dsh))
c     +           )
         K_vphi  = K_vphi  + detmapJ*w(intPt)*
     +           (
     +           +sigma*matmul(matmul(transpose(Nvec),transpose(dvdX)),transpose(dsh))
     +           )


         K_phiu = zero
         K_phiv = zero


         ! The displacement residual and tangent
         !
c         do i=1,nNode
c            Ru(i,1) = zero
c            Kuu(i,i) = one
c         enddo
         Ru = zero
         Kuu = one



      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Start loop over surface flux terms
      !
      if(ndload.gt.0) then
         !
         ! loop over faces and make proper modifications to
         !  residuals and tangents if needed
         !
         do i=1,ndload
            !
            ! based on my convention the face which the flux
            !  acts on is the flux ``label''
            !
            face = jdltyp(i,1)
            flux = adlmag(i,1)


            if((face.ge.1).and.(face.le.6)) then
               !
               ! fluid flux applied
               !
               if(face.eq.1) then
                  faceFlag = 1
               elseif(face.eq.2) then
                  faceFlag = 2
               elseif(face.eq.3) then
                  faceFlag = 3
               elseif(face.eq.4) then
                  faceFlag = 4
               elseif(face.eq.5) then
                  faceFlag = 5
               else
                  faceFlag = 6
               endif
               !
               if(nIntS.eq.1) then
                  call xintSurf3D1pt(faceFlag,xLocal,yLocal,zLocal,wS)
               elseif(nIntS.eq.4) then
                  call xintSurf3D4pt(faceFlag,xLocal,yLocal,zLocal,wS)
               else
                  write(*,*) 'Invalid nIntS points, nIntS=',nIntS
                  write(80,*) 'Invalid nIntS points, nIntS=',nIntS
                  call xit
               endif
               !
               ! loop over integ points on this element face
               !
               do ii=1,nIntS

                  ! Compute shape functions, derivatives, and the
                  !  mapping jacobian (dA)
                  !
                  call computeSurf3D(xLocal(ii),yLocal(ii),zLocal(ii),
     +                 faceFlag,coordsC,sh,dA)
                  !
                  ! Modify the chemical potential residual, loop over nodes
                  !
                  do n=1,nNode
                     Rc(n,1) = Rc(n,1) - wS(ii)*dA*sh(n)*flux
                  enddo
                  !
                  ! No change to the tangent matrix
                  !
               enddo ! end loop over integ points
               !
            else
               write(*,*) 'Unknown face=',face
               write(80,*) 'Unknown face=',face
               call xit
            endif

         enddo ! loop over ndload
      endif ! ndload.gt.0 or not
      !
      ! End loop over surface flux terms
      !----------------------------------------------------------------





      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the tangent matrix
      !

      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode

!      write(*,*) 'nDofN=',nDofN

      ! The residual vector (RHS)
      !
      do i=1,nNode
         A11 = nDofN*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! displacement
         !
         rhs(A11,1) = Ru(A12,1)
         rhs(A11+1,1) = Ru(A12+1,1)
         rhs(A11+2,1) = Ru(A12+2,1)
         !
         ! activator
         !
         rhs(A11+3,1) = R_u(i,1)
         !
         ! suppressor
         !
         rhs(A11+4,1) = R_v(i,1)
         !
         ! potential
         !
         rhs(A11+5,1) = R_phi(i,1)
         !
         !
      enddo



      ! Assemble the element tangent matrix (AMATRX)
      !
      do i=1,nNode
         do j=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            B11 = nDofN*(j-1)+1
            B12 = nDim*(j-1)+1
            !
            ! displacement
            !
            amatrx(A11,B11)     = Kuu(A12,B12)
            amatrx(A11,B11+1)   = Kuu(A12,B12+1)
            amatrx(A11,B11+2)   = Kuu(A12,B12+2)
            amatrx(A11+1,B11)   = Kuu(A12+1,B12)
            amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
            amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
            amatrx(A11+2,B11)   = Kuu(A12+2,B12)
            amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
            amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
            !
            !  k_uu
            !
            amatrx(A11+3,B11+3) = K_uu(i,j)
            !
            !  k_vv
            !
            amatrx(A11+4,B11+4) = K_vv(i,j)
            !
            !  k_phiphi
            !
            amatrx(A11+5,B11+5) = K_phiphi(i,j)
            !
            !  K_uv
            !
            amatrx(A11+3,B11+4) = K_uv(i,j)
            !
            !  K_uphi
            !
            amatrx(A11+3,B11+5) = K_uphi(i,j)
            !
            !  K_vu
            !
            amatrx(A11+4,B11+3) = K_vu(i,j)
            !
            !  K_vphi
            !
            amatrx(A11+4,B11+5) = K_vphi(i,j)
            !
            !  K_phiu
            !
            amatrx(A11+5,B11+3) = K_phiu(i,j)
            !
            !  K_phiv
            !
            amatrx(A11+5,B11+4) = K_phiv(i,j)
            !
            enddo
         enddo
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------






      return
      end subroutine U3D8

************************************************************************
      subroutine integ(props,nprops,dtime,
     +     u_tau,u_t,v_tau,v_t,
     +     Diff_u,Diff_v,F_const,k_const,sigma,permittivity,
     +     c_t,c_tau)


      implicit none

      integer i,j,k,l,m,n,nprops,stat,nargs
      parameter(nargs=5)

      real*8 props(nprops),dtime,theta_tau,theta_t,theta0,Ktherm,dKdT
      real*8 Cheat,dCdT,rho,heatGen,dHdT
      real*8 krate,Diffv,dDvdv
      real*8 u_tau,u_t,v_tau,v_t
      real*8 mobility,eps,dedphi,dmdmu
      real*8 a0,Stefan,C1,C2
      real*8 temp_tau,temp_t
      real*8 gamma,diffusivity,pa,pb,pe,ph,Tr,au
      real*8 args(nargs),c_t,c_tau
      real*8 Diff_u,Diff_v,F_const,k_const
      real*8 sigma,permittivity

      real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)




      ! Obtain material properties
      !
      Diff_u         = props(1)
      Diff_v         = props(2)
      F_const        = props(3)
      k_const        = props(4)
      ph             = props(5)
      Tr             = props(6)
      gamma          = props(7)
      sigma          = props(8)
      permittivity   = props(9)

c      write(*,*),'Diff_u=',Diff_u
c      write(*,*),'Diff_v=',Diff_v
c      write(*,*),'F_const=',F_const
c      write(*,*),'k_const=',k_const
c      write(*,*),'ph=',ph
c      write(*,*),'Tr=',Tr
c      write(*,*),'gamma=',gamma
c      write(*,*),'sigma=',sigma
c      write(*,*),'permittivity=',permittivity

      ! perform time integration here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! compare u_tau with Tr
      !
      if(v_tau .le. Tr) then
        au = 0.49
      elseif(v_tau .gt. Tr) then
        au = 0.49 - 2.5*(v_tau - Tr)
      endif

      ! call Rtsafe to obtain the current
      args(1) = c_t
      args(2) = dtime
      args(3) = gamma
      args(4) = ph
      args(5) = au

c      c_tau = c_t
      call solvec(c_tau,args,nargs)
      !
      !
      ! perform time integration here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      return
      end subroutine integ
************************************************************************

      subroutine AssembleElement(nDim,nNode,ndofel,
     +     Ru,Rc,Kuu,Kuc,Kcu,Kcc,
     +     rhs,amatrx)

      !
      ! Subroutine to assemble the local elements residual and tangent
      !

      implicit none

      integer i,j,k,l,m,n,A11,A12,B11,B12,nDim,nNode,nDofEl,nDofN

      real*8 Ru(nDim*nNode,1),Rc(nNode,1),Kuu(nDim*nNode,nDim*nNode)
      real*8 Kcc(nNode,nNode),Kuc(nDim*nNode,nNode),rhs(ndofel,1)
      real*8 Kcu(nNode,nDim*nNode),amatrx(ndofel,ndofel)


      ! Total number of degrees of freedom per node
      !
      nDofN = nDofEl/nNode


      ! init
      !
      rhs(:,1) = 0.d0
      amatrx = 0.d0

      if(nDim.eq.2) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1) = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            !
            ! chemical potential
            !
            rhs(A11+2,1) = Rc(i,1)
         enddo
         !
         ! Assemble the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11) = Kuu(A12,B12)
               amatrx(A11,B11+1) = Kuu(A12,B12+1)
               amatrx(A11+1,B11) = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               !
               ! chemical potential
               !
               amatrx(A11+2,B11+2) = Kcc(i,j)
               !
               ! displacement - chemical potential
               !
               amatrx(A11,B11+2) = Kuc(A12,j)
               amatrx(A11+1,B11+2) = Kuc(A12+1,j)
               !
               ! chemical potential - displacement
               !
               amatrx(A11+2,B11) = Kcu(i,B12)
               amatrx(A11+2,B11+1) = Kcu(i,B12+1)
               !
            enddo
         enddo
         !
      elseif(nDim.eq.3) then
         !
         ! Assemble the element level residual
         !
         do i=1,nNode
            A11 = nDofN*(i-1)+1
            A12 = nDim*(i-1)+1
            !
            ! displacement
            !
            rhs(A11,1)   = Ru(A12,1)
            rhs(A11+1,1) = Ru(A12+1,1)
            rhs(A11+2,1) = Ru(A12+2,1)
            !
            ! chemical potential
            !
            rhs(A11+3,1) = Rc(i,1)
            !
         enddo
         !
         ! Assembly the element level tangent matrix
         !
         do i=1,nNode
            do j=1,nNode
               A11 = nDofN*(i-1)+1
               A12 = nDim*(i-1)+1
               B11 = nDofN*(j-1)+1
               B12 = nDim*(j-1)+1
               !
               ! displacement
               !
               amatrx(A11,B11)     = Kuu(A12,B12)
               amatrx(A11,B11+1)   = Kuu(A12,B12+1)
               amatrx(A11,B11+2)   = Kuu(A12,B12+2)
               amatrx(A11+1,B11)   = Kuu(A12+1,B12)
               amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
               amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
               amatrx(A11+2,B11)   = Kuu(A12+2,B12)
               amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
               amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
               !
               ! chemical potential
               !
               amatrx(A11+3,B11+3) = Kcc(i,j)
               !
               ! displacement - chemical potential
               !
               amatrx(A11,B11+3) = Kuc(A12,j)
               amatrx(A11+1,B11+3) = Kuc(A12+1,j)
               amatrx(A11+2,B11+3) = Kuc(A12+2,j)
               !
               ! chemical potential - displacement
               !
               amatrx(A11+3,B11) = Kcu(i,B12)
               amatrx(A11+3,B11+1) = Kcu(i,B12+1)
               amatrx(A11+3,B11+2) = Kcu(i,B12+2)
               !
            enddo
         enddo
         !
      else
         write(*,*) 'How did you get nDim=',nDim
         call xit
      endif

      return
      end subroutine AssembleElement

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0


      ! Gauss pt location in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0


      return
      end subroutine xint2D1pt

!************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0


      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt

************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 2 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer nIntPt,nDim

      real*8 xi(1,3),w(1)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0


      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0

      return
      end subroutine xint3D1pt

************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer nIntPt,nDim

      real*8 xi(8,3),w(8)


      ! Init
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0


      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

************************************************************************

      subroutine xintSurf2D1pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations

      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),w(1),zero,one,two
      parameter(zero=0.d0,one=1.d0,two=2.d0)


      ! Gauss weights
      !
      w(1) = two


      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = zero
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D1pt

************************************************************************

      subroutine xintSurf2D2pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations

      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(2),yLocal(2),w(2),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one


      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = -dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D2pt

************************************************************************

      subroutine xintSurf2D3pt(face,xLocal,yLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 2 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(3),yLocal(3),w(3),zero,one,two,three,five,eight,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,five=5.d0,
     +     eight=8.d0,nine=9.d0)


      ! Gauss weights
      !
      w(1) = five/nine
      w(2) = eight/nine
      w(3) = five/nine


      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = -one
         xLocal(2) = zero
         yLocal(2) = -one
         xLocal(2) = dsqrt(three/five)
         yLocal(2) = -one
      elseif(face.eq.2) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(three/five)
         xLocal(2) = one
         yLocal(2) = zero
         xLocal(3) = one
         yLocal(3) = dsqrt(three/five)
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(three/five)
         yLocal(1) = one
         xLocal(2) = zero
         yLocal(2) = one
         xLocal(3) = dsqrt(three/five)
         yLocal(3) = one
      elseif(face.eq.4) then
         xLocal(1) = -one
         yLocal(1) = dsqrt(three/five)
         xLocal(2) = -one
         yLocal(2) = zero
         xLocal(3) = -one
         yLocal(3) = -dsqrt(three/five)
      else
         write(*,*) 'face.ne.1,2,3,4'
         write(80,*) 'face.ne.1,2,3,4'
         call xit
      endif

      end subroutine xintSurf2D3pt

************************************************************************

      subroutine xintSurf3D1pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 1 gauss point for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  zLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(1),yLocal(1),zLocal(1),w(1),zero,one,four
      parameter(zero=0.d0,one=1.d0,four=4.d0)


      ! Gauss weights
      !
      w(1) = four


      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = -one
      elseif(face.eq.2) then
         xLocal(1) = zero
         yLocal(1) = zero
         zLocal(1) = one
      elseif(face.eq.3) then
         xLocal(1) = zero
         yLocal(1) = -one
         zLocal(1) = zero
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = zero
         zLocal(1) = zero
      elseif(face.eq.5) then
         xLocal(1) = zero
         yLocal(1) = one
         zLocal(1) = zero
      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = zero
         zLocal(1) = zero
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D1pt

************************************************************************

      subroutine xintSurf3D4pt(face,xLocal,yLocal,zLocal,w)

      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 4 gauss points for surface integration
      !
      !  xLocal(nIntPt): x coordinates for the integration pts
      !  yLocal(nIntPt): y coordinates for the integration pts
      !  yLocal(nIntPt): z coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights

      implicit none

      integer face

      real*8 xLocal(4),yLocal(4),zLocal(4),w(4),one,three
      parameter(one=1.d0,three=3.d0)


      ! Gauss weights
      !
      w(1) = one
      w(2) = one
      w(3) = one
      w(4) = one


      ! Gauss pt locations in master element
      !
      if(face.eq.1) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = -one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = -one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = -one
      elseif(face.eq.2) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = one
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -dsqrt(one/three)
         zLocal(2) = one
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = one
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = dsqrt(one/three)
         zLocal(4) = one
      elseif(face.eq.3) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = -one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = -one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = -one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = -one
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.4) then
         xLocal(1) = one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      elseif(face.eq.5) then
         xLocal(1) = -dsqrt(one/three)
         yLocal(1) = one
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = dsqrt(one/three)
         yLocal(2) = one
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = dsqrt(one/three)
         yLocal(3) = one
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -dsqrt(one/three)
         yLocal(4) = one
         zLocal(4) = dsqrt(one/three)

      elseif(face.eq.6) then
         xLocal(1) = -one
         yLocal(1) = -dsqrt(one/three)
         zLocal(1) = -dsqrt(one/three)
         xLocal(2) = -one
         yLocal(2) = dsqrt(one/three)
         zLocal(2) = -dsqrt(one/three)
         xLocal(3) = -one
         yLocal(3) = dsqrt(one/three)
         zLocal(3) = dsqrt(one/three)
         xLocal(4) = -one
         yLocal(4) = -dsqrt(one/three)
         zLocal(4) = dsqrt(one/three)
      else
         write(*,*) 'face.ne.1,2,3,4,5,6'
         write(80,*) 'face.ne.1,2,3,4,5,6'
         call xit
      endif

      end subroutine xintSurf3D4pt

!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element


      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      !                          eta
      !   4-----------3          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          |
      !   |           |          O--------- xi
      !   1-----------2        origin at center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)


      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)


      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)


      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)


      return
      end subroutine calcShape2DLinear

************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! This subroutine uses a 8-node linear 3D element as shown
      !
      !      8-----------7
      !     /|          /|       zeta
      !    / |         / |
      !   5-----------6  |       |     eta
      !   |  |        |  |       |   /
      !   |  |        |  |       |  /
      !   |  4--------|--3       | /
      !   | /         | /        |/
      !   |/          |/         O--------- xi
      !   1-----------2        origin at cube center
      !
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      ! d2shxi(i,j,k) = derivatives wrt j and k of shape fn of node i

      implicit none

      integer intpt,nDim,nIntPt,i,j

      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3)
      real*8 d2shxi(8,3,3),xi,eta,zeta

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      !
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      !
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      !
      ! The second derivatives
      !
      d2shxi = zero
      d2shxi(1,1,2) = eighth*(one - zeta)
      d2shxi(1,2,1) = d2shxi(1,1,2)
      d2shxi(1,1,3) = eighth*(one - eta)
      d2shxi(1,3,1) = d2shxi(1,1,3)
      d2shxi(1,2,3) = eighth*(one - xi)
      d2shxi(1,3,2) = d2shxi(1,2,3)
      d2shxi(2,1,2) = -eighth*(one - zeta)
      d2shxi(2,2,1) = d2shxi(2,1,2)
      d2shxi(2,1,3) = -eighth*(one - eta)
      d2shxi(2,3,1) = d2shxi(2,1,3)
      d2shxi(2,2,3) = eighth*(one + xi)
      d2shxi(2,3,2) = d2shxi(2,2,3)
      d2shxi(3,1,2) = eighth*(one - zeta)
      d2shxi(3,2,1) = d2shxi(2,1,2)
      d2shxi(3,1,3) = -eighth*(one + eta)
      d2shxi(3,3,1) = d2shxi(2,1,3)
      d2shxi(3,2,3) = -eighth*(one + xi)
      d2shxi(3,3,2) = d2shxi(2,2,3)
      d2shxi(4,1,2) = -eighth*(one - zeta)
      d2shxi(4,2,1) = d2shxi(2,1,2)
      d2shxi(4,1,3) = eighth*(one + eta)
      d2shxi(4,3,1) = d2shxi(2,1,3)
      d2shxi(4,2,3) = -eighth*(one - xi)
      d2shxi(4,3,2) = d2shxi(2,2,3)
      d2shxi(5,1,2) = eighth*(one + zeta)
      d2shxi(5,2,1) = d2shxi(2,1,2)
      d2shxi(5,1,3) = -eighth*(one - eta)
      d2shxi(5,3,1) = d2shxi(2,1,3)
      d2shxi(5,2,3) = -eighth*(one - xi)
      d2shxi(5,3,2) = d2shxi(2,2,3)
      d2shxi(6,1,2) = eighth*(one + zeta)
      d2shxi(6,2,1) = d2shxi(2,1,2)
      d2shxi(6,1,3) = eighth*(one - eta)
      d2shxi(6,3,1) = d2shxi(2,1,3)
      d2shxi(6,2,3) = -eighth*(one + xi)
      d2shxi(6,3,2) = d2shxi(2,2,3)
      d2shxi(7,1,2) = eighth*(one + zeta)
      d2shxi(7,2,1) = d2shxi(2,1,2)
      d2shxi(7,1,3) = eighth*(one + eta)
      d2shxi(7,3,1) = d2shxi(2,1,3)
      d2shxi(7,2,3) = eighth*(one + xi)
      d2shxi(7,3,2) = d2shxi(2,2,3)
      d2shxi(8,1,2) = -eighth*(one + zeta)
      d2shxi(8,2,1) = d2shxi(2,1,2)
      d2shxi(8,1,3) = -eighth*(one + eta)
      d2shxi(8,3,1) = d2shxi(2,1,3)
      d2shxi(8,2,3) = eighth*(one - xi)
      d2shxi(8,3,2) = d2shxi(2,2,3)

      return
      end subroutine calcShape3DLinear

!************************************************************************


      subroutine computeSurf(nNode,xLocal,yLocal,face,coords,sh,ds)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the length ds, so that one can
      !  do the numerical integration on the boundary for fluxes
      !  on the 4-node or 8-node quadrilateral elements, and also
      !  the 6-node triangle

      implicit none

      integer face,nNode,k

      real*8 xLocal,yLocal,ds,dshxi(nNode,2),sh(nNode),dXdXi,dXdEta,half
      real*8 dYdEta,one,coords(2,nNode),fourth,shape,normal(2,1),dYdXi
      real*8 lam,two,eight,three,four,zero
      parameter(one=1.d0,fourth=1.d0/4.d0,half=0.5d0,two=2.d0,
     +     eight=8.d0,three=3.d0,four=4.d0,zero=0.d0)

      !
      ! Compute shape functions and derivatives
      !
      if(nNode.eq.4) then
         !
         ! We are dealing with a 4-node quadrilateral
         !
         sh(1) = fourth*(one - xLocal)*(one - yLocal)
         sh(2) = fourth*(one + xLocal)*(one - yLocal)
         sh(3) = fourth*(one + xLocal)*(one + yLocal)
         sh(4) = fourth*(one - xLocal)*(one + yLocal)

         dshxi(1,1) = -fourth*(one - yLocal)
         dshxi(1,2) = -fourth*(one - xLocal)
         dshxi(2,1) = fourth*(one - yLocal)
         dshxi(2,2) = -fourth*(one + xLocal)
         dshxi(3,1) = fourth*(one + yLocal)
         dshxi(3,2) = fourth*(one + xLocal)
         dshxi(4,1) = -fourth*(one + yLocal)
         dshxi(4,2) = fourth*(one - xLocal)

         dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +        + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
         dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +        + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
         dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +        + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
         dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +        + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)
         !
      elseif(nNode.eq.8) then
         !
         ! We are dealing with an 8-node quad
         !
         sh(1) = -fourth*(one-xLocal)*(one-yLocal)*(one+xLocal+yLocal)
         sh(2) = -fourth*(one+xLocal)*(one-yLocal)*(one-xLocal+yLocal)
         sh(3) = -fourth*(one+xLocal)*(one+yLocal)*(one-xLocal-yLocal)
         sh(4) = -fourth*(one-xLocal)*(one+yLocal)*(one+xLocal-yLocal)
         sh(5) = half*(one-xLocal)*(one+xLocal)*(one-yLocal)
         sh(6) = half*(one-yLocal)*(one+yLocal)*(one+xLocal)
         sh(7) = half*(one-xLocal)*(one+xLocal)*(one+yLocal)
         sh(8) = half*(one-yLocal)*(one+yLocal)*(one-xLocal)
         !
         dshxi(1,1) =  fourth*(one - yLocal)*(2.d0*xLocal + yLocal)
         dshxi(1,2) =  fourth*(one - xLocal)*(xLocal + 2.d0*yLocal)
         dshxi(2,1) =  fourth*(2.d0*xLocal - yLocal)*(one - yLocal)
         dshxi(2,2) = -fourth*(one + xLocal)*(xLocal - 2.d0*yLocal)
         dshxi(3,1) =  fourth*(one + yLocal)*(2.d0*xLocal + yLocal)
         dshxi(3,2) =  fourth*(one + xLocal)*(xLocal + 2.d0*yLocal)
         dshxi(4,1) =  fourth*(2.d0*xLocal - yLocal)*(one + yLocal)
         dshxi(4,2) = -fourth*(one - xLocal)*(xLocal - 2.d0*yLocal)
         dshxi(5,1) = -xLocal*(one - yLocal)
         dshxi(5,2) = -half*(one - xLocal*xLocal)
         dshxi(6,1) =  half*(one - yLocal*yLocal)
         dshxi(6,2) = -yLocal*(one + xLocal)
         dshxi(7,1) = -xLocal*(one + yLocal)
         dshxi(7,2) =  half*(one - xLocal*xLocal)
         dshxi(8,1) = -half*(one - yLocal*yLocal)
         dshxi(8,2) = -yLocal*(one - xLocal)
         !
         dXdXi = dshxi(1,1)*coords(1,1)+dshxi(2,1)*coords(1,2)
     +        + dshxi(3,1)*coords(1,3)+dshxi(4,1)*coords(1,4)
     +        + dshxi(5,1)*coords(1,5)+dshxi(6,1)*coords(1,6)
     +        + dshxi(7,1)*coords(1,7)+dshxi(8,1)*coords(1,8)
         dXdEta = dshxi(1,2)*coords(1,1)+dshxi(2,2)*coords(1,2)
     +        + dshxi(3,2)*coords(1,3)+dshxi(4,2)*coords(1,4)
     +        + dshxi(5,2)*coords(1,5)+dshxi(6,2)*coords(1,6)
     +        + dshxi(7,2)*coords(1,7)+dshxi(8,2)*coords(1,8)
         dYdXi = dshxi(1,1)*coords(2,1)+dshxi(2,1)*coords(2,2)
     +        + dshxi(3,1)*coords(2,3)+dshxi(4,1)*coords(2,4)
     +        + dshxi(5,1)*coords(2,5)+dshxi(6,1)*coords(2,6)
     +        + dshxi(7,1)*coords(2,7)+dshxi(8,1)*coords(2,8)
         dYdEta = dshxi(1,2)*coords(2,1)+dshxi(2,2)*coords(2,2)
     +        + dshxi(3,2)*coords(2,3)+dshxi(4,2)*coords(2,4)
     +        + dshxi(5,2)*coords(2,5)+dshxi(6,2)*coords(2,6)
     +        + dshxi(7,2)*coords(2,7)+dshxi(8,2)*coords(2,8)
         !
      elseif(nNode.eq.6) then
         !
         ! We are dealing with a 6-node triangle
         !
         ! lam is a parametric coordinate going from
         !  node 2 to node 3, and useful in what follows
         !
         lam = one - xLocal - yLocal
         !
         ! The shape functions
         !
         sh(1) = lam*(two*lam - one)
         sh(2) = xLocal*(two*xLocal - one)
         sh(3) = yLocal*(two*yLocal - one)
         sh(4) = four*xLocal*lam
         sh(5) = four*xLocal*yLocal
         sh(6) = four*yLocal*lam
         !
         ! The first derivatives
         !
         dshxi(1,1) = one - four*lam
         dshxi(1,2) = one - four*lam
         dshxi(2,1) = four*xLocal - one
         dshxi(2,2) = zero
         dshxi(3,1) = zero
         dshxi(3,2) = four*yLocal - one
         dshxi(4,1) = four*(lam - xLocal)
         dshxi(4,2) = -four*xLocal
         dshxi(5,1) = four*yLocal
         dshxi(5,2) = four*xLocal
         dshxi(6,1) = -four*yLocal
         dshxi(6,2) = four*(lam - yLocal)
         !
         dXdXi = zero
         dXdEta = zero
         dYdXi = zero
         dYdEta = zero
         do k=1,6
            dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
            dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
            dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
            dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
         enddo
         !
      else
         write(*,*) 'In computeSurf, nNode=',nNode
         write(80,*) 'In computeSurf, nNode=',nNode
         call xit
      endif


      ! Compute the jacobian of the mapping and the normal
      !
      if((nNode.eq.4).or.(nNode.eq.8)) then
         !
         ! We have either a 4-node or 8-node quadrilateral
         !
         ! Jacobian of the mapping
         !
         if((face.eq.2).or.(face.eq.4)) then
            ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         elseif((face.eq.1).or.(face.eq.3)) then
            ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         else
            write(*,*) 'never should get here'
            write(80,*) 'never should get here'
            call xit
         endif


         ! Surface normal, outward pointing in this case. Useful for
         !  ``follower'' type loads. The normal is referential or spatial
         !  depending on which coords were supplied to this subroutine
         !  (NOT fully tested)
         !
         if((face.eq.2).or.(face.eq.4)) then
            normal(1,1) = dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
            normal(2,1) = -dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
            if(face.eq.4) normal = -normal
         elseif((face.eq.1).or.(face.eq.3)) then
            normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
            normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
            if(face.eq.3) normal = -normal
         else
            write(*,*) 'never should get here'
            write(80,*) 'never should get here'
            call xit
         endif
         !
      elseif(nNode.eq.6) then
         !
         ! We have a 6-node triangle
         !
         ! Jacobian of the mapping (factors of one half are
         !  due to the fact that this master element spans
         !  from [0 1] and not [-1 1])
         !
         if(face.eq.1) then
            ! eta=constant on this face
            ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)/two
         elseif(face.eq.2) then
            ! eta = 1-xi on this face, needs special attention
            dXdXi = (four*xLocal - one)*coords(1,2)
     +           + (four*xLocal - three)*coords(1,3)
     +           + (four - eight*xLocal)*coords(1,5)
            dYdXi = (four*xLocal - one)*coords(2,2)
     +           + (four*xLocal - three)*coords(2,3)
     +           + (four - eight*xLocal)*coords(2,5)
            ds = dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)/two
         elseif(face.eq.3) then
            ! xi=constant on this face
            ds = dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)/two
         else
            write(*,*) 'never should get here'
            write(80,*) 'never should get here'
            call xit
         endif


         ! Surface normal, outward pointing in this case. Useful for
         !  ``follower'' type loads. The normal is referential or spatial
         !  depending on which coords were supplied to this subroutine
         !  (NOT TESTED!!!!!!)
         !
         if(face.eq.1) then
            normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
            normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         elseif(face.eq.2) then
            ! recall that eta=1-xi on this face
            normal(1,1) = dYdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
            normal(2,1) = -dXdXi/dsqrt(dXdXi*dXdXi + dYdXi*dYdXi)
         elseif(face.eq.3) then
            normal(1,1) = -dYdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
            normal(2,1) = dXdEta/dsqrt(dXdEta*dXdEta + dYdEta*dYdEta)
         else
         endif
         !
      else
         write(*,*) 'In computeSurf, nNode=',nNode
         write(80,*) 'In computeSurf, nNode=',nNode
         call xit
      endif

      return
      end subroutine computeSurf

************************************************************************

      subroutine computeSurf3D(xLocal,yLocal,zLocal,face,coords,sh,dA)

      ! This subroutine computes the shape functions, derivatives
      !  of shape functions, and the area dA, so that one can
      !  do the numerical integration on the boundary for fluxes
      !  on the 8-node brick elements

      implicit none

      integer face,stat,i,j,k

      real*8 xLocal,yLocal,zLocal,dA,dshxi(8,3),sh(8),zero,dsh(8,3),one
      real*8 coords(3,8),two,eighth,mapJ(3,3),mag,normal(3,1)

      real*8 dXdXi,dXdEta,dXdZeta,dYdXi,dYdEta,dYdZeta,dZdXi,dZdEta
      real*8 dZdZeta

      parameter(one=1.d0,two=2.d0,eighth=1.d0/8.d0,zero=0.d0)

      ! The shape functions
      !
      sh(1) = eighth*(one - xLocal)*(one - yLocal)*(one - zLocal)
      sh(2) = eighth*(one + xLocal)*(one - yLocal)*(one - zLocal)
      sh(3) = eighth*(one + xLocal)*(one + yLocal)*(one - zLocal)
      sh(4) = eighth*(one - xLocal)*(one + yLocal)*(one - zLocal)
      sh(5) = eighth*(one - xLocal)*(one - yLocal)*(one + zLocal)
      sh(6) = eighth*(one + xLocal)*(one - yLocal)*(one + zLocal)
      sh(7) = eighth*(one + xLocal)*(one + yLocal)*(one + zLocal)
      sh(8) = eighth*(one - xLocal)*(one + yLocal)*(one + zLocal)


      ! Shape function derivatives
      !
      dshxi(1,1) = -eighth*(one - yLocal)*(one - zLocal)
      dshxi(1,2) = -eighth*(one - xLocal)*(one - zLocal)
      dshxi(1,3) = -eighth*(one - xLocal)*(one - yLocal)
      dshxi(2,1) = eighth*(one - yLocal)*(one - zLocal)
      dshxi(2,2) = -eighth*(one + xLocal)*(one - zLocal)
      dshxi(2,3) = -eighth*(one + xLocal)*(one - yLocal)
      dshxi(3,1) = eighth*(one + yLocal)*(one - zLocal)
      dshxi(3,2) = eighth*(one + xLocal)*(one - zLocal)
      dshxi(3,3) = -eighth*(one + xLocal)*(one + yLocal)
      dshxi(4,1) = -eighth*(one + yLocal)*(one - zLocal)
      dshxi(4,2) = eighth*(one - xLocal)*(one - zLocal)
      dshxi(4,3) = -eighth*(one - xLocal)*(one + yLocal)
      dshxi(5,1) = -eighth*(one - yLocal)*(one + zLocal)
      dshxi(5,2) = -eighth*(one - xLocal)*(one + zLocal)
      dshxi(5,3) = eighth*(one - xLocal)*(one - yLocal)
      dshxi(6,1) = eighth*(one - yLocal)*(one + zLocal)
      dshxi(6,2) = -eighth*(one + xLocal)*(one + zLocal)
      dshxi(6,3) = eighth*(one + xLocal)*(one - yLocal)
      dshxi(7,1) = eighth*(one + yLocal)*(one + zLocal)
      dshxi(7,2) = eighth*(one + xLocal)*(one + zLocal)
      dshxi(7,3) = eighth*(one + xLocal)*(one + yLocal)
      dshxi(8,1) = -eighth*(one + yLocal)*(one + zLocal)
      dshxi(8,2) = eighth*(one - xLocal)*(one + zLocal)
      dshxi(8,3) = eighth*(one - xLocal)*(one + yLocal)


      dXdXi = zero
      dXdEta = zero
      dXdZeta = zero
      dYdXi = zero
      dYdEta = zero
      dYdZeta = zero
      dZdXi = zero
      dZdEta = zero
      dZdZeta = zero
      do k=1,8
         dXdXi = dXdXi + dshxi(k,1)*coords(1,k)
         dXdEta = dXdEta + dshxi(k,2)*coords(1,k)
         dXdZeta = dXdZeta + dshxi(k,3)*coords(1,k)
         dYdXi = dYdXi + dshxi(k,1)*coords(2,k)
         dYdEta = dYdEta + dshxi(k,2)*coords(2,k)
         dYdZeta = dYdZeta + dshxi(k,3)*coords(2,k)
         dZdXi = dZdXi + dshxi(k,1)*coords(3,k)
         dZdEta = dZdEta + dshxi(k,2)*coords(3,k)
         dZdZeta = dZdZeta + dshxi(k,3)*coords(3,k)
      enddo


      ! Jacobian of the mapping
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdEta - dYdEta*dZdXi)**two
     +        + (dXdXi*dZdEta - dXdEta*dZdXi)**two
     +        + (dXdXi*dYdEta - dXdEta*dYdXi)**two
     +        )
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         dA = dsqrt(
     +          (dYdXi*dZdZeta - dYdZeta*dZdXi)**two
     +        + (dXdXi*dZdZeta - dXdZeta*dZdXi)**two
     +        + (dXdXi*dYdZeta - dXdZeta*dYdXi)**two
     +        )
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         dA = dsqrt(
     +          (dYdEta*dZdZeta - dYdZeta*dZdEta)**two
     +        + (dXdEta*dZdZeta - dXdZeta*dZdEta)**two
     +        + (dXdEta*dYdZeta - dXdZeta*dYdEta)**two
     +        )
         else
            write(*,*) 'never should get here'
            call xit
      endif


      ! Surface normal, outward pointing in this case. Useful for
      !  ``follower'' type loads. The normal is referential or spatial
      !  depending on which coords were supplied to this subroutine
      !  (NOT fully tested)
      !
      if((face.eq.1).or.(face.eq.2)) then
         ! zeta = constant on this face
         normal(1,1) = dYdXi*dZdEta - dYdEta*dZdXi
         normal(2,1) = dXdXi*dZdEta - dXdEta*dZdXi
         normal(3,1) = dXdXi*dYdEta - dXdEta*dYdXi
         if(face.eq.1) normal = -normal
      elseif((face.eq.3).or.(face.eq.5)) then
         ! eta = constant on this face
         normal(1,1) = dYdXi*dZdZeta - dYdZeta*dZdXi
         normal(2,1) = dXdXi*dZdZeta - dXdZeta*dZdXi
         normal(3,1) = dXdXi*dYdZeta - dXdZeta*dYdXi
         if(face.eq.5) normal = -normal
      elseif((face.eq.4).or.(face.eq.6)) then
         ! xi = constant on this face
         normal(1,1) = dYdEta*dZdZeta - dYdZeta*dZdEta
         normal(2,1) = dXdEta*dZdZeta - dXdZeta*dZdEta
         normal(3,1) = dXdEta*dYdZeta - dXdZeta*dYdEta
         if(face.eq.6) normal = -normal
      else
         write(*,*) 'never should get here'
         call xit
      endif
      mag = dsqrt(normal(1,1)**two+normal(2,1)**two+normal(3,1)**two)
      normal(1,1) = normal(1,1)/mag
      normal(2,1) = normal(2,1)/mag
      normal(3,1) = normal(3,1)/mag

      end subroutine computeSurf3D

************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      return
      end subroutine mapShape2Da

************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.  This subroutine works for both 8-node
      !  linear and 20-node quadratic 3D elements.
      !
      implicit none

      integer i,j,k,nNode,ieror,stat

      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode)
      real*8 mapJ(3,3),mapJ_inv(3,3),detmapJ

      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      ! The second derivatives may be calculated.
      !

      return
      end subroutine mapShape3D

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1

      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))

      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv3D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if

      det_A_inv = 1.d0/det_A

      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))


      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv


      istat = 1

      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)

      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if

      det_A_inv = 1.d0/det_A

      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3)
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet

!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

****************************************************************************
      subroutine solvec(root,args,nargs)

      ! This subroutine will numerically solve for c_tau

      implicit none

      integer maxit,j,nargs

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin
      real*8 LambdaBar_tau,LambdaBar_t

      parameter(maxit=100)
      parameter(xacc=1.d-5,zero=0.d0,one=1.d0)
!      parameter(zero=0.d0,one=1.d0)




      rootMax =  1.d0
      rootMin =  0.d0

      x1 = rootMin
      x2 = rootMax
      call gfunc(x1,FL,DF,args,nargs)
      call gfunc(x2,FH,DF,args,nargs)

      if(fl*fh.ge.zero) then
         write(*,*) 'FYI, root not bracketed on mu'
         write(*,*) 'fl=',fl
         write(*,*) 'x1=',x1
         write(*,*) 'fh=',fh
         write(*,*) 'x2=',x2
         call xit
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      ROOT = 0.5D0 *( X1 + X2)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD

      call gfunc(root,F,DF,args,nargs)

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
         IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call gfunc(root,F,DF,args,nargs)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'solveRate EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solveRate EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solvec
****************************************************************************
      subroutine gfunc(root,f,df,args,nargs)

      implicit none

      integer nargs

      real*8 args(nargs),mu,f,df,mu_t,Amp
      real*8 mu_ss,lamL,LambdaBar_tau,LambdaBar_t,Achainmax
      real*8 X_tau,X_t,Achain_tau,Achain_t,aux1,aux2
      real*8 c_t,dtime,gamma,ph,au,root


      real*8 zero,one,two,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)




      ! Obtain relevant quantities
      !
      c_t    = args(1)
      dtime  = args(2)
      gamma  = args(3)
      ph     = args(4)
      au     = args(5)



      ! Compute the residual
      !
      f = c_t - root + dtime*gamma*ph*(au*(root**2.0)
     +    -au*root-(root**3.0) + (root**2.0))




      ! Compute the tangent
      !
      df = -1.0 + dtime*gamma*ph*(2.0*au*root -au -3.0*root**2.0 + 2.0*root)


      return
      end subroutine gfunc
****************************************************************************
