! Routine dvmlm, a limited memory quasy-Netwon method, 
! by B.M. Averick, R.G. Carter, J.J. More
! obtained from http://ftp.eq.uc.pt/software/math/minpack-2/vmlm/dvmlm.f

! *****************************************************************
! 
!            COPYRIGHT NOTIFICATION
! 
! This program discloses material protectable under copyright laws of
! the United States. Permission to copy and modify this software and its
! documentation for internal research use is hereby granted, provided
! that this notice is retained thereon and on all copies or modifications. 
! The University of Chicago makes no representations as to the suitability 
! and operability of this software for any purpose. 
! It is provided "as is" without express or implied warranty.
! 
! Use of this software for commercial purposes is expressly prohibited
! without contacting 
!
!    Jorge J. More'
!    Mathematics and Computer Science Division
!    Argonne National Laboratory
!    9700 S. Cass Ave.
!    Argonne, Illinois 60439-4844
!    e-mail: more@mcs.anl.gov
!
! Argonne National Laboratory with facilities in the states of
! Illinois and Idaho, is owned by The United States Government, and
! operated by the University of Chicago under provision of a contract
! with the Department of Energy.
!
! *****************************************************************
! 
!            ADDITIONAL INFORMATION
! 
! D. C. Liu and J. Nocedal, 
! On the limited memory BFGS method for large scale optimization, 
! Math. Programming, 45 (1989), pp. 503--528.
! 
! J. Nocedal, 
! The performance of several algorithms for large-scale 
! unconstrained optimization, 
! in Large-Scale Numerical Optimization, 
! T. F.  Coleman and Y. Li, eds., 
! Society for Industrial and Applied Mathematics, 1991, pp. 138--151.
! 
! B. M. Averick and J. J. More',
! Evaluation of large-scale optimization problems on vector 
! and parallel architectures,
! SIAM J. Optimization 4, (1994), 708-721.
! 
! *****************************************************************


      subroutine dvmlm(n,x,f,fgrad,frtol,fatol,fmin,task,m,s,y,rho,
     +                 isave,dsave,wa1,wa2)
      character*(*) task
      integer n, m
      integer isave(5)
      double precision f, frtol, fatol, fmin
      double precision x(n), fgrad(n), s(n,m), y(n,m), rho(m), wa1(n),
     +                 wa2(m), dsave(24)
c     **********
c
c     Subroutine dvmlm
c
c     This subroutine computes a local minimizer of a function
c     of n variables by a limited memory variable metric method.
c     The user must evaluate the function and the gradient.
c
c     This subroutine uses reverse communication.
c     The user must choose an initial approximation x to the
c     minimizer, evaluate the function and the gradient at x,
c     and make the initial call with task set to 'START'.
c     On exit task indicates the required action.
c
c     A typical invocation of dvmlm has the following outline:
c
c     Choose a starting vector x.
c     Evaluate the function at x; store in f.
c     Evaluate the gradient at x; store in fgrad.
c
c     task = 'START'
c  10 continue
c        call dvmlm(n,x,f,fgrad,frtol,fatol,fmin,task,m,s,y,rho,
c                   isave,dsave,wa1,wa2)
c        if (task .eq. 'FG') then
c           Evaluate the function at x; store in f.
c           Evaluate the gradient at x; store in fgrad.
c           go to 10
c        else if (task .eq. 'NEWX') then
c           The approximation x, function f, and gradient fgrad
c           are available for inspection.
c           go to 10
c        end if
c
c     NOTE: The user must not alter work arrays between calls.
c
c     The subroutine statement is
c
c       subroutine dvmlm(n,x,f,fgrad,frtol,fatol,fmin,task,m,s,y,rho,
c                        isave,dsave,wa1,wa2)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x is an approximation to the solution.
c         On exit x is the current approximation.
c
c       f is a double precision variable.
c         On entry f is the value of the function at x.
c         On final exit f is the value of the function at x.
c
c       fgrad is a double precision array of dimension n.
c         On entry fgrad is the value of the gradient at x.
c         On final exit fgrad is the value of the gradient at x.
c
c       frtol is a double precision variable.
c         On entry frtol specifies the relative error desired in the
c            function. Convergence occurs if the estimate of the
c            relative error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than frtol.
c         On exit frtol is unchanged.
c
c       fatol is a double precision variable.
c         On entry fatol specifies the absolute error desired in the
c            function. Convergence occurs if the estimate of the
c            absolute error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than fatol.
c         On exit fatol is unchanged.
c
c       fmin is a double precision variable.
c         On entry fmin specifies a lower bound for the function.
c            The subroutine exits with a warning if f < fmin.
c         On exit fmin is unchanged.
c
c       task is a character variable of length at least 60.
c         On initial entry task must be set to 'START'.
c         On exit task indicates the required action:
c
c            If task(1:2) = 'FG' then evaluate the function and
c            gradient at x and call dvmlm again.
c
c            If task(1:4) = 'NEWX' then a new iterate has been
c            computed. The approximation x, function f, and
c            gradient fgrad are available for examination.
c
c            If task(1:4) = 'CONV' then the search is successful.
c
c            If task(1:4) = 'WARN' then the subroutine is not able
c            to satisfy the convergence conditions. The exit value
c            of x contains the best approximation found.
c
c            If task(1:5) = 'ERROR' then there is an error in the
c            input arguments.
c
c         On exit with convergence, a warning or an error, the
c            variable task contains additional information.
c
c       m is an integer variable.
c          On entry m specifies the amount of storage.
c          On exit m is unchanged.
c
c       s is a double precision work array of dimension (n,m).
c
c       y is a double precision work array of dimension (n,m).
c
c       rho is a double precision work array of dimension m.
c
c       isave is an integer work array of dimension 5.
c
c       dsave is a double precision work array of dimension 24.
c
c       wa1 is a double precision work array of dimension n.
c
c       wa2 is a double precision work array of dimension m.
c
c     Subprograms called
c
c       MINPACK-2 ... dcsrch, dlmmv
c
c       Level 1 BLAS ... daxpy, dcopy, ddot, dnrm2, dscal
c
c     MINPACK-2 Project. April 1995.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
c
c     **********
      double precision zero, one
      parameter (zero=0.0d0,one=1.0d0)

      character*30 work
      integer i, iter, mark
      double precision gd, gd0, f0, scale, sftol, sgtol, stp, stpmin,
     +                 stpmax, sxtol

      double precision dnrm2, ddot
      external dcsrch, dlmmv, daxpy, dcopy, ddot, dnrm2, dscal

      if (task(1:5) .eq. 'START') then

c        Check the input arguments for errors.

         if (n .le. 0) task = 'ERROR: N .LE. 0'
         if (m .le. 0) task = 'ERROR: M .LE. 0'
         if (frtol .le. zero) task = 'ERROR: FRTOL .LE. 0'
         if (fatol .le. zero) task = 'ERROR: FATOL .LE. 0'
         if (f .le. fmin) task = 'ERROR: INITIAL F .LE. FMIN'

c        Exit if there are errors on input.

         if (task(1:5) .eq. 'ERROR') return

c        Initialize local variables.

         iter = 1
         mark = 1

c        Initialize step information.

         scale = dnrm2(n,fgrad,1)
         do 10 i = 1, n
            s(i,1) = fgrad(i)/scale
   10    continue

c        Initialize line search parameters.

         sftol = 1.0d-3
         sgtol = 0.9d0
         sxtol = 0.1d0

c        Set work to start the search.

         work = 'START SEARCH'

      else

c        Restore local variables.

         if (isave(1) .eq. 1) work = 'SEARCH'
         if (isave(1) .eq. 2) work = 'SEARCH DIRECTION'
         iter = isave(2)
         mark = isave(3)
         sftol = dsave(1)
         sgtol = dsave(2)
         sxtol = dsave(3)
         f0 = dsave(4)
         gd = dsave(5)
         gd0 = dsave(6)
         stp = dsave(7)
         stpmin = dsave(8)
         stpmax = dsave(9)
         scale = dsave(10)
      end if

   20 continue

      if (work .eq. 'START SEARCH') then

c        Initialize the line search subroutine.

         f0 = f
         stp = one
         gd0 = -ddot(n,fgrad,1,s(1,mark),1)
         stpmin = zero
         stpmax = (fmin-f0)/(sgtol*gd0)
         stp = min(stp,stpmax)
         call dcopy(n,x,1,wa1,1)
         call dcopy(n,fgrad,1,y(1,mark),1)
         task = 'START SEARCH'
         work = 'SEARCH'

      end if

      if (work .eq. 'SEARCH') then

c        Determine the line search parameter.

         if (f .lt. fmin) then
            task = 'WARNING: F .LT. FMIN'
            go to 30
         end if
         gd = -ddot(n,fgrad,1,s(1,mark),1)

         call dcsrch(stp,f,gd,sftol,sgtol,sxtol,task,stpmin,stpmax,
     +               isave(4),dsave(11))

c        Compute the new iterate.

         call dcopy(n,wa1,1,x,1)
         call daxpy(n,-stp,s(1,mark),1,x,1)

c        Continue if the line search has converged.

         if (task(1:4) .ne. 'CONV' .and.
     +       task .ne. 'WARNING: XTOL TEST SATISFIED') go to 30

c        Compute the step and gradient change.

         iter = iter + 1
         call daxpy(n,-one,fgrad,1,y(1,mark),1)
         call dscal(n,stp,s(1,mark),1)
         rho(mark) = ddot(n,y(1,mark),1,s(1,mark),1)

c        Compute the scale.

         if (rho(mark) .gt. zero) then
            scale = rho(mark)/ddot(n,y(1,mark),1,y(1,mark),1)
         else
            scale = one
         end if

c        Set task to signal a new iterate.
c        Set work to compute a new search direction.

         task = 'NEWX'
         work = 'SEARCH DIRECTION'

c        Test for convergence.

         if (abs(f-f0) .le. fatol .and. stp*abs(gd0) .le. fatol)
     +       task = 'CONVERGENCE: FATOL TEST SATISFIED'
         if (abs(f-f0) .le. frtol*abs(f0) .and.
     +       stp*abs(gd0) .le. frtol*abs(f0))
     +       task = 'CONVERGENCE: FRTOL TEST SATISFIED'

         go to 30

      end if

      if (work .eq. 'SEARCH DIRECTION') then

c        Compute -H*g.

         call dcopy(n,fgrad,1,wa1,1)

         call dlmmv(n,min(m,iter-1),s,y,rho,scale,mark,wa1,wa2)

         mark = mark + 1
         if (mark .eq. m+1) mark = 1
         call dcopy(n,wa1,1,s(1,mark),1)

c        Set task and work to initialize the line search.

         task = 'START SEARCH'
         work = 'START SEARCH'

         go to 20

      end if

   30 continue

c     Save local variables.

      if (work .eq. 'SEARCH') isave(1) = 1
      if (work .eq. 'SEARCH DIRECTION') isave(1) = 2
      isave(2) = iter
      isave(3) = mark
      dsave(1) = sftol
      dsave(2) = sgtol
      dsave(3) = sxtol
      dsave(4) = f0
      dsave(5) = gd
      dsave(6) = gd0
      dsave(7) = stp
      dsave(8) = stpmin
      dsave(9) = stpmax
      dsave(10) = scale

      end


      subroutine dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,
     +                  isave,dsave)
      character*(*) task
      integer isave(2)
      double precision f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      double precision dsave(13)
c     **********
c
c     Subroutine dcsrch
c
c     This subroutine finds a step that satisfies a sufficient
c     decrease condition and a curvature condition.
c
c     Each call of the subroutine updates an interval with
c     endpoints stx and sty. The interval is initially chosen
c     so that it contains a minimizer of the modified function
c
c           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
c
c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c     interval is chosen so that it contains a minimizer of f.
c
c     The algorithm is designed to find a step that satisfies
c     the sufficient decrease condition
c
c           f(stp) <= f(0) + ftol*stp*f'(0),
c
c     and the curvature condition
c
c           abs(f'(stp)) <= gtol*abs(f'(0)).
c
c     If ftol is less than gtol and if, for example, the function
c     is bounded below, then there is always a step which satisfies
c     both conditions.
c
c     If no step can be found that satisfies both conditions, then
c     the algorithm stops with a warning. In this case stp only
c     satisfies the sufficient decrease condition.
c
c     A typical invocation of dcsrch has the following outline:
c
c     Evaluate the function at stp = 0.0d0; store in f.
c     Evaluate the gradient at stp = 0.0d0; store in g.
c     Choose a starting step stp.
c
c     task = 'START'
c  10 continue
c        call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,
c    +               isave,dsave)
c        if (task .eq. 'FG') then
c           Evaluate the function and the gradient at stp
c           go to 10
c           end if
c
c     NOTE: The user must not alter work arrays between calls.
c
c     The subroutine statement is
c
c       subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
c                         task,isave,dsave)
c     where
c
c       stp is a double precision variable.
c         On entry stp is the current estimate of a satisfactory
c            step. On initial entry, a positive initial estimate
c            must be provided.
c         On exit stp is the current estimate of a satisfactory step
c            if task = 'FG'. If task = 'CONV' then stp satisfies
c            the sufficient decrease and curvature condition.
c
c       f is a double precision variable.
c         On initial entry f is the value of the function at 0.
c            On subsequent entries f is the value of the
c            function at stp.
c         On exit f is the value of the function at stp.
c
c       g is a double precision variable.
c         On initial entry g is the derivative of the function at 0.
c            On subsequent entries g is the derivative of the
c            function at stp.
c         On exit g is the derivative of the function at stp.
c
c       ftol is a double precision variable.
c         On entry ftol specifies a nonnegative tolerance for the
c            sufficient decrease condition.
c         On exit ftol is unchanged.
c
c       gtol is a double precision variable.
c         On entry gtol specifies a nonnegative tolerance for the
c            curvature condition.
c         On exit gtol is unchanged.
c
c       xtol is a double precision variable.
c         On entry xtol specifies a nonnegative relative tolerance
c            for an acceptable step. The subroutine exits with a
c            warning if the relative difference between sty and stx
c            is less than xtol.
c         On exit xtol is unchanged.
c
c       task is a character variable of length at least 60.
c         On initial entry task must be set to 'START'.
c         On exit task indicates the required action:
c
c            If task(1:2) = 'FG' then evaluate the function and
c            derivative at stp and call dcsrch again.
c
c            If task(1:4) = 'CONV' then the search is successful.
c
c            If task(1:4) = 'WARN' then the subroutine is not able
c            to satisfy the convergence conditions. The exit value of
c            stp contains the best point found during the search.
c
c            If task(1:5) = 'ERROR' then there is an error in the
c            input arguments.
c
c         On exit with convergence, a warning or an error, the
c            variable task contains additional information.
c
c       stpmin is a double precision variable.
c         On entry stpmin is a nonnegative lower bound for the step.
c         On exit stpmin is unchanged.
c
c       stpmax is a double precision variable.
c         On entry stpmax is a nonnegative upper bound for the step.
c         On exit stpmax is unchanged.
c
c       isave is an integer work array of dimension 2.
c
c       dsave is a double precision work array of dimension 13.
c
c     Subprograms called
c
c       MINPACK-2 ... dcstep
c
c     MINPACK-1 Project. June 1983.
c     Argonne National Laboratory.
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
c
c     **********
      double precision zero, p5, p66
      parameter (zero=0.0d0,p5=0.5d0,p66=0.66d0)
      double precision xtrapl, xtrapu
      parameter (xtrapl=1.1d0,xtrapu=4.0d0)

      logical brackt
      integer stage
      double precision finit, ftest, fm, fx, fxm, fy, fym, ginit, gtest,
     +                 gm, gx, gxm, gy, gym, stx, sty, stmin, stmax,
     +                 width, width1

      external dcstep

c     Initialization block.

      if (task(1:5) .eq. 'START') then

c        Check the input arguments for errors.

         if (stp .lt. stpmin) task = 'ERROR: STP .LT. STPMIN'
         if (stp .gt. stpmax) task = 'ERROR: STP .GT. STPMAX'
         if (g .ge. zero) task = 'ERROR: INITIAL G .GE. ZERO'
         if (ftol .lt. zero) task = 'ERROR: FTOL .LT. ZERO'
         if (gtol .lt. zero) task = 'ERROR: GTOL .LT. ZERO'
         if (xtol .lt. zero) task = 'ERROR: XTOL .LT. ZERO'
         if (stpmin .lt. zero) task = 'ERROR: STPMIN .LT. ZERO'
         if (stpmax .lt. stpmin) task = 'ERROR: STPMAX .LT. STPMIN'

c        Exit if there are errors on input.

         if (task(1:5) .eq. 'ERROR') return

c        Initialize local variables.

         brackt = .false.
         stage = 1
         finit = f
         ginit = g
         gtest = ftol*ginit
         width = stpmax - stpmin
         width1 = width/p5

c        The variables stx, fx, gx contain the values of the step,
c        function, and derivative at the best step.
c        The variables sty, fy, gy contain the value of the step,
c        function, and derivative at sty.
c        The variables stp, f, g contain the values of the step,
c        function, and derivative at stp.

         stx = zero
         fx = finit
         gx = ginit
         sty = zero
         fy = finit
         gy = ginit
         stmin = zero
         stmax = stp + xtrapu*stp
         task = 'FG'

         go to 10

      else

c        Restore local variables.

         if (isave(1) .eq. 1) then
            brackt = .true.
         else
            brackt = .false.
         end if
         stage = isave(2)
         ginit = dsave(1)
         gtest = dsave(2)
         gx = dsave(3)
         gy = dsave(4)
         finit = dsave(5)
         fx = dsave(6)
         fy = dsave(7)
         stx = dsave(8)
         sty = dsave(9)
         stmin = dsave(10)
         stmax = dsave(11)
         width = dsave(12)
         width1 = dsave(13)

      end if

c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c     algorithm enters the second stage.

      ftest = finit + stp*gtest
      if (stage .eq. 1 .and. f .le. ftest .and. g .ge. zero) stage = 2

c     Test for warnings.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     +    task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      if (brackt .and. stmax-stmin .le. xtol*stmax)
     +    task = 'WARNING: XTOL TEST SATISFIED'
      if (stp .eq. stpmax .and. f .le. ftest .and. g .le. gtest)
     +    task = 'WARNING: STP = STPMAX'
      if (stp .eq. stpmin .and. (f .gt. ftest .or. g .ge. gtest))
     +    task = 'WARNING: STP = STPMIN'

c     Test for convergence.

      if (f .le. ftest .and. abs(g) .le. gtol*(-ginit))
     +    task = 'CONVERGENCE'

c     Test for termination.

      if (task(1:4) .eq. 'WARN' .or. task(1:4) .eq. 'CONV') go to 10

c     A modified function is used to predict the step during the
c     first stage if a lower function value has been obtained but
c     the decrease is not sufficient.

      if (stage .eq. 1 .and. f .le. fx .and. f .gt. ftest) then

c        Define the modified function and derivative values.

         fm = f - stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest

c        Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,
     +               stmax)

c        Reset the function and derivative values for f.

         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      else

c       Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax)

      end if

c     Decide if a bisection step is needed.

      if (brackt) then
         if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty-stx)
         width1 = width
         width = abs(sty-stx)
      end if

c     Set the minimum and maximum steps allowed for stp.

      if (brackt) then
         stmin = min(stx,sty)
         stmax = max(stx,sty)
      else
         stmin = stp + xtrapl*(stp-stx)
         stmax = stp + xtrapu*(stp-stx)
      end if

c     Force the step to be within the bounds stpmax and stpmin.

      stp = max(stp,stpmin)
      stp = min(stp,stpmax)

c     If further progress is not possible, let stp be the best
c     point obtained during the search.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax) .or.
     +    (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx

c     Obtain another function and derivative.

      task = 'FG'

   10 continue

c     Save local variables.

      if (brackt) then
         isave(1) = 1
      else
         isave(1) = 0
      end if
      isave(2) = stage
      dsave(1) = ginit
      dsave(2) = gtest
      dsave(3) = gx
      dsave(4) = gy
      dsave(5) = finit
      dsave(6) = fx
      dsave(7) = fy
      dsave(8) = stx
      dsave(9) = sty
      dsave(10) = stmin
      dsave(11) = stmax
      dsave(12) = width
      dsave(13) = width1

      end

      subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,
     +                  stpmax)
      logical brackt
      double precision stx, fx, dx, sty, fy, dy, stp, fp, dp, stpmin,
     +                 stpmax
c     **********
c
c     Subroutine dcstep
c
c     This subroutine computes a safeguarded step for a search
c     procedure and updates an interval that contains a step that
c     satisfies a sufficient decrease and a curvature condition.
c
c     The parameter stx contains the step with the least function
c     value. If brackt is set to .true. then a minimizer has
c     been bracketed in an interval with endpoints stx and sty.
c     The parameter stp contains the current step.
c     The subroutine assumes that if brackt is set to .true. then
c
c           min(stx,sty) < stp < max(stx,sty),
c
c     and that the derivative at stx is negative in the direction
c     of the step.
c
c     The subroutine statement is
c
c       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
c                         stpmin,stpmax)
c
c     where
c
c       stx is a double precision variable.
c         On entry stx is the best step obtained so far and is an
c            endpoint of the interval that contains the minimizer.
c         On exit stx is the updated best step.
c
c       fx is a double precision variable.
c         On entry fx is the function at stx.
c         On exit fx is the function at stx.
c
c       dx is a double precision variable.
c         On entry dx is the derivative of the function at
c            stx. The derivative must be negative in the direction of
c            the step, that is, dx and stp - stx must have opposite
c            signs.
c         On exit dx is the derivative of the function at stx.
c
c       sty is a double precision variable.
c         On entry sty is the second endpoint of the interval that
c            contains the minimizer.
c         On exit sty is the updated endpoint of the interval that
c            contains the minimizer.
c
c       fy is a double precision variable.
c         On entry fy is the function at sty.
c         On exit fy is the function at sty.
c
c       dy is a double precision variable.
c         On entry dy is the derivative of the function at sty.
c         On exit dy is the derivative of the function at the exit sty.
c
c       stp is a double precision variable.
c         On entry stp is the current step. If brackt is set to .true.
c            then on input stp must be between stx and sty.
c         On exit stp is a new trial step.
c
c       fp is a double precision variable.
c         On entry fp is the function at stp
c         On exit fp is unchanged.
c
c       dp is a double precision variable.
c         On entry dp is the the derivative of the function at stp.
c         On exit dp is unchanged.
c
c       brackt is an logical variable.
c         On entry brackt specifies if a minimizer has been bracketed.
c            Initially brackt must be set to .false.
c         On exit brackt specifies if a minimizer has been bracketed.
c            When a minimizer is bracketed brackt is set to .true.
c
c       stpmin is a double precision variable.
c         On entry stpmin is a lower bound for the step.
c         On exit stpmin is unchanged.
c
c       stpmax is a double precision variable.
c         On entry stpmax is an upper bound for the step.
c         On exit stpmax is unchanged.
c
c     MINPACK-1 Project. June 1983
c     Argonne National Laboratory.
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision zero, p66, two, three
      parameter (zero=0.0d0,p66=0.66d0,two=2.0d0,three=3.0d0)

      double precision gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta

      sgnd = dp*(dx/abs(dx))

c     First case: A higher function value. The minimum is bracketed.
c     If the cubic step is closer to stx than the quadratic step, the
c     cubic step is taken, otherwise the average of the cubic and
c     quadratic steps is taken.

      if (fp .gt. fx) then
         theta = three*(fx-fp)/(stp-stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2-(dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma-dx) + theta
         q = ((gamma-dx)+gamma) + dp
         r = p/q
         stpc = stx + r*(stp-stx)
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/two)*(stp-stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
            stpf = stpc + (stpq-stpc)/two
         end if
         brackt = .true.

c     Second case: A lower function value and derivatives of opposite
c     sign. The minimum is bracketed. If the cubic step is farther from
c     stp than the secant step, the cubic step is taken, otherwise the
c     secant step is taken.

      else if (sgnd .lt. zero) then
         theta = three*(fx-fp)/(stp-stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2-(dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma-dp) + theta
         q = ((gamma-dp)+gamma) + dx
         r = p/q
         stpc = stp + r*(stx-stp)
         stpq = stp + (dp/(dp-dx))*(stx-stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
         end if
         brackt = .true.

c     Third case: A lower function value, derivatives of the same sign,
c     and the magnitude of the derivative decreases.

      else if (abs(dp) .lt. abs(dx)) then

c        The cubic step is computed only if the cubic tends to infinity
c        in the direction of the step or if the minimum of the cubic
c        is beyond stp. Otherwise the cubic step is defined to be the
c        secant step.

         theta = three*(fx-fp)/(stp-stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))

c        The case gamma = 0 only arises if the cubic does not tend
c        to infinity in the direction of the step.

         gamma = s*sqrt(max(zero,(theta/s)**2-(dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma-dp) + theta
         q = (gamma+(dx-dp)) + gamma
         r = p/q
         if (r .lt. zero .and. gamma .ne. zero) then
            stpc = stp + r*(stx-stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
         end if
         stpq = stp + (dp/(dp-dx))*(stx-stp)

         if (brackt) then

c           A minimizer has been bracketed. If the cubic step is
c           closer to stp than the secant step, the cubic step is
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .lt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            end if
            if (stp .gt. stx) then
               stpf = min(stp+p66*(sty-stp),stpf)
            else
               stpf = max(stp+p66*(sty-stp),stpf)
            end if
         else

c           A minimizer has not been bracketed. If the cubic step is
c           farther from stp than the secant step, the cubic step is
c           taken, otherwise the secant step is taken.

            if (abs(stpc-stp) .gt. abs(stpq-stp)) then
               stpf = stpc
            else
               stpf = stpq
            end if
            stpf = min(stpmax,stpf)
            stpf = max(stpmin,stpf)
         end if

c     Fourth case: A lower function value, derivatives of the same sign,
c     and the magnitude of the derivative does not decrease. If the
c     minimum is not bracketed, the step is either stpmin or stpmax,
c     otherwise the cubic step is taken.

      else
         if (brackt) then
            theta = three*(fp-fy)/(sty-stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2-(dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma-dp) + theta
            q = ((gamma-dp)+gamma) + dy
            r = p/q
            stpc = stp + r*(sty-stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
         end if
      end if

c     Update the interval which contains a minimizer.

      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. zero) then
            sty = stx
            fy = fx
            dy = dx
         end if
         stx = stp
         fx = fp
         dx = dp
      end if

c     Compute the new step.

      stp = stpf

      end


      subroutine dlmmv(n,m,s,y,rho,scale,mark,v,wa)
      integer n, m, mark
      double precision scale
      double precision s(n,m), y(n,m), rho(m), v(n), wa(m)
c     **********
c
c     This subroutine computes the matrix-vector product H*v
c     where H is the inverse BFGS approximation.
c
c     The matrix H depends on an initial matrix H0, m steps
c     s(1),...,s(m), and m gradient differences y(1),...,y(m).
c     These vectors are stored in the arrays s and y.
c     The most recent step and gradient difference are stored
c     in columns s(1,mark) and y(1,mark), respectively.
c     The initial matrix H0 is assumed to be scale*I.
c
c     The subroutine statement is
c
c       subroutine dlmmv(n,m,s,y,rho,scale,mark,v,wa)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       m is an integer variable.
c         On entry m specifies the number of steps and gradient
c            differences that are stored.
c         On exit m is unchanged.
c
c       s is a double precision array of dimension (n,m)
c         On entry s contains the m steps.
c         On exit s is unchanged.
c
c       y is a double precision array of dimension (n,m)
c         On entry y contains the m gradient differences.
c         On exit y is unchanged.
c
c       rho is a double precision array of dimension m
c         On entry rho contains the m innerproducts (s(i),y(i)).
c         On exit rho is unchanged.
c
c       scale is a double precision variable
c         On entry scale specifies the initial matrix H0 = scale*I.
c         On exit scale is unchanged.
c
c       mark is an integer variable.
c         On entry mark points to the current s(i) and y(i).
c         On exit mark is unchanged.
c
c       v is a double precision array of dimension n.
c         On entry v contains the vector v.
c         On exit v contains the matrix-vector product H*v.
c
c       wa is a double precision work array of dimension m.
c
c     Subprograms called
c
c       Level 1 BLAS ... daxpy, ddot, dscal
c
c     MINPACK-2 project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
c
c     **********
      integer i, k
      double precision beta

      double precision ddot
      external daxpy, ddot, dscal

      k = mark + 1
      do 10 i = 1, m
         k = k - 1
         if (k .eq. 0) k = m
         wa(k) = ddot(n,s(1,k),1,v,1)/rho(k)
         call daxpy(n,-wa(k),y(1,k),1,v,1)
   10 continue
      call dscal(n,scale,v,1)
      do 20 i = 1, m
         beta = wa(k) - ddot(n,y(1,k),1,v,1)/rho(k)
         call daxpy(n,beta,s(1,k),1,v,1)
         k = k + 1
         if (k .eq. m+1) k = 1
   20 continue

      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*), dy(*), da
      integer i, incx, incy, ix, iy, m, mp1, n
c
      if (n .le. 0) return
      if (da .eq. 0.0d0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         dy(iy) = dy(iy) + da*dx(ix)
         ix = ix + incx
         iy = iy + incy
   10 continue

      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 continue
      m = mod(n,4)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         dy(i) = dy(i) + da*dx(i)
   30 continue
      if (n .lt. 4) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 4
         dy(i) = dy(i) + da*dx(i)
         dy(i+1) = dy(i+1) + da*dx(i+1)
         dy(i+2) = dy(i+2) + da*dx(i+2)
         dy(i+3) = dy(i+3) + da*dx(i+3)
   50 continue

      return

      end

      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*), dy(*)
      integer i, incx, incy, ix, iy, m, mp1, n
c
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         dy(iy) = dx(ix)
         ix = ix + incx
         iy = iy + incy
   10 continue

      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 continue
      m = mod(n,7)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         dy(i) = dx(i)
   30 continue
      if (n .lt. 7) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 7
         dy(i) = dx(i)
         dy(i+1) = dx(i+1)
         dy(i+2) = dx(i+2)
         dy(i+3) = dx(i+3)
         dy(i+4) = dx(i+4)
         dy(i+5) = dx(i+5)
         dy(i+6) = dx(i+6)
   50 continue

      return

      end

      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*), dy(*), dtemp
      integer i, incx, incy, ix, iy, m, mp1, n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if (n .le. 0) return
      if (incx .eq. 1 .and. incy .eq. 1) go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1, n
         dtemp = dtemp + dx(ix)*dy(iy)
         ix = ix + incx
         iy = iy + incy
   10 continue
      ddot = dtemp

      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 continue
      m = mod(n,5)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if (n .lt. 5) go to 60
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 5
         dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     +           dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
   50 continue
   60 continue
      ddot = dtemp

      return

      end

      double precision function dnrm2 ( n, dx, incx)
      integer i, incx, ix, j, n, next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c     modified to correct failure to update ix, 1/25/92.
c     modified 3/93 to return if incx .le. 0.
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0 .and. incx.gt.0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      i = 1
      ix = 1
c                                                 begin main loop
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 continue
      ix = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j = ix,n
      if(dabs(dx(i)) .ge. hitest) go to 100
         sum = sum + dx(i)**2
         i = i + incx
   95 continue
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      ix = ix + 1
      i = i + incx
      if( ix .le. n ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end

      subroutine dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da, dx(1)
      integer i, incx, m, mp1, n, nincx
c
      if (n .le. 0 .or. incx .le. 0) return
      if (incx .eq. 1) go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1, nincx, incx
         dx(i) = da*dx(i)
   10 continue

      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 continue
      m = mod(n,5)
      if (m .eq. 0) go to 40
      do 30 i = 1, m
         dx(i) = da*dx(i)
   30 continue
      if (n .lt. 5) return
   40 continue
      mp1 = m + 1
      do 50 i = mp1, n, 5
         dx(i) = da*dx(i)
         dx(i+1) = da*dx(i+1)
         dx(i+2) = da*dx(i+2)
         dx(i+3) = da*dx(i+3)
         dx(i+4) = da*dx(i+4)
   50 continue

      return

      end
