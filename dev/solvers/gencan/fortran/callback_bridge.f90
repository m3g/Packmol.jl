!
! Part of the Packmol.jl GENCAN FFI wrapper (src/gencan/).
!
! Bridges GENCAN's evalal/evalnal (in gencan_core.f) to a single combined
! f+g callback supplied by Julia via a C function pointer. GENCAN's own
! active-set/line-search/CG internals are untouched; only the two small
! user-objective hooks are redirected here.
!
! Packmol.jl's fg! always computes f and g together in one pass, so both
! evalal (needs only f) and evalnal (needs only g) call the same callback
! and discard the unneeded output.
!
module callback_bridge_mod

   use iso_c_binding
   implicit none

   abstract interface
      function fg_callback_iface(n, x, f, g) result(flag) bind(c)
         import :: c_int, c_double
         integer(c_int), value :: n
         real(c_double), intent(in) :: x(n)
         real(c_double), intent(out) :: f
         real(c_double), intent(out) :: g(n)
         integer(c_int) :: flag
      end function fg_callback_iface
   end interface

   procedure(fg_callback_iface), pointer :: fg_callback_ptr => null()

   ! Module-level (not automatic/stack) scratch array for the gradient
   ! discarded by evalfg_bridge. At large n (e.g. ~900,000 for a
   ! 150,000-molecule system) an automatic array here would risk stack
   ! overflow since evalal is called repeatedly inside hot line-search
   ! loops; this is allocated once and only grown, never shrunk.
   double precision, allocatable :: g_scratch(:)

contains

   subroutine set_fg_callback(cptr) bind(c, name="set_fg_callback")
      type(c_funptr), value :: cptr
      call c_f_procpointer(cptr, fg_callback_ptr)
   end subroutine set_fg_callback

   ! Used by evalal: needs f, discards g into g_scratch.
   subroutine evalfg_bridge(n, x, f, flag)
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision, intent(out) :: f
      integer, intent(out) :: flag

      if (.not. allocated(g_scratch)) then
         allocate(g_scratch(n))
      else if (size(g_scratch) < n) then
         deallocate(g_scratch)
         allocate(g_scratch(n))
      end if

      flag = fg_callback_ptr(n, x, f, g_scratch(1:n))
   end subroutine evalfg_bridge

   ! Used by evalnal: needs g, discards f (a scalar, safe as a local var).
   subroutine evalgf_bridge(n, x, g, flag)
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      double precision, intent(out) :: g(n)
      integer, intent(out) :: flag
      double precision :: f_scratch

      flag = fg_callback_ptr(n, x, f_scratch, g)
   end subroutine evalgf_bridge

end module callback_bridge_mod

!
! Packmol's own gencan.f adds a call to an external `packmolprecision(n,x)`
! logical function (see gencan_core.f, "LM: Added packmolprecision function
! test") as an extra early-exit shortcut on top of GENCAN's normal stopping
! criteria. It has no equivalent split of Packmol's combined penalty into
! separate distance/constraint components in this wrapper, and it is purely
! an optional shortcut, not a correctness requirement — Packmol.jl's own
! outer loop (packmol_main.jl) already re-checks its own tol/constraint
! convergence between GENCAN calls. Stub it to never trigger; GENCAN falls
! through to its own standard stopping criteria immediately afterward.
!
function packmolprecision(n, x)
   implicit none
   integer :: n
   double precision :: x(n)
   logical :: packmolprecision
   packmolprecision = .false.
end function packmolprecision
