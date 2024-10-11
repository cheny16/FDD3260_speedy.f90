module physics
    use types, only: p
    use params

    implicit none

    private
    public initialize_physics, get_physical_tendencies

contains
    ! Initialize physical parametrization routines
    subroutine initialize_physics
        use physical_constants, only: grav, cp, p0, sigl, sigh, grdsig, grdscp, wvi
        use geometry, only: hsg, fsg, dhs

        integer :: k

        ! 1.2 Functions of sigma and latitude
        sigh(0) = hsg(1)

        do k = 1, kx
            sigl(k) = log(fsg(k))
            sigh(k) = hsg(k+1)
            grdsig(k) = grav/(dhs(k)*p0)
            grdscp(k) = grdsig(k)/cp
        end do

        ! Weights for vertical interpolation at half-levels(1,kx) and surface
        ! Note that for phys.par. half-lev(k) is between full-lev k and k+1
        ! Fhalf(k) = Ffull(k)+WVI(K,2)*(Ffull(k+1)-Ffull(k))
        ! Fsurf = Ffull(kx)+WVI(kx,2)*(Ffull(kx)-Ffull(kx-1))
        do k = 1, kx-1
            wvi(k,1) = 1./(sigl(k+1)-sigl(k))
            wvi(k,2) = (log(sigh(k))-sigl(k))*wvi(k,1)
        end do

        wvi(kx,1) = 0.
        wvi(kx,2) = (log(0.99)-sigl(kx))*wvi(kx-1,1)
    end

    !> Compute physical parametrization tendencies for u, v, t, q and add them
    !  to the dynamical grid-point tendencies
    subroutine get_physical_tendencies(vor, div, t, q, phi, psl, utend, vtend, ttend, qtend)
        use auxiliaries, only: precnv, precls, cbmf, tsr, ssrd, ssr, slrd, slr, olr, slru, ustr, &
            & vstr, shf, evap, hfluxn
        use physical_constants, only: sigh, grdsig, grdscp, cp
        use geometry, only: fsg
        use boundaries, only: phis0
        use land_model, only: fmask_l
        use sea_model, only: sst_am, ssti_om, sea_coupling_flag
        use sppt, only: mu, gen_sppt
        use convection, only: get_convection_tendencies
        use large_scale_condensation, only: get_large_scale_condensation_tendencies
        use shortwave_radiation, only: get_shortwave_rad_fluxes, clouds, compute_shortwave
        use longwave_radiation, only: &
                get_downward_longwave_rad_fluxes, get_upward_longwave_rad_fluxes
        use surface_fluxes, only: get_surface_fluxes
        use vertical_diffusion, only: get_vertical_diffusion_tend
        use humidity, only: spec_hum_to_rel_hum
        use spectral, only: spec_to_grid, uvspec

        complex(p), intent(in) :: vor(mx,nx,kx) !! Vorticity
        complex(p), intent(in) :: div(mx,nx,kx) !! Divergence
        complex(p), intent(in) :: t(mx,nx,kx)   !! Temperature
        complex(p), intent(in) :: q(mx,nx,kx)   !! Specific Humidity
        complex(p), intent(in) :: phi(mx,nx,kx) !! Geopotential
        complex(p), intent(in) :: psl(mx,nx)    !! ln(Surface pressure)

        real(p), intent(inout) :: utend(ix,il,kx) !! Zonal velocity tendency
        real(p), intent(inout) :: vtend(ix,il,kx) !! Meridional velocity tendency
        real(p), intent(inout) :: ttend(ix,il,kx) !! Temperature tendency
        real(p), intent(inout) :: qtend(ix,il,kx) !! Specific humidity tendency

        complex(p), dimension(mx,nx) :: ucos, vcos
        real(p), dimension(ix,il) :: pslg, rps, gse
        real(p), dimension(ix,il,kx) :: ug, vg, tg, qg, phig, utend_dyn, vtend_dyn, ttend_dyn, qtend_dyn
        real(p), dimension(ix,il,kx) :: se, rh, qsat
        real(p), dimension(ix,il) :: psg, ts, tskin, u0, v0, t0, cloudc, clstr, cltop, prtop
        real(p), dimension(ix,il,kx) :: tt_cnv, qt_cnv, tt_lsc, qt_lsc, tt_rsw, tt_rlw, ut_pbl, vt_pbl,&
            & tt_pbl, qt_pbl
        integer :: iptop(ix,il), icltop(ix,il,2), icnv(ix,il), i, j, k
        real(p) :: sppt_pattern(ix,il,kx)

        ! Keep a copy of the original (dynamics only) tendencies
        utend_dyn = utend
        vtend_dyn = vtend
        ttend_dyn = ttend
        qtend_dyn = qtend

        ! =========================================================================
        ! Compute grid-point fields
        ! =========================================================================
        
        ! Convert model spectral variables to grid-point variables
        
        !$omp parallel default(shared) private(k, i, j, ucos, vcos)

        !$omp do schedule(guided)
        do k = 1, kx
            call uvspec(vor(:,:,k), div(:,:,k), ucos, vcos)
            ug(:,:,k)   = spec_to_grid(ucos, 2)
            vg(:,:,k)   = spec_to_grid(vcos, 2)
            tg(:,:,k)   = spec_to_grid(t(:,:,k), 1)
            qg(:,:,k)   = spec_to_grid(q(:,:,k), 1)
            phig(:,:,k) = spec_to_grid(phi(:,:,k), 1)
        end do
        !$omp end do
        
        !$omp single
        pslg = spec_to_grid(psl, 1)
        !$omp end single

        ! =========================================================================
        ! Compute thermodynamic variables
        ! =========================================================================

        !$omp single
        psg = exp(pslg)
        rps = 1.0/psg

        qg = max(qg, 0.0)
        se = cp*tg + phig
        !$omp end single

        !$omp do schedule(static)
        do k = 1, kx
            call spec_hum_to_rel_hum(tg(:,:,k), psg, fsg(k), qg(:,:,k), rh(:,:,k), qsat(:,:,k))
        end do
        !$omp end do

        ! =========================================================================
        ! Precipitation
        ! =========================================================================

        ! Deep convection

        !$omp single
        call get_convection_tendencies(psg, se, qg, qsat, iptop, cbmf, precnv, tt_cnv, qt_cnv)
        !$omp end single

        !$omp do schedule(static)
        do k = 2, kx
            tt_cnv(:,:,k) = tt_cnv(:,:,k)*rps*grdscp(k)
            qt_cnv(:,:,k) = qt_cnv(:,:,k)*rps*grdsig(k)
        end do
        !$omp end do

        icnv = kx - iptop

        ! Large-scale condensation

        !$omp single
        call get_large_scale_condensation_tendencies(psg, qg, qsat, iptop, precls, tt_lsc, qt_lsc)
        !$omp end single

        !$omp single
            !$omp task
            ttend = ttend + tt_cnv + tt_lsc
            !$omp end task

            !$omp task
            qtend = qtend + qt_cnv + qt_lsc
            !$omp end task
        !$omp end single nowait

        ! =========================================================================
        ! Radiation (shortwave and longwave) and surface fluxes
        ! =========================================================================

        ! Compute shortwave tendencies and initialize lw transmissivity
        ! The shortwave radiation may be called at selected time steps
        if (compute_shortwave) then
            !$omp single
            gse = (se(:,:,kx-1) - se(:,:,kx))/(phig(:,:,kx-1) - phig(:,:,kx))

            call clouds(qg, rh, precnv, precls, iptop, gse, fmask_l, icltop, cloudc, clstr)

            do i = 1, ix
                do j = 1, il
                    cltop(i,j) = sigh(icltop(i,j,1) - 1)*psg(i,j)
                    prtop(i,j) = float(iptop(i,j))
                end do
            end do

            call get_shortwave_rad_fluxes(psg, qg, icltop, cloudc, clstr, ssrd, ssr, tsr, tt_rsw)

            do k = 1, kx
                tt_rsw(:,:,k) = tt_rsw(:,:,k)*rps*grdscp(k)
            end do
            !$omp end single
        end if

        ! Compute downward longwave fluxes
        
        !$omp single
        call get_downward_longwave_rad_fluxes(tg, slrd, tt_rlw)
        !$omp end single

        ! Compute surface fluxes and land skin temperature
        
        !$omp single
        call get_surface_fluxes(psg, ug, vg, tg, qg, rh, phig, phis0, fmask_l, sst_am, &
                & ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, .true.)
        !$omp end single
        
        ! Recompute sea fluxes in case of anomaly coupling
        if (sea_coupling_flag > 0) then
            !$omp single
            call get_surface_fluxes(psg, ug, vg, tg, qg, rh, phig, phis0, fmask_l, ssti_om, &
                   & ssrd, slrd, ustr, vstr, shf, evap, slru, hfluxn, ts, tskin, u0, v0, t0, .false.)
            !$omp end single
        end if

        ! Compute upward longwave fluxes, convert them to tendencies and add
        ! shortwave tendencies

        !$omp single
        call get_upward_longwave_rad_fluxes(tg, ts, slrd, slru(:,:,3), slr, olr, tt_rlw)

        do k = 1, kx
            tt_rlw(:,:,k) = tt_rlw(:,:,k)*rps*grdscp(k)
        end do

        ttend = ttend + tt_rsw + tt_rlw
        !$omp end single

        ! =========================================================================
        ! Planetary boundary later interactions with lower troposphere
        ! =========================================================================

        ! Vertical diffusion and shallow convection

        !$omp single
        call get_vertical_diffusion_tend(se, rh, qg, qsat, phig, icnv, ut_pbl, vt_pbl, &
            & tt_pbl, qt_pbl)
        !$omp end single

        ! Add tendencies due to surface fluxes

#ifdef _OPENMP
        !$omp do schedule(static)
        do k = 1, kx
            if (k == kx) then
                ut_pbl(:,:,k) = ut_pbl(:,:,k) + ustr(:,:,3) * rps * grdsig(k)
                vt_pbl(:,:,k) = vt_pbl(:,:,k) + vstr(:,:,3) * rps * grdsig(k)
                tt_pbl(:,:,k) = tt_pbl(:,:,k) + shf(:,:,3)  * rps * grdscp(k)
                qt_pbl(:,:,k) = qt_pbl(:,:,k) + evap(:,:,3) * rps * grdsig(k)
            end if
        end do
        !$omp end do

        !$omp do schedule(static)
        do k = 1, kx
            utend(:,:,k) = utend(:,:,k) + ut_pbl(:,:,k)
            vtend(:,:,k) = vtend(:,:,k) + vt_pbl(:,:,k)
            ttend(:,:,k) = ttend(:,:,k) + tt_pbl(:,:,k)
            qtend(:,:,k) = qtend(:,:,k) + qt_pbl(:,:,k)
        end do
        !$omp end do
#else
        ut_pbl(:,:,kx) = ut_pbl(:,:,kx) + ustr(:,:,3)*rps*grdsig(kx)
        vt_pbl(:,:,kx) = vt_pbl(:,:,kx) + vstr(:,:,3)*rps*grdsig(kx)
        tt_pbl(:,:,kx) = tt_pbl(:,:,kx)  + shf(:,:,3)*rps*grdscp(kx)
        qt_pbl(:,:,kx) = qt_pbl(:,:,kx) + evap(:,:,3)*rps*grdsig(kx)

        utend = utend + ut_pbl
        vtend = vtend + vt_pbl
        ttend = ttend + tt_pbl
        qtend = qtend + qt_pbl
#endif     
           
        ! !$omp single
        !     !$omp task depend(out: ut_pbl)
        !     ut_pbl(:,:,kx) = ut_pbl(:,:,kx) + ustr(:,:,3)*rps*grdsig(kx)
        !     !$omp end task

        !     !$omp task depend(out: vt_pbl)
        !     vt_pbl(:,:,kx) = vt_pbl(:,:,kx) + vstr(:,:,3)*rps*grdsig(kx)
        !     !$omp end task
            
        !     !$omp task depend(out: tt_pbl)
        !     tt_pbl(:,:,kx) = tt_pbl(:,:,kx)  + shf(:,:,3)*rps*grdscp(kx)
        !     !$omp end task

        !     !$omp task depend(out: qt_pbl)
        !     qt_pbl(:,:,kx) = qt_pbl(:,:,kx) + evap(:,:,3)*rps*grdsig(kx)
        !     !$omp end task

        !     !$omp task depend(in: ut_pbl)
        !     utend = utend + ut_pbl
        !     !$omp end task

        !     !$omp task depend(in: vt_pbl)
        !     vtend = vtend + vt_pbl
        !     !$omp end task

        !     !$omp task depend(in: tt_pbl)
        !     ttend = ttend + tt_pbl
        !     !$omp end task

        !     !$omp task depend(in: qt_pbl)
        !     qtend = qtend + qt_pbl
        !     !$omp end task
        ! !$omp end single

        ! Add SPPT noise
        if (sppt_on) then
            !$omp single
            sppt_pattern = gen_sppt()
            !$omp end single

            ! The physical contribution to the tendency is *tend - *tend_dyn, where * is u, v, t, q

            !$omp do schedule(static) private(k)
            do k = 1,kx
                utend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(utend(:,:,k) - utend_dyn(:,:,k)) &
                        & + utend_dyn(:,:,k)
                vtend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(vtend(:,:,k) - vtend_dyn(:,:,k)) &
                        & + vtend_dyn(:,:,k)
                ttend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(ttend(:,:,k) - ttend_dyn(:,:,k)) &
                        & + ttend_dyn(:,:,k)
                qtend(:,:,k) = (1 + sppt_pattern(:,:,k)*mu(k))*(qtend(:,:,k) - qtend_dyn(:,:,k)) &
                        & + qtend_dyn(:,:,k)
            end do
            !$omp end do
        end if
        !$omp end parallel
    end
end module