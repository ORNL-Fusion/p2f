!> Namelist defaults.
!! \author DLG
!<

module read_namelist

implicit none
integer :: MaxSteps = 10000
integer, parameter :: LONG = selected_int_kind(9)
integer(kind=LONG) :: nP = 1000 !< Number of particles to use.
real :: particleSize = 10.0e-04 !< Gaussian particle size in \%c. 
logical :: gParticle = .false. !< Use gaussian particle shape or not.
character(len=100) :: pl_fileName = 'data/fdis_25keV_D_6_flat_D3D.dav' !<Input particle list filename.
character(len=100) :: eqdsk_fileName = 'data/g122080.03100' !<Input gEqDsk filename.
real :: vPerInRange = 2.0 !<vPer grid range in \%c
real :: vParInRange = 2.0 !<vPar grid range in \%c
integer :: nGridPtsGaussianVSpace = 6 !<Number of grid points to span half gaussian in velocity. Note: choosing a value to small may give discontiuities for large particleSize values.
integer :: nGridPtsGaussianXSpace = 4 !<Number of grid points to span half gaussian in real space. Note: choosing a value to small may give discontiuities for large particleSize values.
integer :: R_nBins = 10 !<Number of bins in R (Major radius).
integer :: z_nBins = 20 !<Number of bins in z (Vertical).
integer :: vPer_nBins = 64 !<Number of vPer bins. Note: set this value low (~16) when gParticle = .false.
integer :: vPar_nBins = 128 !<Number of vPar bins. Note: set this value low (~32) when gParticle = .false.
real :: amu = 1.0 !< 1.0 for H, 2.0 for D, etc ...
real :: AtomicZ = 1.0 !< 1.0 for H, 2.0 for D, etc ...
logical :: plotOrbit = .false. !<Use DISLIN to plot the orbits when DistributeAlongOrbit = .false. but this only works on dlg-hp now and is primarily for debugging. So leave it off.
logical :: NFO = .false. !<Compute orbits without finite orbits.
logical :: DistributeAlongOrbit = .false. !<No orbit averaging, just particles where they are.
logical :: pList_is_nCDF = .true. !<If using a DLG netCDF input particle list set this. This is used in the sMC-rf code or when re-mapping the ORBIT-rf output on a high resolution sMap grid.
logical :: pList_is_nubeam1 = .false. !<Use a nubeam file
real :: weightLimit = 1.1e18 !<If the particle weight is above this number the particle is ignored.
real :: traceStepLength = 0.002 !< Approximate step length [m] of integrator in m.
integer :: usePtsPerTrack = 50 !< Approximate number of points to use from each track in the orbit averaging.
real :: energyThreshold = 0.0 !< [keV] lower energy limit of particles not to include
 
namelist / P2FIN / nP,&
    gParticle, &
    particleSize, &
    pl_fileName, &
    eqdsk_fileName, &
    vPerInRange, &
    vParInRange, &
    nGridPtsGaussianVSpace, &
    nGridPtsGaussianXSpace, &
    R_nBins, &
    z_nBins, &
    vPer_nBins, &
    vPar_nBins, &
    amu, &
    AtomicZ, &
    plotOrbit, &
    NFO, &
    DistributeAlongOrbit, &
    pList_is_nCDF, & 
    weightLimit, &
    traceStepLength, &
    usePtsPerTrack, &
    pList_is_nubeam1, &
    energyThreshold, &
    MaxSteps

contains
subroutine init_namelist

    implicit none
    character(len=100) :: nml_fileName

    nml_fileName    = 'p2f.nml'
    
    open ( 7, file = nml_fileName, delim = 'APOSTROPHE' )
    read ( unit = 7, nml = P2FIN )
    close ( 7 )

end subroutine init_namelist

end module read_namelist
