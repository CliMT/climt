! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! @brief Calculate ancillary horizon angles on LFRic grid

program calc_horizon_angles

  use netcdf
  use realtype_rd, only: RealK
  use rad_ccf, only: pi
  use socrates_horizon_2d, only: horizon_2d
  implicit none

  character(*), parameter :: orog_name='unfiltered_surface_altitude'
  character(*), parameter :: alt_orog_name='m01s00i007'
  character(*), parameter :: filt_orog_name='surface_altitude'
  character(*), parameter :: lon_name_ext='_face_x'
  character(*), parameter :: lat_name_ext='_face_y'
  character(*), parameter :: node_lon_name_ext='_node_x'
  character(*), parameter :: node_lat_name_ext='_node_y'
  character(*), parameter :: node_map_name_ext='_face_nodes'
  character(*), parameter :: hor_ang_name='horizon_angle'
  character(*), parameter :: hor_asp_name='horizon_aspect'
  character(*), parameter :: panel_num(6)=['1','2','3','4','5','6']
  character(*), parameter :: dir_num(16)=['01','02','03','04','05','06', &
                      '07','08','09','10','11','12','13','14','15','16']

  character(len=1024) :: file_orog, mesh_name
  character(len=18) :: horizon_angle_str
  character(len=19) :: horizon_aspect_str

  integer :: n_horiz_ang = 16
  integer :: n_horiz_layer = 1
  integer :: horiz_limit
  integer :: num_h_ang, num_h_asp

  integer :: i, j, k, l, n, x, y
  integer :: node, face, face_next, quadrant, corner_node, panel
  integer :: node_x1, node_x2, node_y1, node_y2, first_node_y2
  integer :: next_x1, next_x2, next_y1, next_y2
  integer :: x_inc, y_inc, x_offset, y_offset
  integer :: n_unmapped_faces
  integer :: x_map_missing

  integer :: status, ncid, varid
  integer :: dimid0, dimid1, dimid2, dimid3, dimid4, dimid5, dimid6, dimid7
  integer :: dimids(7)
  integer :: n_profile, n_node, n_edge, n_x, n_y, n_panel
  integer :: corner_quadrant, node_quadrant(4), lower
  integer, allocatable :: x_map(:), y_map(:), x_panel_map(:), y_panel_map(:)
  integer, allocatable :: node_map(:, :), face_map(:, :)
  integer, allocatable :: n_connected_faces(:)
  real(RealK), allocatable :: surface_altitude(:), latitude(:), longitude(:)
  real(RealK), allocatable :: node_latitude(:), node_longitude(:)
  real(RealK), allocatable :: surface_altitude_2d(:, :)
  real(RealK), allocatable :: latitude_2d(:, :), longitude_2d(:, :)
  real(RealK), allocatable :: surface_altitude_work(:, :)
  real(RealK), allocatable :: latitude_work(:, :), longitude_work(:, :)
  real(RealK), allocatable :: horizon_angle(:, :)
  real(RealK), allocatable :: horizon_aspect(:, :)
  real(RealK), allocatable :: horizon_angle_2d(:, :, :, :)
  real(RealK), allocatable :: horizon_aspect_2d(:, :, :)
  real(RealK) :: node_aspect(4), node_lat, node_lon, face_lat, face_lon
  real(RealK) :: eqt_res, tan_i, tan_j, halo_i, weight_lower, weight_upper

  real(RealK), parameter :: pi_over_2 = pi/2.0_RealK
  real(RealK), parameter :: pi_over_4 = pi/4.0_RealK
  real(RealK), parameter :: pi_over_180 = pi/180.0_RealK
  real(RealK), parameter :: twopi = pi*2.0_RealK

  logical :: node_mask(4)
  logical :: l_local, l_final_panel, l_output_panel_fields, l_open

  ! Open orography ancillary file read/write
  l_open = .false.
  l_output_panel_fields = .false.
  do i = 1, command_argument_count()
    call get_command_argument(i, file_orog)
    select case (file_orog)
    case ('-x')
      l_output_panel_fields = .true.
    case ('-h', '--help')
      write(*, *) 'Usage: calc_horizon_angles [-x] <orog netcdf ancillary>'
      write(*, *) '  -x outputs diagnostic 2D fields to netcdf file'
      stop
    case default
      status = nf90_open(trim(file_orog), NF90_WRITE, ncid)
      if (status == NF90_NOERR) then
        l_open = .true.
      else
        write(*,*) 'ERROR reading netCDF file: ',nf90_strerror(status)
        stop
      end if
    end select
  end do
  if (.not.l_open) then
    write(*, *) 'Usage: calc_horizon_angles <orog netcdf ancillary>'
    stop
  end if

  ! See if horizon angles are already present
  status = nf90_inq_varid(ncid, hor_ang_name, varid)
  if (status == NF90_NOERR) then
    call nf(nf90_close(ncid))
    stop 'Horizon angles already present: will not be overwritten'
  end if

  ! Get unfiltered surface altitude
  status = nf90_inq_varid(ncid, orog_name, varid)
  if (status /= NF90_NOERR) then
    status = nf90_inq_varid(ncid, alt_orog_name, varid)
    if (status /= NF90_NOERR) then
      write(*, *) "Can't find unfiltered surface altitude field"
      status = nf90_inq_varid(ncid, filt_orog_name, varid)
      if (status /= NF90_NOERR) then
        stop "Can't find filtered surface altitude field"
      else
        write(*, *) "Using filtered surface altitude field"
      end if
    end if
  end if
  ! Get number of points
  call nf(nf90_inquire_variable(ncid, varid, dimids=dimids))
  dimid1 = dimids(1)
  call nf(nf90_inquire_dimension(ncid, dimid1, len=n_profile))
  allocate(surface_altitude(n_profile))
  status = nf90_get_var(ncid, varid, surface_altitude)

  ! Get name of mesh
  status = nf90_get_att(ncid, varid, 'mesh', mesh_name)
  if (status /= NF90_NOERR) then
    stop "Can't find mesh for surface altitude"
  end if

  ! Get longitude of faces
  status = nf90_inq_varid(ncid, trim(mesh_name)//lon_name_ext, varid)
  if (status /= NF90_NOERR) then
    stop "Can't find longitude field"
  end if
  allocate(longitude(n_profile))
  call nf(nf90_get_var(ncid, varid, longitude))

  ! Get latitude of faces
  status = nf90_inq_varid(ncid, trim(mesh_name)//lat_name_ext, varid)
  if (status /= NF90_NOERR) then
    stop "Can't find latitude field"
  end if
  allocate(latitude(n_profile))
  call nf(nf90_get_var(ncid, varid, latitude))

  ! Get longitude of nodes
  status = nf90_inq_varid(ncid, trim(mesh_name)//node_lon_name_ext, varid)
  if (status /= NF90_NOERR) then
    stop "Can't find node longitude field"
  end if
  ! Get number of nodes
  call nf(nf90_inquire_variable(ncid, varid, dimids=dimids))
  dimid0 = dimids(1)
  call nf(nf90_inquire_dimension(ncid, dimid0, len=n_node))
  allocate(node_longitude(n_node))
  call nf(nf90_get_var(ncid, varid, node_longitude))

  ! Get latitude of nodes
  status = nf90_inq_varid(ncid, trim(mesh_name)//node_lat_name_ext, varid)
  if (status /= NF90_NOERR) then
    stop "Can't find node latitude field"
  end if
  allocate(node_latitude(n_node))
  call nf(nf90_get_var(ncid, varid, node_latitude))

  ! Get node map (map of nodes connected to each face)
  status = nf90_inq_varid(ncid, trim(mesh_name)//node_map_name_ext, varid)
  if (status /= NF90_NOERR) then
    stop "Can't find node map (face_node_connectivity)"
  end if
  allocate(node_map(4, n_profile))
  status = nf90_get_var(ncid, varid, node_map)

  ! Calculate face map (map of faces connected to each node)
  allocate(face_map(4, n_node))
  allocate(n_connected_faces(n_node))
  face_map(1:4, 1:n_node) = 0
  n_connected_faces(1:n_node) = 0
  do l=1, n_profile
    do n=1, 4
      node = node_map(n, l)
      n_connected_faces(node) = n_connected_faces(node) + 1
      face_map(n_connected_faces(node), node) = l
    end do
  end do

  ! Construct 2D map of faces
  allocate(x_map(n_profile))
  allocate(y_map(n_profile))
  x_map_missing = n_profile + 1
  x_map(1:n_profile) = n_profile + 1
  y_map(1:n_profile) = n_profile + 1

  ! Check if this is a local area grid or cubed sphere
  if (n_node - n_profile > 2) then
    ! Assume 2D local area grid
    l_local = .true.
    n_edge = n_node - n_profile - 1
    n_panel = 1
    ! Initialise n_x, n_y to max values to be reset later
    n_x = n_edge
    n_y = n_edge
  else
    ! Assume cubed sphere grid
    l_local = .false.
    l_final_panel = .false.
    n_edge = nint(sqrt(real(n_profile/6, RealK)))
    eqt_res = pi / real(2*n_edge, RealK)
    n_x = n_edge*4
    n_y = n_edge*3
    n_panel = 6
  end if

  ! Find first corner node and face
  do node=1, n_node
    if (n_connected_faces(node) == 1 .or. &
        n_connected_faces(node) == 3) then
      face=face_map(1, node)
      x_map(face) = 0
      y_map(face) = 0
      corner_node = node
      exit
    end if
    if (node == n_node) stop "Can't find corner node."
  end do

  ! Orient nodes of first face
  face_lat = latitude(face) * pi_over_180
  face_lon = longitude(face) * pi_over_180
  do n=1, 4
    node = node_map(n, face)
    node_lat = node_latitude(node) * pi_over_180
    node_lon = node_longitude(node) * pi_over_180
    ! Find compass bearing of nodes from face centre
    node_aspect(n) = modulo( atan2( &
        cos(node_lat) * sin(node_lon-face_lon), &
        sin(node_lat) * cos(face_lat) &
      - cos(node_lat) * sin(face_lat) * cos(node_lon-face_lon) ), twopi )
  end do
  deallocate(node_latitude)
  deallocate(node_longitude)

  ! Map nodes of first face onto four quadrants
  !     N
  !   4 | 1
  ! W --|-- E
  !   3 | 2
  !     S
  node_mask(1:4) = .true.
  do quadrant=1, 4
    n = minloc(node_aspect, 1, node_mask)
    node = node_map(n, face)
    node_quadrant(quadrant) = node
    if (node == corner_node) corner_quadrant = quadrant
    node_mask(n) = .false.
  end do

  ! Set increments and connecting nodes to loop away from corner
  node_y1 = corner_node
  select case (corner_quadrant)
  case(1)
    node_y2 = node_quadrant(4)
    next_y1 = node_quadrant(2)
    next_y2 = node_quadrant(3)
    x_inc = -1
    y_inc = -1
  case(2)
    node_y2 = node_quadrant(3)
    next_y1 = node_quadrant(1)
    next_y2 = node_quadrant(4)
    x_inc = -1
    y_inc =  1
  case(3)
    node_y2 = node_quadrant(2)
    next_y1 = node_quadrant(4)
    next_y2 = node_quadrant(1)
    x_inc =  1
    y_inc =  1
  case(4)
    node_y2 = node_quadrant(1)
    next_y1 = node_quadrant(3)
    next_y2 = node_quadrant(2)
    x_inc =  1
    y_inc = -1
  end select
  first_node_y2 = node_y2

  ! Loop in the N-S direction
  y = 0
  y_loop: do
    if (y /= 0) then

      ! Find next face from connecting nodes
      if (node_y1 == 0 .or. node_y2 == 0) stop "Can't find connecting y-nodes"
      y_faces: do l=1, n_connected_faces(node_y1)
        face_next = face_map(l, node_y1)
        if (any(face_next &
             == face_map(1:n_connected_faces(node_y2), node_y2))) then
          ! See if common face of these nodes has already been mapped
          if (x_map(face_next) == n_profile + 1) then
            ! If not then this is the next face along
            face = face_next
            x_map(face) = 0
            y_map(face) = y
          end if
        else if (x_map(face_next) == n_profile + 1) then
          ! This must be the last face of the row for a cubed sphere
          ! Set the y_map for this face so we can test on mapped faces later
          y_map(face_next) = y
        end if
      end do y_faces

      ! Map next connecting nodes
      next_y1 = 0
      next_y2 = 0
      do n=1, 4
        node = node_map(n, face)
        if (node /= node_y1 .and. node /= node_y2) then
          if (l_local) then
            if (n_connected_faces(node) == 1) then
              ! Corner node
              next_y1 = node
            else if (abs(y+y_inc) == n_y) then
              ! End of column but not corner
              next_y2 = node
            else if (n_connected_faces(node) == 2) then
              ! Edge node
              next_y1 = node
            else
              ! Internal node
              next_y2 = node
            end if
          else
            ! Cubed sphere
            if (n_connected_faces(node) == 3) then
              ! Corner node
              next_y1 = node
            else if (abs(y+y_inc) == n_edge*2 .or. &
                    (abs(y+y_inc) == n_edge+1 .and. l_final_panel)) then
              ! End of column but not corner
              next_y2 = node
            else
              n_unmapped_faces=0
              do l=1, n_connected_faces(node)
                if (y_map(face_map(l, node)) == n_profile + 1) then
                  n_unmapped_faces=n_unmapped_faces+1
                end if
              end do
              if (n_unmapped_faces <= 2) then
                ! Edge node
                next_y1 = node
              else
                ! Internal node
                next_y2 = node
              end if
            end if
          end if
        end if
      end do

    end if

    ! Loop in the E-W direction
    node_x1 = node_y2
    node_x2 = next_y2
    x = 0
    x_loop: do

      ! Find next face from connecting nodes
      if (node_x1 == 0 .or. node_x2 == 0) stop "Can't find connecting x-nodes"
      x_faces: do l=1, n_connected_faces(node_x1)
        face_next = face_map(l, node_x1)
        if (any(face_next &
             == face_map(1:n_connected_faces(node_x2), node_x2))) then
          ! See if common face of these nodes has already been mapped
          if (x_map(face_next) == n_profile + 1) then
            ! If not then this is the next face along
            face = face_next
            x_map(face) = x+x_inc
            y_map(face) = y
            exit x_faces
          end if
        end if
      end do x_faces

      ! Map next connecting nodes
      next_x1 = 0
      next_x2 = 0
      do n=1, 4
        node = node_map(n, face)
        if (node /= node_x1 .and. node /= node_x2) then
          ! It doesn't matter which way round the next nodes are
          if (next_x1 == 0) then
            next_x1 = node
          else
            next_x2 = node
          end if
          if (n_connected_faces(node) == 1 .and. y == 0) then
            ! Reached the end of the first row for the local area grid
            ! We can now define the x, y dimensions
            n_x = abs(x+x_inc*2)
            n_y = n_edge - n_x
          end if
        end if
      end do

      x = x + x_inc
      if (abs(x+x_inc) == n_x) exit x_loop
      node_x1 = next_x1
      node_x2 = next_x2
    end do x_loop

    y = y + y_inc
    if (l_local) then
      if (abs(y) == n_y) exit y_loop
    else
      ! For the cubed sphere we first loop across four panels,
      ! reduce the x-dimension to loop across the fifth,
      ! then adjust the starting point and increment for the final panel.
      if (abs(y) == n_edge) n_x = n_edge
      if (abs(y) == n_edge*2) then
        l_final_panel = .true.
        y = -y_inc
        y_inc = -y_inc
        next_y1 = corner_node
        next_y2 = first_node_y2
      end if
      if (l_final_panel .and. abs(y) == n_edge+1) then
        ! Finished cubed-sphere. Reset grid size to single panel and exit loop.
        n_x = n_edge
        n_y = n_edge
        exit y_loop
      end if
    end if
    node_y1 = next_y1
    node_y2 = next_y2
  end do y_loop
  deallocate(face_map)
  deallocate(node_map)

  ! Offset maps so that all values are positive
  x_offset = 1 - minval(x_map)
  y_offset = 1 - minval(y_map)
  x_map(1:n_profile) = x_map(1:n_profile) + x_offset
  y_map(1:n_profile) = y_map(1:n_profile) + y_offset
  x_map_missing = x_map_missing + x_offset

  ! Set horiz_limit (should ideally be done using max orography and distance)
  if (l_local) then
    horiz_limit=min(n_x, n_y, 100)
  else
    horiz_limit=n_edge/10
  end if

  ! Create fields for horizon angles and aspects
  num_h_ang = n_horiz_ang*n_horiz_layer
  num_h_asp = n_horiz_ang
  allocate(horizon_angle(n_profile, num_h_ang))
  allocate(horizon_aspect(n_profile, num_h_asp))
  horizon_angle = 0.0_RealK
  horizon_aspect = -1.0_RealK

  ! Add dimensions to the orog ancillary file
  call nf(nf90_def_dim(ncid, 'num_h_ang', num_h_ang, dimid2))
  call nf(nf90_def_dim(ncid, 'num_h_asp', num_h_asp, dimid3))
  if (l_output_panel_fields) then
    call nf(nf90_def_dim(ncid, 'num_x', n_x, dimid4))
    call nf(nf90_def_dim(ncid, 'num_y', n_y, dimid5))
    call nf(nf90_def_dim(ncid, 'num_x_halo', n_x+2*horiz_limit, dimid6))
    call nf(nf90_def_dim(ncid, 'num_y_halo', n_y+2*horiz_limit, dimid7))
  end if

  ! Allocate 2D fields
  allocate(surface_altitude_2d(1-horiz_limit:n_x+horiz_limit, &
                               1-horiz_limit:n_y+horiz_limit))
  allocate(longitude_2d(1-horiz_limit:n_x+horiz_limit, &
                        1-horiz_limit:n_y+horiz_limit))
  allocate(latitude_2d(1-horiz_limit:n_x+horiz_limit, &
                       1-horiz_limit:n_y+horiz_limit))
  allocate(horizon_angle_2d(n_x, n_y, n_horiz_layer, n_horiz_ang))
  allocate(horizon_aspect_2d(n_x, n_y, n_horiz_ang))

  ! Allocate maps for the separate panels
  allocate(x_panel_map(n_profile))
  allocate(y_panel_map(n_profile))

  ! Calculate horizon angles for each panel separately
  panel_loop: do panel = 1, n_panel

    ! Initialise 2D fields and halos to zero for each panel
    surface_altitude_2d=0.0_RealK
    longitude_2d=0.0_RealK
    latitude_2d=0.0_RealK

    ! Set mapping for each panel
    if (l_local) then
      x_panel_map = x_map
      y_panel_map = y_map
    else
      x_panel_map(1:n_profile) = x_map_missing
      select case (panel)
      ! Equatorial panels 1-4
      case (1)
        do l=1, n_profile
          if (x_map(l) <= n_x+horiz_limit .and. &
              y_map(l) > n_y-horiz_limit .and. &
              y_map(l) <= n_y*2+horiz_limit) then
            x_panel_map(l) = x_map(l)
            y_panel_map(l) = y_map(l) - n_y
          else if (x_map(l) > n_x*4-horiz_limit) then
            ! Wrap last equatorial panel into halo
            x_panel_map(l) = x_map(l) - n_x*4
            y_panel_map(l) = y_map(l) - n_y
          end if
        end do
      case (2)
        do l=1, n_profile
          if (x_map(l) > n_x-horiz_limit .and. &
              x_map(l) <= n_x*2+horiz_limit .and. &
              y_map(l) > n_y .and. y_map(l) <= n_y*2) then
            x_panel_map(l) = x_map(l) - n_x
            y_panel_map(l) = y_map(l) - n_y
          else if (x_map(l) > n_x-horiz_limit .and. &
                   y_map(l) > n_y*2) then
            ! Wrap north polar panel into halo
            x_panel_map(l) = y_map(l) - n_y*2
            y_panel_map(l) = n_y + n_x + 1 - x_map(l)
          else if (x_map(l) > n_x-horiz_limit .and. &
                   y_map(l) <= n_y) then
            ! Wrap south polar panel into halo
            x_panel_map(l) = n_y + 1 - y_map(l)
            y_panel_map(l) = x_map(l) - n_x
          end if
        end do
      case (3)
        do l=1, n_profile
          if (x_map(l) > n_x*2-horiz_limit .and. &
              x_map(l) <= n_x*3+horiz_limit) then
            x_panel_map(l) = x_map(l) - n_x*2
            y_panel_map(l) = y_map(l) - n_y
          else if (y_map(l) > n_y*3-horiz_limit) then
            ! Wrap north polar panel into halo
            x_panel_map(l) = n_x + 1 - x_map(l)
            y_panel_map(l) = n_y*4 + 1 - y_map(l)
          else if (y_map(l) <= horiz_limit) then
            ! Wrap south polar panel into halo
            x_panel_map(l) = n_x + 1 - x_map(l)
            y_panel_map(l) = 1 - y_map(l)
          end if
        end do
      case (4)
        do l=1, n_profile
          if (x_map(l) > n_x*3-horiz_limit) then
            x_panel_map(l) = x_map(l) - n_x*3
            y_panel_map(l) = y_map(l) - n_y
          else if (x_map(l) <= horiz_limit .and. &
                   y_map(l) > n_y .and. y_map(l) <= n_y*2) then
            ! Wrap first equatorial panel into halo
            x_panel_map(l) = x_map(l) + n_x
            y_panel_map(l) = y_map(l) - n_y
          else if (x_map(l) <= horiz_limit .and. &
                   y_map(l) > n_y*2) then
            ! Wrap north polar panel into halo
            x_panel_map(l) = n_y*3 + 1 - y_map(l)
            y_panel_map(l) = n_y + x_map(l)
          else if (x_map(l) <= horiz_limit .and. &
                   y_map(l) <= n_y) then
            ! Wrap south polar panel into halo
            x_panel_map(l) = y_map(l)
            y_panel_map(l) = 1 - x_map(l)
          end if
        end do
      case (5)
        ! South pole
        do l=1, n_profile
          if (y_map(l) <= n_y+horiz_limit) then
            if (x_map(l) <= n_x) then
              x_panel_map(l) = x_map(l)
              y_panel_map(l) = y_map(l)
            ! Wrap equatorial points into halo
            else if (x_map(l) > n_x .and. &
                     x_map(l) <= n_x*2) then
              x_panel_map(l) = y_map(l)
              y_panel_map(l) = n_x*2 + 1 - x_map(l)
            else if (x_map(l) > n_x*2 .and. &
                     x_map(l) <= n_x*3) then
              x_panel_map(l) = n_x*3 + 1 - x_map(l)
              y_panel_map(l) = n_y + 1 - y_map(l)
            else if (x_map(l) > n_x*3) then
              x_panel_map(l) = n_y + 1 - y_map(l)
              y_panel_map(l) = x_map(l) - n_x*3
            end if
          end if
        end do
      case (6)
        ! North pole
        do l=1, n_profile
          if (y_map(l) > n_y*2-horiz_limit) then
            if (x_map(l) <= n_x) then
              x_panel_map(l) = x_map(l)
              y_panel_map(l) = y_map(l) - n_y*2
            ! Wrap equatorial points into halo
            else if (x_map(l) > n_x .and. &
                     x_map(l) <= n_x*2) then
              x_panel_map(l) = n_x + n_y*2 + 1 - y_map(l)
              y_panel_map(l) = x_map(l) - n_x
            else if (x_map(l) > n_x*2 .and. &
                     x_map(l) <= n_x*3) then
              x_panel_map(l) = n_x*3 + 1 - x_map(l)
              y_panel_map(l) = n_y*3 + 1 - y_map(l)
            else if (x_map(l) > n_x*3) then
              x_panel_map(l) = y_map(l) - n_y*2
              y_panel_map(l) = n_x*4 + 1 - x_map(l)
            end if
          end if
        end do
      end select
    end if

    ! Map surface altitude, latitude and longitude, to 2D grid
    do l=1, n_profile
      if (x_panel_map(l) /= x_map_missing) then
        surface_altitude_2d(x_panel_map(l), y_panel_map(l)) &
          = surface_altitude(l)
        longitude_2d(x_panel_map(l), y_panel_map(l)) &
          = longitude(l) * pi / 180.0_RealK
        latitude_2d(x_panel_map(l), y_panel_map(l)) &
          = latitude(l) * pi / 180.0_RealK
      end if
    end do

    if (l_output_panel_fields .and. .not.l_local) then
      ! Output diagnostic fields of panel grids before stretching halos
      call nf(nf90_def_var(ncid, 'surface_altitude_halo_'//panel_num(panel), &
        NF90_DOUBLE, (/ dimid6, dimid7 /), varid))
      call nf(nf90_put_att(ncid, varid, 'standard_name', &
                                 'surface_altitude_halo_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'long_name', &
                                 'surface_altitude_halo_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'units', 'm'))
      call nf(nf90_enddef(ncid))
      call nf(nf90_put_var(ncid, varid, surface_altitude_2d(:,:)))
      call nf(nf90_redef(ncid))

      call nf(nf90_def_var(ncid, 'longitude_halo_'//panel_num(panel), &
        NF90_DOUBLE, (/ dimid6, dimid7 /), varid))
      call nf(nf90_put_att(ncid, varid, 'standard_name', &
                                 'longitude_halo_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'long_name', &
                                 'longitude_halo_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'units', 'degrees_east'))
      call nf(nf90_enddef(ncid))
      call nf(nf90_put_var(ncid, varid, longitude_2d(:,:)))
      call nf(nf90_redef(ncid))

      call nf(nf90_def_var(ncid, 'latitude_halo_'//panel_num(panel), &
        NF90_DOUBLE, (/ dimid6, dimid7 /), varid))
      call nf(nf90_put_att(ncid, varid, 'standard_name', &
                                 'latitude_halo_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'long_name', &
                                 'latitude_halo_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'units', 'degrees_north'))
      call nf(nf90_enddef(ncid))
      call nf(nf90_put_var(ncid, varid, latitude_2d(:,:)))
      call nf(nf90_redef(ncid))
    end if

    if (.not.l_local) then
      ! For the cubed-sphere, the panel halos need to be stretched into
      ! the corners so that x and y directions are along great circles.
      do j=1, horiz_limit
        ! Corners points bisect halo edge points
        surface_altitude_2d(n_edge+j, n_edge+j) &
          = 0.5_RealK*surface_altitude_2d(n_edge, n_edge+j) &
          + 0.5_RealK*surface_altitude_2d(n_edge+j, n_edge)
        surface_altitude_2d(n_edge+j, 1-j) &
          = 0.5_RealK*surface_altitude_2d(n_edge, 1-j) &
          + 0.5_RealK*surface_altitude_2d(n_edge+j, 1)
        surface_altitude_2d(1-j, 1-j) &
          = 0.5_RealK*surface_altitude_2d(1, 1-j) &
          + 0.5_RealK*surface_altitude_2d(1-j, 1)
        surface_altitude_2d(1-j, n_edge+j) &
          = 0.5_RealK*surface_altitude_2d(1-j, n_edge) &
          + 0.5_RealK*surface_altitude_2d(1, n_edge+j)
        call interp_latlon( 0.5_RealK, &
          latitude_2d(n_edge, n_edge+j), longitude_2d(n_edge, n_edge+j), &
          latitude_2d(n_edge+j, n_edge), longitude_2d(n_edge+j, n_edge), &
          latitude_2d(n_edge+j, n_edge+j), longitude_2d(n_edge+j, n_edge+j) )
        call interp_latlon( 0.5_RealK, &
          latitude_2d(n_edge, 1-j), longitude_2d(n_edge, 1-j), &
          latitude_2d(n_edge+j, 1), longitude_2d(n_edge+j, 1), &
          latitude_2d(n_edge+j, 1-j), longitude_2d(n_edge+j, 1-j) )
        call interp_latlon( 0.5_RealK, &
          latitude_2d(1, 1-j), longitude_2d(1, 1-j), &
          latitude_2d(1-j, 1), longitude_2d(1-j, 1), &
          latitude_2d(1-j, 1-j), longitude_2d(1-j, 1-j) )
        call interp_latlon( 0.5_RealK, &
          latitude_2d(1-j, n_edge), longitude_2d(1-j, n_edge), &
          latitude_2d(1, n_edge+j), longitude_2d(1, n_edge+j), &
          latitude_2d(1-j, n_edge+j), longitude_2d(1-j, n_edge+j) )
        ! Other halo points are interpolated
        tan_j = tan(eqt_res*(real(j, RealK)-0.5_RealK) + pi_over_4)
        allocate(surface_altitude_work(1-j:n_edge+j, 4))
        allocate(latitude_work(1-j:n_edge+j, 4))
        allocate(longitude_work(1-j:n_edge+j, 4))
        do i=2-j, n_edge+j-1
          tan_i = tan(eqt_res*(real(i, RealK)-0.5_RealK) - pi_over_4)
          halo_i = (atan(tan_i / tan_j) + pi_over_4) / eqt_res + 0.5_RealK
          weight_upper = halo_i - aint(halo_i)
          weight_lower = 1.0_RealK - weight_upper
          lower = int(halo_i)
          surface_altitude_work(i, 1) &
            = weight_lower*surface_altitude_2d(lower, n_edge+j) &
            + weight_upper*surface_altitude_2d(lower+1, n_edge+j)
          surface_altitude_work(i, 2) &
            = weight_lower*surface_altitude_2d(lower, 1-j) &
            + weight_upper*surface_altitude_2d(lower+1, 1-j)
          surface_altitude_work(i, 3) &
            = weight_lower*surface_altitude_2d(n_edge+j, lower) &
            + weight_upper*surface_altitude_2d(n_edge+j, lower+1)
          surface_altitude_work(i, 4) &
            = weight_lower*surface_altitude_2d(1-j, lower) &
            + weight_upper*surface_altitude_2d(1-j, lower+1)
          call interp_latlon( weight_upper, &
            latitude_2d(lower+1, n_edge+j), longitude_2d(lower+1, n_edge+j), &
            latitude_2d(lower, n_edge+j), longitude_2d(lower, n_edge+j), &
            latitude_work(i, 1), longitude_work(i, 1) )
          call interp_latlon( weight_upper, &
            latitude_2d(lower+1, 1-j), longitude_2d(lower+1, 1-j), &
            latitude_2d(lower, 1-j), longitude_2d(lower, 1-j), &
            latitude_work(i, 2), longitude_work(i, 2) )
          call interp_latlon( weight_upper, &
            latitude_2d(n_edge+j, lower+1), longitude_2d(n_edge+j, lower+1), &
            latitude_2d(n_edge+j, lower), longitude_2d(n_edge+j, lower), &
            latitude_work(i, 3), longitude_work(i, 3) )
          call interp_latlon( weight_upper, &
            latitude_2d(1-j, lower+1), longitude_2d(1-j, lower+1), &
            latitude_2d(1-j, lower), longitude_2d(1-j, lower), &
            latitude_work(i, 4), longitude_work(i, 4) )
        end do
        do i=2-j, n_edge+j-1
          surface_altitude_2d(i, n_edge+j) = surface_altitude_work(i, 1)
          surface_altitude_2d(i, 1-j) = surface_altitude_work(i, 2)
          surface_altitude_2d(n_edge+j, i) = surface_altitude_work(i, 3)
          surface_altitude_2d(1-j, i) = surface_altitude_work(i, 4)
          latitude_2d(i, n_edge+j) = latitude_work(i, 1)
          latitude_2d(i, 1-j) = latitude_work(i, 2)
          latitude_2d(n_edge+j, i) = latitude_work(i, 3)
          latitude_2d(1-j, i) = latitude_work(i, 4)
          longitude_2d(i, n_edge+j) = longitude_work(i, 1)
          longitude_2d(i, 1-j) = longitude_work(i, 2)
          longitude_2d(n_edge+j, i) = longitude_work(i, 3)
          longitude_2d(1-j, i) = longitude_work(i, 4)
        end do
        deallocate(longitude_work)
        deallocate(latitude_work)
        deallocate(surface_altitude_work)
      end do
    end if

    ! Calculate horizon angles and aspects
    call horizon_2d(n_x, n_y, n_horiz_layer, n_horiz_ang, horiz_limit, &
      latitude_2d, longitude_2d, surface_altitude_2d, &
      horizon_angle_2d, horizon_aspect_2d)

    if (l_local) then
      ! Correct horizon aspects at boundary that have been calculated
      ! using zero lat/lon.
      do k = 1, n_horiz_ang
        do y=3, n_y-2
          horizon_aspect_2d(1, y, k) = horizon_aspect_2d(3, y, k)
          horizon_aspect_2d(2, y, k) = horizon_aspect_2d(3, y, k)
          horizon_aspect_2d(n_x-1, y, k) = horizon_aspect_2d(n_x-2, y, k)
          horizon_aspect_2d(n_x, y, k)   = horizon_aspect_2d(n_x-2, y, k)
        end do
        do x=1, n_x
          horizon_aspect_2d(x, 1, k) = horizon_aspect_2d(x, 3, k)
          horizon_aspect_2d(x, 2, k) = horizon_aspect_2d(x, 3, k)
          horizon_aspect_2d(x, n_y-1, k) = horizon_aspect_2d(x, n_y-2, k)
          horizon_aspect_2d(x, n_y, k)   = horizon_aspect_2d(x, n_y-2, k)
        end do
      end do
    end if

    ! Map 2D horizon angles and aspects back to original mesh
    do k = 1, n_horiz_ang
      do l=1, n_profile
        if (x_panel_map(l) > 0 .and. x_panel_map(l) <= n_x .and. &
            y_panel_map(l) > 0 .and. y_panel_map(l) <= n_y) then
          horizon_angle(l, k) &
            = horizon_angle_2d(x_panel_map(l), y_panel_map(l), 1, k)
          horizon_aspect(l, k) &
            = horizon_aspect_2d(x_panel_map(l), y_panel_map(l), k)
        end if
      end do
    end do

    if (l_output_panel_fields) then
      ! Output diagnostic fields of panel grids
      call nf(nf90_def_var(ncid, 'surface_altitude_2d_'//panel_num(panel), &
        NF90_DOUBLE, (/ dimid6, dimid7 /), varid))
      call nf(nf90_put_att(ncid, varid, 'standard_name', &
                                 'surface_altitude_2d_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'long_name', &
                                 'surface_altitude_2d_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'units', 'm'))
      call nf(nf90_enddef(ncid))
      call nf(nf90_put_var(ncid, varid, surface_altitude_2d(:,:)))
      call nf(nf90_redef(ncid))

      call nf(nf90_def_var(ncid, 'longitude_2d_'//panel_num(panel), &
        NF90_DOUBLE, (/ dimid6, dimid7 /), varid))
      call nf(nf90_put_att(ncid, varid, 'standard_name', &
                                 'longitude_2d_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'long_name', &
                                 'longitude_2d_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'units', 'degrees_east'))
      call nf(nf90_enddef(ncid))
      call nf(nf90_put_var(ncid, varid, longitude_2d(:,:)))
      call nf(nf90_redef(ncid))

      call nf(nf90_def_var(ncid, 'latitude_2d_'//panel_num(panel), &
        NF90_DOUBLE, (/ dimid6, dimid7 /), varid))
      call nf(nf90_put_att(ncid, varid, 'standard_name', &
                                 'latitude_2d_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'long_name', &
                                 'latitude_2d_'//panel_num(panel)))
      call nf(nf90_put_att(ncid, varid, 'units', 'degrees_north'))
      call nf(nf90_enddef(ncid))
      call nf(nf90_put_var(ncid, varid, latitude_2d(:,:)))
      call nf(nf90_redef(ncid))

      do k = 1, n_horiz_ang
        horizon_angle_str='horizon_angle_'//dir_num(k)//'_'//panel_num(panel)
        call nf(nf90_def_var(ncid, horizon_angle_str, NF90_DOUBLE, &
          (/ dimid4, dimid5 /), varid))
        call nf(nf90_put_att(ncid, varid, 'standard_name', horizon_angle_str))
        call nf(nf90_put_att(ncid, varid, 'long_name', horizon_angle_str))
        call nf(nf90_put_att(ncid, varid, 'units', 'radian'))
        call nf(nf90_enddef(ncid))
        call nf(nf90_put_var(ncid, varid, horizon_angle_2d(:, :, 1, k)))
        call nf(nf90_redef(ncid))

        horizon_aspect_str='horizon_aspect_'//dir_num(k)//'_'//panel_num(panel)
        call nf(nf90_def_var(ncid, horizon_aspect_str, NF90_DOUBLE, &
          (/ dimid4, dimid5 /), varid))
        call nf(nf90_put_att(ncid, varid, 'standard_name', horizon_aspect_str))
        call nf(nf90_put_att(ncid, varid, 'long_name', horizon_aspect_str))
        call nf(nf90_put_att(ncid, varid, 'units', 'radian'))
        call nf(nf90_enddef(ncid))
        call nf(nf90_put_var(ncid, varid, horizon_aspect_2d(:, :, k)))
        call nf(nf90_redef(ncid))
      end do
    end if

  end do panel_loop
  deallocate(y_panel_map)
  deallocate(x_panel_map)
  deallocate(y_map)
  deallocate(x_map)
  deallocate(longitude_2d)
  deallocate(latitude_2d)
  deallocate(surface_altitude_2d)
  deallocate(horizon_aspect_2d)
  deallocate(horizon_angle_2d)
  deallocate(longitude)
  deallocate(latitude)
  deallocate(surface_altitude)

  ! Add horizon angle fields to the orog ancillary file
  call nf(nf90_def_var(ncid, 'horizon_angle', NF90_DOUBLE, &
    (/ dimid1, dimid2 /), varid))
  call nf(nf90_put_att(ncid, varid, 'standard_name', 'horizon_angle'))
  call nf(nf90_put_att(ncid, varid, 'long_name', 'horizon_angle'))
  call nf(nf90_put_att(ncid, varid, 'units', 'radian'))
  call nf(nf90_put_att(ncid, varid, 'location', 'face'))
  call nf(nf90_put_att(ncid, varid, 'mesh', trim(mesh_name)))
  call nf(nf90_put_att(ncid, varid, 'online_operation', 'once'))
  call nf(nf90_put_att(ncid, varid, 'coordinates', &
    trim(mesh_name)//lon_name_ext//' '//trim(mesh_name)//lat_name_ext))
  call nf(nf90_enddef(ncid))
  call nf(nf90_put_var(ncid, varid, horizon_angle))
  call nf(nf90_redef(ncid))

  call nf(nf90_def_var(ncid, 'horizon_aspect', NF90_DOUBLE, &
    (/ dimid1, dimid3 /), varid))
  call nf(nf90_put_att(ncid, varid, 'standard_name', 'horizon_aspect'))
  call nf(nf90_put_att(ncid, varid, 'long_name', 'horizon_aspect'))
  call nf(nf90_put_att(ncid, varid, 'units', 'radian'))
  call nf(nf90_put_att(ncid, varid, 'location', 'face'))
  call nf(nf90_put_att(ncid, varid, 'mesh', trim(mesh_name)))
  call nf(nf90_put_att(ncid, varid, 'online_operation', 'once'))
  call nf(nf90_put_att(ncid, varid, 'coordinates', &
    trim(mesh_name)//lon_name_ext//' '//trim(mesh_name)//lat_name_ext))
  call nf(nf90_enddef(ncid))
  call nf(nf90_put_var(ncid, varid, horizon_aspect))
  call nf(nf90_close(ncid))

  deallocate(horizon_aspect)
  deallocate(horizon_angle)

  write (*, '(3(a,i0))') ' n_x: ', n_x, ' | n_y: ', n_y, ' | n_panel: ', n_panel

contains

  subroutine interp_latlon(weight_upper, &
    lat_upper, lon_upper, lat_lower, lon_lower, lat, lon)
    implicit none
    real(RealK), intent(in) :: weight_upper
    real(RealK), intent(in) :: lat_upper, lon_upper, lat_lower, lon_lower
    real(RealK), intent(out) :: lat, lon
    real(RealK) :: weight_lower, min_lon

    ! A straight interpolation is done for efficiency which should be
    ! sufficiently accurate away from the poles.
    ! (This could be extended to a great circle calculation if required.)
    weight_lower = 1.0_RealK - weight_upper
    lat = weight_upper*lat_upper + weight_lower*lat_lower
    if (abs(lon_upper - lon_lower) > pi) then
      min_lon = -pi
      if (lon_upper > lon_lower) then
        lon = weight_upper*lon_upper + weight_lower*(lon_lower+twopi)
        if (lon_lower >= 0.0_RealK) min_lon = 0.0_RealK
      else
        lon = weight_upper*(lon_upper+twopi) + weight_lower*lon_lower
        if (lon_upper >= 0.0_RealK) min_lon = 0.0_RealK
      end if
      lon = modulo(lon - min_lon, twopi) + min_lon
    else
      lon = weight_upper*lon_upper + weight_lower*lon_lower
    end if
  end subroutine interp_latlon

  subroutine nf(status)
    use netcdf
    implicit none
    integer, intent(in):: status
    if (status /= NF90_NOERR) then
       write(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
       stop 'Stopped while performing netCDF operation'
    end if
  end subroutine nf

end program calc_horizon_angles
