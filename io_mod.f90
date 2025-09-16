Module io_mod
    use netcdf
    use base
    contains
        subroutine read_bd(sst, divm, uam, vam, uom, vom, wem)
            Implicit None
            Real, dimension(nx_o, ny_o) :: sst, divm, uam, vam, uom, vom, wem
            Integer :: ncid, varid
            call check(nf90_open('vars_bd.nc', nf90_nowrite, ncid))
            call check(nf90_inq_varid(ncid, 'sst', varid))
            call check(nf90_get_var(ncid, varid, sst, start=[1, 1, 1], count=[nx_o, ny_o, 1]))

        end subroutine read_bd

        subroutine read_vars(TO)
            Implicit None
            Real, dimension(nx_o, ny_o) :: TO
            Integer :: ncid, varid

            call check(nf90_open('vars_init.nc', nf90_nowrite, ncid))
            call check(nf90_inq_varid(ncid, 'TO', varid))
            call check(nf90_get_var(ncid, varid, TO))
            end subroutine read_vars

        subroutine check(status)
            Integer :: status
            if(status /= 0)then
                print*, nf90_strerror(status)
                stop
            endif
        end subroutine check

End Module io_mod
