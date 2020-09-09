module parameter_module
integer, parameter                    ::            n    = 128
integer, parameter                    ::            niu  = 6
integer, parameter                    ::            ib   = 0
real(8), parameter                    ::            pi   = 3.1415926
real(8), parameter                    ::            rho  = 1.0
real(8), parameter                    ::            beta = 1.0
real(8), parameter                    ::            di   = 1.0
real(8), dimension(niu, ib:n-1)       ::            p, q
real(8)                               ::            dr, dk
real(8), dimension(ib:n-1)            ::            cr, ck, fmayer
real(8), dimension(ib:n-1)            ::            r, k, hk, sk
real(8), dimension(ib:n-1)            ::            gammar, gammak, pgammar
real(8), dimension(ib:n-1)            ::            dgamma, pdgamma
real(8), dimension(niu)               ::            aa, paa, da
real(8), dimension(niu, niu)          ::            jac, ijac
real(8)                               ::            ksi, gamma_1, gamma_2, delta_2
