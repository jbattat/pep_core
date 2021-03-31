      subroutine NDSRAD(ncall,solcon,x,rmag,bcor,rrsmag,radcon,
     .                  shape, fi, qrad, q)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      a3, a5, a8, aa, aas, ab, abs, ac, acs, ad, 
     .          ads, ae, aes, af, afs, ags, as
      real      fd1, fd2, fd3, fd4, fd5, fd6c, fd6e, fd7, fd8, fda, fdb
      real      fdc, fdd, fde, fdf, fn1, fn2, fn3, fn4, fn5, fn6c, fn6e,
     .          fn7, fn8, fna, fnb, fnc, fnd, fne, fnf, fs1
      real      fs2, fs3, fs4, fs5, fs6c, fs6e, fs7, fs8, fsa, fsb, fsc,
     .          fsd, fse, fsf, fx1, fx2, fx3, fx4, fx5, fx6
      real      fx7, fx8, fx9, fxa, fxb, fxc, fxd, fxe, fxf, fz1, fz2, 
     .          fz3, fz4, fz5, fz6, fz7, fz8, fz9, fza, fzb
      real      fzc, fzd, fze, fzf, pa8, paa, pab, pac, pad, pae, paf, 
     .          pag, pas, pd1, pd2, pd3, pd4, pd5, pd6c, pd6e
      real      pd7, pd8, pda, pdb, pdc, pdd, pde, pdf, pl, pn1, pn2, 
     .          pn3, pn4, pn5, pn6c, pn6e, pn7, pn8, pna, pnb
      real      pnc, pnd, pne, pnf, ps1, ps2, ps3, ps4, ps5, ps6c, ps6e,
     .          ps7, ps8, psa, psb, psc, psd, pse, psf, psxa
      real      psxb, psxc, psxd, psxe, psxf, psxt, pu, px1, px2, px3, 
     .          px4, px5, px6, px7, px8, px9, pxa, pxb, pxc, pxd
      real      pxe, pxf, pz1, pz2, pz3, pz4, pz5, pz6, pz7, pz8, pz9, 
     .          pza, pzb, pzc, pzd, pze, pzf, sla, sxa, sxb
      real      sxc, sxd, sxe, sxf, sxt, t10, t11, t12, t4, 
     .          t45, t6, t7, t8mt11, t8mt6, t8mt9
      real      t9, u, vz1gz1
      integer   i, j, k, l, lmk, lxk, m, ncall
 
c*** end of declarations inserted by spag
 
 
c     r.king jun 1982   pep version of gps radiation model nds
c     coded from celeste version obtained from j.carr of nswc
c
c         ncall= -1  calculate forces for equations of motion
c              = -2  calculate partial derivatives
c         solcon = pep solar radiation pressure constant for gps
c                  satellites (assumes cross section of 6.0088 sq m)
c         x(6) =  sbcor(6) = s/c w.r.t. earth (au)
c         rmag= rsb = magnitude of x
c         bcor(6) = s/c w.r.t. sun (au)
c         rsmag= rsb = magnitude of bcor
c         radcon(3)= con(9), con(10), con(11) in ertorb
c            radcon(1)= oi39 = coefficient of radiation presssure in
c                              direction perpendicular to solar panels
c            radcon(2)= oi40 = coef. of rad. press. in y-axis direction
c            radcon(3)= oi41 = angle between panel axis and y-axis (rad)
c         shape= lambda = fraction of sun not occulted by earth
c         fi(3) =  output force along inertial x,y,z axes
c         qrad(3,3)= explicit partials of fi w.r.t. radcon
c         q(3,3)= partials of fi w.r.t. coordinates ( x )
c
      real*10 solcon, x(6),rmag,bcor(6),rrsmag,radcon(3),shape,
     .          fi(3),cc,rs(3),rsmag2,rrs(3),urrs(3),ur(3),ux(3),
     .          uy(3),uz(3),pak(3,3),fr(3),fx,fy,fz,uxmag,
     .          uymag, xys, xx, srm, DOT, oi39, oi40, oi41, vz1, vz12, 
     .          wz1, wz12, gz1, vz1wz1, dz1, xiz1, hlk1(3,3),hlk2(3,
     .          3), hlk3(3,3),x2,xy,yy,xpy,tmj(3,3),xpys,xpyss,
     .          ybmk3(3,3),mk2zb(3,3),cckr1,elk(3,3),eilk(3,3),
     .          h(3),elks,p(3,3),qrad(3,3),q(3,3),px,py,pz,d,
     .          d1
 
      equivalence (tmj(1,1),ux),(tmj(1,2),uy),
     .            (tmj(1,3),uz)
c
c constants of nds model
      real*4 r1/.7686/, r2/.4700/, r3/.0317/, r4/.0752/, r42/.1504/,
     .     t2/.6066/, t44/.0329/, t5/.8739/, t8a/.1905/, t8e/.2583/,
     .     b2/-.8972/, b8/.1351/, b9/.2758/, b10/.4431/, as1/.3221/,
     .     as2/.2807/, a33/.8350/, a8ax/1.3124/, sn1/1.3204/,
     .     ss1/.6266/, sd1/.0631/, sn2/1.1421/, ss2/.3337/, sd2/.1058/,
     .     sn3/1.6880/, ss3/.3120/, sd3/.1147/, sn4/.6940/, ss4/.6940/,
     .     sd4/0./, sn5/1.0143/, ss5/.9570/, sd5/.4278/, sn6c/0./,
     .     ss6c/0./, sd6c/0.0/, sn6e/.0050/, ss6e/.0013/, sd6e/.0002/,
     .     sn6/1.215/, ss6/.3550/, sd6/.1126/, sn7/.0283/, ss7/.0072/,
     .     sd7/.0012/, sn8/1.6450/, ss8/.3550/, sd8/.1433/,
     .     ss9/-7.1276/, s2ps3x/.6587/, s4ms5x/.1848/, s6ms7x/.0653/,
     .     b4mr4x/.6206/, b5mr4x/.4803/, b6mr4x/.7880/, b7mr3x/.8099/,
     .     va/.2778/, vb/.7178/, vc/.2934/, vd/.3687/, ve/.9695/,
     .     vf/.3115/, vt/.1599/, slt/0./
c
c
      if(ncall.lt.-1) then
c
c---------------------------------------------------------------
c
c calculate partial derivatives
         px = 0._10
         py = 0._10
         pz = 0._10
c
c main body side
c
         pn1 = -2.*vz1*sn1
         ps1 = ss1*(wz1 - vz1gz1)
         pd1 = -sd1*gz1
         px1 = -(pn1 + pd1)
         pz1 = -ps1
         px  = px + px1
         pz  = pz + pz1
c
c apogee engine side
c
         if(sla.gt.0.) then
            pn2 = sn2*(wz1*pl - gz1*sla)
            ps2 = ss2*(sla + vz1*pl)
            pd2 = sd2*pl
            px2 = -(pn2 + pd2)
            pz2 = -ps2
            px  = px + px2
            pz  = pz + pz2
         endif
c
c main body aft end
c
         if(vz1.lt.0.) then
            ps3 = ss3*((wz1-vz1gz1)*a3 - pas*vz1wz1)
            pd3 = sd3*(a3 - pas*vz1)
            pn3 = sn3*(2.*vz1*a3 - pas*vz12)
            px3 = ps3
            pz3 = pn3 - pd3
            px  = px + px3
            pz  = pz + pz3
         endif
c
c apogee engine aft end
c
         if(vz1.lt.0.) then
            pn4 = 2.*sn4*vz1
            ps4 = ss4*(wz1 - vz1gz1)
            pd4 = sd4
 
            px4 = ps4
            pz4 = pn4 - pd4
            px  = px + px4
            pz  = pz + pz4
         endif
c
c tt and c antenna side
c
         if(t45.gt.0.) then
            pn5 = sn5*(wz12*pag - 2.*vz1*a5)
            ps5 = ss5*((wz1-vz1gz1)*a5 + vz1wz1*pag)
            pd5 = sd5*(wz1*pag - gz1*a5)
            px5 = -(pn5 + pd5)
            pz5 = -ps5
            px  = px + px5
            pz  = pz + pz5
         endif
c
c tt and c antenna end
c
         if(vz1.gt.0.) then
            pn6e = 2.*vz1*sn6e
            ps6e = ss6e*(wz1 - vz1gz1)
            pd6e = sd6e
         endif
         if(t45.ge.slt) then
            pn6c = -2.*vz1*sn6c
            ps6c = ss6c*(wz1 - vz1gz1)
            pd6c = -sd6c*gz1
         endif
         px6 = -(ps6e + pn6c + pd6c)
         pz6 = -(pn6e + pd6e + ps6c)
         px  = px + px6
         pz  = pz + pz6
c
c navigation antenna pair a - 6
c
         if(t8mt6.gt.0.) then
            pna = sn6*(wz12*paa - 2.*vz1*aa)
            psa = ss6*(vz1wz1*paa + aa*(wz1-vz1gz1))
            pda = sd6*(wz1*paa - aa*gz1)
            pxa = -2.*(pna + pda)
            pza = -2.*psa
            px  = px + pxa
            pz  = pz + pza
         endif
c
c navigation antenna pair d - 7
c
         if(vz1.gt.0.) then
            pnd = sn6*(wz12*pad - 2.*vz1*ad)
            psd = ss6*(vz1wz1*pad + ad*(wz1-vz1gz1))
            pdd = sd6*(wz1*pad - ad*gz1)
            pxd = -2.*(pnd + pdd)
            pzd = -2.*psd
            px  = px + pxd
            pz  = pz + pzd
         endif
c
c antenna pair b - 9
c
         if(t8mt9.gt.0.) then
            pnb = sn6*(wz12*pab - 2.*vz1*ab)
            psb = ss6*(vz1wz1*pab + ab*(wz1-vz1gz1))
            pdb = sd6*(wz1*pab - ab*gz1)
            pxb = -2.*(pnb + pdb)
            pzb = -2.*psb
            px  = px + pxb
            pz  = pz + pzb
         endif
c
c antenna pair c - 10
c
         if(vz1.gt.0.) then
            pnc = sn6*(wz12*pac - 2.*vz1*ac)
            psc = ss6*(vz1wz1*pac + ac*(wz1-vz1gz1))
            pdc = sd6*(wz1*pac - ac*gz1)
            pxc = -2.*(pnc + pdc)
            pzc = -2.*psc
            px  = px + pxc
            pz  = pz + pzc
         endif
c
c antenna pair e - 10
c
         if(t8mt11.gt.0.) then
            pne = sn6*(wz12*pae - 2.*vz1*ae)
            pse = ss6*(vz1wz1*pae + ae*(wz1-vz1gz1))
            pde = sd6*(wz1*pae - ae*gz1)
            pxe = -2.*(pne + pde)
            pze = -2.*pse
            px  = px + pxe
            pz  = pz + pze
         endif
c
c antenna pair f - 12
c
         if(vz1.gt.0.) then
            pnf = sn6*(wz12*paf - 2.*vz1*af)
            psf = ss6*(vz1wz1*paf + af*(wz1-vz1gz1))
            pdf = sd6*(wz1*paf - af*gz1)
            pxf = -2.*(pnf + pdf)
            pzf = -2.*psf
            px  = px + pxf
            pz  = pz + pzf
         endif
c
c forces on antenna end
c
         if(vz1.ge.0.) then
            pn7 = 2.*vz1*sn7
            ps7 = ss7*(wz1 - vz1gz1)
            pd7 = sd7
            px7 = -12.*ps7
            pz7 = -12.*(pn7 + pd7)
            px  = px + px7
            pz  = pz + pz7
         endif
c
c main body forward end - 6
c
         if(vz1.gt.0) then
            pn8 = sn8*(vz12*pa8 + 2.*vz1*a8)
            ps8 = ss8*(vz1wz1*pa8 + a8*(wz1-vz1gz1))
            pd8 = sd8*(pa8*vz1 + a8)
            px8 = -ps8
            pz8 = -(pn8 + pd8)
            px  = px + px8
            pz  = pz + pz8
         endif
c
c solar arrays  9
c
         px9 = ss9*(-gz1)
         pz9 = ss9
         px  = px + px9
         pz  = pz + pz9
c
c compute total derivatives
         do l = 1, 3
            d = 0._10
            do k = 1, 3
               if(l.eq.k) d = 1._10
               elk(l,k)  = (ur(l)*ur(k) - d)/rmag
               eilk(l,k) = (d - urrs(l)*urrs(k))/rrsmag
            end do
         end do
         do k = 1, 3
            srm = 0._10
            do l = 1, 3
               srm = srm + ur(l)*eilk(l,k) + urrs(l)*elk(l,k)
            end do
            h(k) = srm
         end do
 
         do k = 1, 3
            p(1,k) = px*h(k)
            p(2,k) = py*h(k)
            p(3,k) = pz*h(k)
         end do
 
         do l = 1, 3
            do k = 1, 3
               d1  = 0._10
               lmk = iabs(l - k)
               if(lmk.eq.0) then
                  elks = 0._10
                  lxk  = 1
               else if(lmk.eq.2) then
                  elks = (k - l)/2.
                  lxk  = 6 - l - k
               else
                  elks = l - k
                  lxk  = 6 - l - k
               endif
               if(l.eq.k) d1 = 1._10
               hlk3(l,k) = (ur(l)*ur(k) - d1)/rmag
               if(uymag.ne.0.) hlk2(l,k)
     .             = (elks*rrs(lxk) - tmj(l,2)
     .             *(x(k)*rsmag2-rs(k)*(x(1)*rs(1)+x(2)*rs(2)+x(3)
     .             *rs(3)))/uymag)/uymag
            end do
         end do
 
         if(uymag.eq.0.) then
            x2    = x(1)*x(1)
            xy    = x(1)*x(2)
            yy    = x(2)*x(2)
            xpy   = x2 + yy
            xpys  = SQRT(xpy)
            xpyss = xpy*xpys
            hlk2(1,1) = x(1)*x(2)/xpyss
            hlk2(1,2) = -x2/xpyss
            hlk2(1,3) = 0._10
            hlk2(2,1) = yy/xpyss
            hlk2(2,2) = -xy/xpyss
            hlk2(2,3) = 0._10
            hlk2(3,1) = 0._10
            hlk2(3,2) = 0._10
            hlk2(3,3) = 0._10
         endif
 
         do k = 1, 3
            ybmk3(k,1) = uy(2)*hlk3(3,k) - uy(3)*hlk3(2,k)
            ybmk3(k,2) = uy(3)*hlk3(1,k) - uy(1)*hlk3(3,k)
            ybmk3(k,3) = uy(1)*hlk3(2,k) - uy(2)*hlk3(1,k)
            mk2zb(k,1) = hlk2(2,k)*ur(3) - hlk2(3,k)*ur(2)
            mk2zb(k,2) = hlk2(3,k)*ur(1) - hlk2(1,k)*ur(3)
            mk2zb(k,3) = hlk2(1,k)*ur(2) - hlk2(2,k)*ur(1)
         end do
 
         do k = 1, 3
            do l = 1, 3
               hlk1(l,k) = ybmk3(k,l) - mk2zb(k,l)
            end do
         end do
         cckr1 = cc*oi39
         do m = 1, 3
            do k = 1, 3
               q(m,k) = cckr1*(tmj(m,1)*p(1,k) + tmj(m,2)*p(2,k)
     .                   + tmj(m,3)*p(3,k)) + hlk1(m,k)*fr(1)
     .                   + hlk2(m,k)*fr(2) + hlk3(m,k)*fr(3)
            end do
         end do
         do j = 1, 3
            do i = 1, 3
               qrad(i,j) = tmj(i,1)*pak(1,j) + tmj(i,2)*pak(2,j)
     .                      + tmj(i,3)*pak(3,j)
            end do
         end do
      else
c
c
c determine radiation constant for nds model
         cc   = solcon/6.0088_10
         oi39 = radcon(1)
         oi40 = radcon(2)
         oi41 = radcon(3)
c
c initialize forces along body axes
         fx = 0._10
         fy = 0._10
         fz = 0._10
c
c setup sun/satellite/earth vectors
         do i = 1, 3
            rrs(i)  = -bcor(i)
            rs(i)   = rrs(i) + x(i)
            ur(i)   = x(i)/rmag
            urrs(i) = rrs(i)/rrsmag
         end do
         rsmag2 = DOT(rs,rs)
c
c calculate geometric quantities needed by nds model
         vz1    = -DOT(ur,urrs)
         vz12   = vz1*vz1
         wz1    = SQRT(1. - vz12)
         wz12   = wz1*wz1
         gz1    = vz1/wz1
         vz1gz1 = vz1*gz1
         vz1wz1 = vz1*wz1
         dz1    = -1./(vz12*wz1)
         xiz1   = 1./(wz12*wz1)
c
c main body side
c
         fn1 = sn1*wz12
         fs1 = ss1*vz1wz1
         fd1 = sd1*wz1
c
c forces
c
         fx1 = -(fn1 + fd1)
         fz1 = -fs1
         fx  = fx + fx1
         fz  = fz + fz1
c
c apogee engine side
c
         u  = 0.
         pu = 0.
         if(vz1.gt.0.) then
            pu = r1 - r2
            u  = pu*vz1
         endif
         sla = t2*wz1 - u
         pl  = -t2*gz1 - pu
         if(sla.gt.0) then
            fn2 = sn2*wz1*sla
            fs2 = ss2*vz1*sla
            fd2 = sd2*sla
            fx2 = -(fn2 + fd2)
            fz2 = -fs2
            fx  = fx + fx2
            fz  = fz + fz2
         endif
c
c main body aft end
c
         if(vz1.lt.0.) then
            pas = 0.
            if((0.gt.vz1) .and. (vz1.gt.b2)) then
               as = as1
            else if(vz1.eq.b2) then
               as = as2
            else if((b2.gt.vz1) .and. (vz1.gt.-1.)) then
               as  = -2.*r2*t2/gz1
               pas = -2.*r2*t2*dz1
            else if(vz1.eq.-1.) then
               as = 0.
            endif
            a3  = a33 - as
            fn3 = sn3*vz12*a3
            fs3 = ss3*vz1wz1*a3
            fd3 = sd3*vz1*a3
            fx3 = fs3
            fz3 = fn3 - fd3
            fx  = fx + fx3
            fz  = fz + fz3
         endif
c
c apogee engine aft end
c
         if(vz1.lt.0.) then
            fn4 = sn4*vz12
            fs4 = ss4*vz1wz1
            fd4 = sd4*vz1
            fx4 = fs4
            fz4 = fn4 - fd4
            fx  = fx + fx4
            fz  = fz + fz4
         endif
c
c tt and c antenna side
c
         t4  = 0.
         pag = 0.
         if(vz1.lt.0.) t4 = t44*gz1
         t45 = t5 - t4
         pag = t44*xiz1
         if(t45.gt.0.) then
            a5  = 2.*r3*t45
            pag = -2.*r3*pag
            fn5 = sn5*wz12*a5
            fs5 = ss5*vz1wz1*a5
            fd5 = sd5*wz1*a5
            fx5 = -(fn5 + fd5)
            fz5 = -fs5
            fx  = fx + fx5
            fz  = fz + fz5
         endif
c
c tt and c antenna end
c
         fn6e = 0.
         fd6e = 0.
         fs6e = 0.
         pn6e = 0.
         pd6e = 0.
         ps6e = 0.
         if(vz1.gt.0.) then
            fn6e = sn6e*vz12
            fs6e = ss6e*vz1wz1
            fd6e = sd6e*vz1
         endif
         t45  = t4 + t5
         fn6c = 0.
         fs6c = 0.
         fd6c = 0.
         pn6c = 0.
         ps6c = 0.
         pd6c = 0.
         if(t45.ge.slt) then
            fn6c = sn6c*wz12
            fs6c = ss6c*vz1wz1
            fd6c = sd6c*wz1
         endif
         fx6 = -(fs6e + fn6c + fd6c)
         fz6 = -(fn6e + fd6e + fs6c)
         fx  = fx + fx6
         fz  = fz + fz6
c
c navigation antenna pair(6a)
c
         t6  = 0.
         paa = 0.
         if(vz1.le.0.) then
            t6  = gz1*(r4 - b8)
            paa = xiz1*(b8 - r4)
         endif
         t8mt6 = t8a - t6
         if(t8mt6.gt.0.) then
            aa  = 2.*r4*(t8mt6)
            paa = 2.*r4*paa
            fna = sn6*aa*wz12
            fsa = ss6*aa*vz1wz1
            fda = sd6*aa*wz1
            fxa = -2.*(fna + fda)
            fza = -2.*fsa
            fx  = fx + fxa
            fz  = fz + fza
         endif
c
c antenna pair d - 7
c
         if(vz1.gt.0.) then
            t7  = 0.
            pad = 0.
            if(va.gt.vz1) then
               t7  = -s2ps3x*gz1 + t8a
               pad = s2ps3x*xiz1
            endif
            ad  = 2.*r4*(t8a - t7)
            pad = 2.*r4*pad
            fnd = sn6*ad*wz12
            fsd = ss6*ad*vz1wz1
            fdd = sd6*ad*wz1
            fxd = -2.*(fnd + fdd)
            fzd = -2.*fsd
            fx  = fx + fxd
            fz  = fz + fzd
         endif
c
c antenna pair b - 9
c
         t9  = 0.
         pab = 0.
         if(vz1.le.0) then
            t9  = (r4 - b9)*gz1
            pab = (b9 - r4)*xiz1
         endif
         t8mt9 = t8a - t9
         pab   = 2.*r4*pab
         if(t8mt9.gt.0.) then
            ab  = 2.*r4*(t8mt9)
            fnb = sn6*ab*wz12
            fsb = ss6*ab*vz1wz1
            fdb = sd6*ab*wz1
            fxb = -2.*(fnb + fdb)
            fzb = -2.*fsb
            fx  = fx + fxb
            fz  = fz + fzb
         endif
c
c antenna pair c - 10
c
         if(vz1.gt.0.) then
            t10 = 0.
            pac = 0.
            if(vb.gt.vz1) then
               t10 = -s4ms5x*gz1 + t8a
               pac = s4ms5x*xiz1
            endif
            ac  = 2.*r4*(t8a - t10)
            pac = 2.*r4*pac
            fnc = sn6*ac*wz12
            fsc = ss6*ac*vz1wz1
            fdc = sd6*ac*wz1
            fxc = -2.*(fnc + fdc)
            fzc = -2.*fsc
            fx  = fx + fxc
            fz  = fz + fzc
         endif
c
c antenna pair  e - 10
c
         t11 = 0.
         pae = 0.
         if(vz1.le.0.) then
            t11 = (r4 - b10)*gz1
            pae = (b10 - r4)*xiz1
         endif
         t8mt11 = t8e - t11
         if(t8mt11.gt.0.) then
            pae = 2.*r4*pae
            ae  = 2.*r4*(t8mt11)
            fne = sn6*ae*wz12
            fse = ss6*ae*vz1wz1
            fde = sd6*ae*wz1
            fxe = -2.*(fne + fde)
            fze = -2.*fse
            fx  = fx + fxe
            fz  = fz + fze
         endif
c
c antenna pair f - 12
c
         if(vz1.gt.0.) then
            t12 = 0.
            paf = 0.
            if(ve.gt.vz1) then
               t12 = -s6ms7x*gz1 + t8e
               paf = s6ms7x*xiz1
            endif
            af  = 2.*r4*(t8e - t12)
            paf = 2.*r4*paf
            fnf = sn6*af*wz12
            fsf = ss6*af*vz1wz1
            fdf = sd6*af*wz1
            fxf = -2.*(fnf + fdf)
            fzf = -2.*fsf
            fx  = fx + fxf
            fz  = fz + fzf
         endif
c
c force on antenna ends
c
         if(vz1.ge.0.) then
            fn7 = sn7*vz12
            fs7 = ss7*vz1wz1
            fd7 = sd7*vz1
            fx7 = -12.*fs7
            fz7 = -12.*(fn7 + fd7)
            fx  = fx + fx7
            fz  = fz + fz7
         endif
c
c main body forward end - 8
c
         if(vz1.gt.0.) then
            aas  = 0.
            abs  = 0.
            acs  = 0.
            ads  = 0.
            aes  = 0.
            afs  = 0.
            ags  = 0.
            psxa = 0.
            psxb = 0.
            psxc = 0.
            psxd = 0.
            psxe = 0.
            psxf = 0.
            psxt = 0.
            if(vz1.ne.1.) then
               if((1..gt.vz1) .and. (vz1.ge.va)) then
                  sxa  = t8a/gz1
                  psxa = t8a*dz1
               else if(vz1.lt.va) then
                  sxa = s2ps3x
               endif
               if((1..gt.vz1) .and. (vz1.ge.vb)) then
                  sxb  = t8a/gz1
                  psxb = t8a*dz1
               else if(vz1.lt.vb) then
                  sxb = s4ms5x
               endif
               if((1..gt.vz1) .and. (vz1.ge.vc)) then
                  sxc  = t8a/gz1
                  psxc = t8a*dz1
               else if(vz1.lt.vc) then
                  sxc = b4mr4x
               endif
               if((1..gt.vz1) .and. (vz1.ge.vd)) then
                  sxd  = t8a/gz1
                  psxd = t8a*dz1
               else if(vz1.lt.vd) then
                  sxd = b5mr4x
               endif
               if((1..gt.vz1) .and. (vz1.ge.ve)) then
                  sxe  = t8e/gz1
                  psxe = t8e*dz1
               else if(vz1.lt.ve) then
                  sxe = s6ms7x
               endif
               if((1..gt.vz1) .and. (vz1.ge.vf)) then
                  sxf  = t8e/gz1
                  psxf = t8e*dz1
               else if(vz1.lt.vf) then
                  sxf = b6mr4x
               endif
               if((1..gt.vz1) .and. (vz1.ge.vt)) then
                  sxt  = (t5 + slt)/gz1
                  psxt = (t5 + slt)*dz1
               else if(vz1.lt.vt) then
                  sxt = b7mr3x
               endif
               aas = r42*sxa
               abs = r42*sxb
               acs = r42*sxc
               ads = r42*sxd
               aes = r42*sxe
               afs = r42*sxf
               ags = 2.*r3*sxt
            endif
            a8  = a8ax - 2.*(aas + abs + acs + ads + aes + afs) - ags
            pa8 = -2.*r42*(psxa + psxb + psxc + psxd + psxe + psxf)
     .            - 2.*r3*psxt
            fn8 = sn8*vz12*a8
            fs8 = ss8*vz1wz1*a8
            fd8 = sd8*vz1*a8
            fx8 = -fs8
            fz8 = -(fn8 + fd8)
            fx  = fx + fx8
            fz  = fz + fz8
         endif
c
c solar arrays  9
c
         fx9 = ss9*wz1
         fz9 = ss9*vz1
         fx  = fx + fx9
         fz  = fz + fz9
c
c compute total forces along body axes
c
         fr(1)     = fx
         fr(2)     = fy
         fr(3)     = fz
         pak(1,1) = cc*fr(1)
         pak(2,1) = 0.
         pak(3,1) = cc*fr(3)
         pak(1,2) = COS(oi41)*shape*1E-12
         pak(2,2) = SIN(oi41)*shape*1E-12
c substitute the following two cards to remove the factor
c "shape*1e-12" from the y axis force model.
c pak(1,2) = COS(oi41)
c pak(2,2) = SIN(oi41)
         pak(3,2) = 0.
         pak(1,3) = -oi40*pak(2,2)
         pak(2,3) = oi40*pak(1,2)
         pak(3,3) = 0.
         fr(1)     = oi39*pak(1,1) + pak(2,3)
         fr(2)     = -pak(1,3)
         fr(3)     = oi39*pak(3,1)
c
c transform forces to inertial axes
         uy(1) = rrs(2)*x(3) - rrs(3)*x(2)
         uy(2) = rrs(3)*x(1) - rrs(1)*x(3)
         uy(3) = rrs(1)*x(2) - rrs(2)*x(1)
 
         uymag = SQRT(uy(1)*uy(1) + uy(2)*uy(2) + uy(3)*uy(3))
         if(uymag.eq.0.) then
            xx    = x(1)*x(1)
            yy    = x(2)*x(2)
            xys   = SQRT(xx + yy)
            uy(1) = -x(2)/xys
            uy(2) = x(1)/xys
            uy(3) = 0.
         else
 
            uy(1) = uy(1)/uymag
            uy(2) = uy(2)/uymag
 
            uy(3) = uy(3)/uymag
         endif
 
         ux(1) = uy(2)*ur(3) - uy(3)*ur(2)
         ux(2) = uy(3)*ur(1) - uy(1)*ur(3)
         ux(3) = uy(1)*ur(2) - uy(2)*ur(1)
 
         uxmag = SQRT(ux(1)*ux(1) + ux(2)*ux(2) + ux(3)*ux(3))
 
         ux(1) = -ux(1)/uxmag
         ux(2) = -ux(2)/uxmag
         ux(3) = -ux(3)/uxmag
 
         uz(1) = -ur(1)
         uz(2) = -ur(2)
         uz(3) = -ur(3)
 
         do i = 1, 3
            srm = 0._10
            do j = 1, 3
               srm = srm + tmj(i,j)*fr(j)
            end do
            fi(i) = srm
 
         end do
      endif
 
      return
      end
