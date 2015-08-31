import numpy as np


def transfnc_eh(k, Om=0.3, h=0.7, T_cmb=2.725, incl_baryons=True):
    """
    Compute the Eisenstein & Hu (1997) Transfer function for wave number k
    given a set of cosmological parameters.
    """

    if not incl_baryons:

        q = (k * (T_cmb / 2.7) ** 2) / (Om * h ** 2)

        Bc = 1. / (1 - 0.949 * fvb)
        av =

        Lq = np.log(np.exp(1) + 1.84 * Bc * np.sqrt(av) * q)
        Cq = 14.4 + 325. / (1 + 60.5 * q ** 1.11)

        return Lq / (Lq + Cq * q ** 2)

    elif incl_baryons:


        # TODO: Fix up this part

        # convert k to Mpc^-1 rather than hMpc^-1
        rk = k * h
        Omh2 = Om * h ** 2

        thet=2.728/2.7;

        # Eqn 4 - redshift of drag epoch
        b1=0.313*pow(om_mhsq,-0.419)*(1.+0.607*pow(om_mhsq,0.674));
        b2=0.238*pow(om_mhsq,0.223);
        zd=1291.*(1.+b1*pow(om_b*hsq,b2))*pow(om_mhsq,0.251)
            /(1.+0.659*pow(om_mhsq,0.828));

        // Equation 2 - redshift of matter-radiation equality
        ze=2.50e4*om_mhsq/thetpf;

        // value of R=(ratio of baryon-photon momentum density) at drag epoch
        rd=31500.*om_b*hsq/(thetpf*zd);

        // value of R=(ratio of baryon-photon momentum density) at epoch of matter-radiation equality
        re=31500.*om_b*hsq/(thetpf*ze);

        // Equation 3 - scale of ptcle horizon at matter-radiation equality
        rke=7.46e-2*om_mhsq/(thetsq);

        // Equation 6 - sound horizon at drag epoch
        s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1.+sqrt(re)));

        // Equation 7 - silk damping scale
        rks=1.6*pow(om_b*hsq,0.52)*pow(om_mhsq,0.73)*(1.+pow(10.4*om_mhsq,-0.95));

        // Equation 10  - define q
        q=rk/13.41/rke;

        / Equations 11 - CDM transfer function fits
        a1=pow(46.9*om_mhsq,0.670)*(1.+pow(32.1*om_mhsq,-0.532));
        a2=pow(12.0*om_mhsq,0.424)*(1.+pow(45.0*om_mhsq,-0.582));
        ac=pow(a1,(-om_b/om_m))*pow(a2,pow(-(om_b/om_m),3.));

        // Equations 12 - CDM transfer function fits
        b1=0.944/(1.+pow(458.*om_mhsq,-0.708));
        b2=pow(0.395*om_mhsq,-0.0266);
        bc=1./(1.+b1*(pow(1.-om_b/om_m,b2)-1.));

        // Equation 18
        f=1./(1.+pow(rk*s/5.4,4.));

        // Equation 20
        c1=14.2 + 386./(1.+69.9*pow(q,1.08));
        c2=14.2/ac + 386./(1.+69.9*pow(q,1.08));

        // Equation 17 - CDM transfer function
        tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) +
            (1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q);

        // Equation 15
        y=(1.+ze)/(1.+zd);
        g=y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)));

        // Equation 14
        ab=g*2.07*rke*s/pow(1.+rd,0.75);

        // Equation 23
        bn=8.41*pow(om_mhsq,0.435);

        // Equation 22
        ss=s/pow(1.+pow(bn/rk/s,3.),1./3.);

        // Equation 24ggkkkkkkkkkkkkkkk
        bb=0.5+(om_b/om_m) + (3.-2.*om_b/om_m)*sqrt(pow(17.2*om_mhsq,2.)+1.);

        // Equations 19 & 21
        tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q*q)/(1+pow(rk*s/5.2,2.));
        tb=(tb+ab*exp(-pow(rk/rks,1.4))/(1.+pow(bb/rk/s,3.)))*sin(rk*ss)/rk/ss;

        // Equation 8
        tk_eh=(om_b/om_m)*tb+(1.-om_b/om_m)*tc;

