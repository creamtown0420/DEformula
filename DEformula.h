// change <mp::cpp_dec_float<16>> to any precision
namespace mp = boost::multiprecision;
using Real16 = mp::number<mp::cpp_dec_float<16>>;
using namespace std;
int kmin = 3;
int kmax = 14;
Real16 hr = 6.0f;
Real16 c0 = 0.01f; // safe value for eps
Real16 pi = boost::math::constants::pi<Real16>();

Real16 pi2 = pi / 2.0f;
template <class T>
Real16 dde3(T integrad, Real16 xa, Real16 xb, Real16 ya, Real16 yb, Real16 za, Real16 zb, Real16 eps, Real16 s, int info)
{ 
    int l, nc, j;
    Real16 h, s0, xt, wt, t, ft, as, shk, mba, pba, err, seps;
    info = 0;
    seps = 0.001f;
    h = hr * 0.5f;
    mba = (xb - xa) * 0.5f;
    pba = (xb + xa) * 0.5f; 
    ft = dde3dsyz(integrad, ya, yb, za, zb, pba, eps, ft, info);

    s0 = ft * h * pi2 * mba;
    cout << "s0" << s0 << "\n";
    nc = 1;
    for (l = 2; l <= kmax; l++)
    {
        s = 0.0f;
        nc = 2 * nc;
        h = h * 0.5f;
        for (j = 1; j <= nc; j++)
        {
            Real16 f50 = 2 * j - nc - 1;
            t = (f50)*h;

            shk = pi2 * mp::sinh(t);
            xt = mp::tanh(shk);

            wt = pi2 * mp::cosh(t) / (mp::cosh(shk) * mp::cosh(shk));

            ft = dde3dsyz(integrad, ya, yb, za, zb, mba * xt + pba, eps, ft, info);

            s = s + ft * wt;
        }
        s = s0 * 0.5f + s * h * mba;

        as = mp::abs(s);
        err = mp::abs(s - s0);
        if (as >= 1)
            err = err / as;

        if (err <= seps & l >= kmin)
            break;
        s0 = s;
    }
    if (l == kmax + 1)
        info += 1;
    return s;
}
template <class T>
Real16 dde3dsyz(T integrad, Real16 ya, Real16 yb, Real16 za, Real16 zb, Real16 xc, Real16 eps, Real16 s, int info)
{
    int j, nc, l;
    Real16 h, s0, xt, wt, t, ft, as, shk, mba, pba, err, seps;

    seps = c0 * sqrt(eps);
    h = hr * 0.5f;
    ft = 0.0f;
    s = 0.0f;
    mba = (yb - ya) * 0.5f;
    pba = (yb + ya) * 0.5f;
    ft = dde3dsz(integrad, za, zb, xc, pba, eps, ft, info);
    s0 = ft * h * pi2 * mba;
    nc = 1;
    for (l = 2; l <= kmax; l++)
    {
        s = 0.0f;
        nc = 2 * nc;
        h = h * 0.5f;
        for (j = 1; j <= nc; j++)
        {
            Real16 f50 = 2 * j - nc - 1;
            t = (f50)*h;
            shk = pi2 * mp::sinh(t);

            xt = mp::tanh(shk);
            wt = pi2 * mp::cosh(t) / (mp::cosh(shk) * mp::cosh(shk));
            ft = dde3dsz(integrad, za, zb, xc, mba * xt + pba, eps, ft, info);

            s = s + ft * wt;
        }
        s = s0 * 0.5f + s * h * mba;
        as = mp::abs(s);
        err = mp::abs(s - s0);
        if (as >= 1)
            err = err / as;
        if (err <= seps & l >= kmin)
            break;
        s0 = s;
    }
    if (l == kmax + 1)
        info += 1;

    return s;
}
template <class T>
Real16 dde3dsz(T integrad, Real16 za, Real16 zb, Real16 xc, Real16 yc, Real16 eps, Real16 s, int info)
{
    int j, nc, l;
    Real16 h, s0, xt, wt, t, as, shk, mba, pba, err, seps;
    seps = c0 * sqrt(eps);
    h = hr * 0.5f;

    mba = (zb - za) * 0.5f;
    pba = (zb + za) * 0.5f;
    s0 = integrad(xc, yc, pba) * h * pi2 * mba;

    s = 0.0f;
    nc = 1;

    for (l = 2; l <= kmax; l++)
    {
        s = 0.0f;
        nc = 2 * nc;
        h = h * 0.5f;
        for (j = 1; j <= nc; j++)
        {
            Real16 f50 = 2 * j - nc - 1;
            t = (f50)*h;
            shk = pi2 * mp::sinh(t);

            xt = mp::tanh(shk);
            wt = pi2 * mp::cosh(t) / (mp::cosh(shk) * mp::cosh(shk));
            s += integrad(xc, yc, mba * xt + pba) * wt;
        }
        s = s0 * 0.5f + s * h * mba;
        as = mp::abs(s);
        err = mp::abs(s - s0);
        if (as >= 1)
            err = err / as;
        if (err <= seps & l >= kmin)
            break;
        s0 = s;
    }
    if (l == kmax + 1)
        info += 1;

    return s;
}
