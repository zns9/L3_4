float rk4(float xin, float xold, float tau, float dt)
{
    float k1, k2, k3, k4, xout;

    k1 = dt*(xin- xold        )/tau;
    k2 = dt*(xin-(xold+k1/2.0))/tau;
    k3 = dt*(xin-(xold+k2/2.0))/tau;
    k4 = dt*(xin-(xold+k3    ))/tau;

    xout = xold + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    return(xout);
}


float rk2(float xin, float xold, float tau, float dt)
{
    float k1, k2, xout;
    k1 = dt*(xin- xold        )/tau;
    k2 = dt*(xin-(xold+k1/2.0))/tau;
    xout = xold + (k1 + k2)/2.0;
    return(xout);
}


float euler(float xin, float xold, float tau, float dt)
{
    float k1, xout;
    k1 = dt*(xin- xold)/tau;
    xout = xold + k1;
    return(xout);
}

float euler_bk(float xin, float xold, float tau, float dt)
{
    float k1, k2, xout;
    k1 = xold + dt*xin/tau;
    k2 = 1.0 + (dt/tau);
    xout = k1 / k2;
    return(xout);
}


float euler_tr(float xin, float xold, float tau, float dt)
{
    float k1, k2, xout;
    k1 = xold*(1.0 - dt/(2.0*tau)) + dt*xin/tau;
    k2 = 1.0 + dt/(2.0*tau);
    xout = k1 / k2;
    return(xout);
}


float rk4_2ndo(float xin, float *work, float a, float b, float dt)
{
    float kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4, az, x, y, xout;
    
    az = a/2.0;
    x  = work[0];
    y  = work[1];

    kx1 = dt*(xin - az*x - b*y        );
    ky1 = dt * x;
    kx2 = dt*(xin - az*(x + kx1/2.0) - b*(y + ky1/2.0));
    ky2 = dt*(x + kx1/2.0);
    kx3 = dt*(xin - az*(x + kx2/2.0) - b*(y + ky2/2.0));
    ky3 = dt*(x + kx2/2.0);
    kx4 = dt*(xin - az*(x + kx3    ) - b*(y + ky3    ));
    ky4 = dt*(x + kx3    );

    x = x + (kx1 + 2.0*kx2 + 2.0*kx3 + kx4)/6.0;
    y = y + (ky1 + 2.0*ky2 + 2.0*ky3 + ky4)/6.0;
    work[0] = x;
    work[1] = y;
    return(y);
}


float diff(float xin, float xinold, float dt)
{
    float xout;
    xout = (xin - xinold)/dt;
    return(xout);
}


float integ(float xin, float xinold, float xold, float dt)
{
    float xout;
    xout = xold + (xin + xinold)/2.0 * dt;
    return(xout);
}

unsigned int rotr(unsigned int intval, int bit)
{
    if(bit>=2 && bit<=16)
    {
        intval = ((intval << (bit-1)) | (intval>>1));
        intval = intval & 0xffff >> (16-bit);
        return(intval);
    }
    else if (bit==1)
    {
        return(intval);
    }
    else
    {
        return(-1);
    }
}

unsigned int rotl(unsigned int intval, int bit)
{
    if(bit>=2 && bit<=16)
    {
        intval = ((intval >> (bit-1)) | (intval<<1));
        intval = intval & 0xffff >> (16-bit);
        return(intval);
    }
    else if (bit==1)
    {
        return(intval);
    }
    else
    {
        return(-1);
    }
}


float limit(float data, float lmax, float lmin)
{
    if ( data > lmax) data = lmax;
    if ( data < lmin) data = lmin;
    return(data);
}


float pid(float xin, float kp, float ki,float kd, float *work, float dt)
{
    float a, b, c, xout;

    a = xin;
    b = work[1] + (xin + work[0])/2.0 * dt;
    c = (xin - work[0])/dt;
    work[0] = xin;
    work[1] = b;

    xout = a*kp + b*ki + c*kd;
    return(xout);
}

float pi_d(float xin, float fb, float kp, float ki,float kd, float *work, float dt)
{
    float a, b, c, xout;

    a = xin;
    b = work[1] + (xin + work[0])/2.0 * dt;
    c = (fb - work[2])/dt;
    work[0] = xin;
    work[1] = b;
    work[2] = fb;

    xout = a*kp + b*ki - c*kd;
    return(xout);
}

float i_pd(float xin, float fb, float kp, float ki,float kd, float *work, float dt)
{
    float a, b, c, xout;

    a = fb;
    b = work[1] + (xin + work[0])/2.0 * dt;
    c = (fb - work[2])/dt;
    work[0] = xin;
    work[1] = b;
    work[2] = fb;

    xout = -a*kp + b*ki - c*kd;
    return(xout);
}



float leadlag(float xin, float *work, float a0, float a1, float b0, float b1, float dt)
{
    float k1, k2, k3, k4, c0, c1, xold, xout;
    
    xold = work[0];
    c0 = a0 - a1*b0/b1;
    c1 = b0/b1;

    k1 = dt*(c0*xin - c1*xold);
    k2 = dt*(c0*xin - c1*(xold+k1/2.0));
    k3 = dt*(c0*xin - c1*(xold+k2/2.0));
    k4 = dt*(c0*xin - c1*(xold+k3    ));

    xold = xold + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    work[0] = xold;
    xout = (xold + xin*a1)/b1;
    return(xout);
}

float ad_x_ad_[100]={0.0};
int ad_x_old_ad_[100]={0,0};

unsigned int adc(int adc_ch, int bits, float x, float xmin, float xmax, float ct, float dt)
{
    float xdiv;
    int   adc1, adresolution;

    adresolution = pow(2,bits);
    xdiv  = xmax - xmin;
    
    if (x > xmax) x = xmax;
    if (x < xmin) x = xmin;

    ad_x_ad_[adc_ch] += dt;
    if (ad_x_ad_[adc_ch] >= ct)
    {
        adc1    = (x - xmin) / xdiv * adresolution;
        if (adc1>(adresolution-1)) adc1 = adresolution-1;
        if (adc1<0)                adc1 = 0; 
        ad_x_ad_[adc_ch] = 0;
    }
    else
    {
        adc1 = ad_x_old_ad_[adc_ch];
    }
    ad_x_old_ad_[adc_ch] = adc1;
    return(adc1);
}

float dac(int bits, int x, float xmin, float xmax)
{
    float xdiv, dac1;
    int   daresolution;

    daresolution = pow(2,bits);
    xdiv  = xmax - xmin;
    
    if (x > daresolution-1) x = daresolution-1;
    if (x < 0) x = 0;

    dac1    = (xdiv * x)/daresolution + xmin;
    return(dac1);
}

